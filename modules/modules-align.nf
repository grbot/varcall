#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref            = file(params.ref, type: 'file')
known_indels_1 = file(params.known_indels_1, type: 'file')
known_indels_2 = file(params.known_indels_2, type: 'file')
dbsnp          = file(params.dbsnp, type: 'file')
project_name   = params.project_name
outdir         = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_bwa {
    tag { "${sample_id}.rBwa" }
    label 'bwa_samtools'
    memory { 64.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    // publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: raw_bam
      
    script:
    // readgroup_info="@RG\\tID:${flowcell}.${lane}\\tLB:LIBA\\tSM:${sample_id}\\tPL:Illumina"
      
    // if (lane == ".") {
    //     sample_id = "${sample_id}"
    // } else {
  	//     sample_id = "${sample_id}-${flowcell}.${lane}"
    // }
    
    nr_threads = task.cpus - 1

    """
    flowcell=`zcat ${fastq_r1} | head -n 1 | awk -F':' '{ print \$3 }'`
    lane=`zcat ${fastq_r1} | head -n 1 | awk -F':' '{ print \$4 }'`
    readgroup_info="@RG\\tID:\$flowcell.\$lane\\tLB:LIBA\\tSM:${sample_id}\\tPL:Illumina"
    bwa mem \
        -R \"\$readgroup_info\" \
        -t ${nr_threads}  \
        -K 100000000 \
        -Y \
        ${ref} \
        ${fastq_r1} \
        ${fastq_r2} | \
        samtools sort \
        -@ ${nr_threads} \
        -m ${params.memory_per_thread} \
        - > ${sample_id}.bam
    """
}

process run_mark_duplicates {
    tag { "${sample_id}.rMD" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    // publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: md_bam

    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms4g -Xmx${mem}g" \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 5000 \
        --INPUT ${bam} \
        --METRICS_FILE ${bam}.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${sample_id}.md.bam
    """
}

process run_create_recalibration_table {
    tag { "${sample_id}.rCRT" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    // publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), path("${sample_id}.recal.table"), emit: recal_table

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        BaseRecalibrator \
        --input ${bam} \
        --output ${sample_id}.recal.table \
        --tmp-dir . \
        -R ${ref} \
        --known-sites ${dbsnp} \
        --known-sites ${known_indels_1} \
        --known-sites ${known_indels_2}
    """
}

process run_recalibrate_bam {
    tag { "${sample_id}.rRB" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    cpus { 2 }
    // publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam), path(index), file(recal_table)

    output:
    tuple val(sample_id), path("${sample_id}.md.recal.bam"), path("${sample_id}.md.recal.bai"), emit: recal_bam

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        ApplyBQSR \
        --input ${bam} \
        --output ${sample_id}.md.recal.bam \
        --tmp-dir . \
        -R ${ref} \
        --create-output-bam-index true \
        --bqsr-recal-file ${recal_table}
    """
}

process bam_to_cram {
    tag { "${sample_id}.btC" }
    label 'bwa_samtools'
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path("${bam.baseName}.cram"), path("${bam.baseName}.crai"), emit: cram_file

    """
    samtools view --reference ${ref} \
        --output-fmt cram,version=3.0 -o ${bam.baseName}.cram ${bam}
    samtools index ${bam.baseName}.cram ${bam.baseName}.crai
    """
}

process run_cram_flagstat {
    tag { "${sample_id}.rCF" }
    label 'bwa_samtools'
    memory { 4.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(cram), path(index)

    output:
    tuple val(sample_id), path(cram), path(index), path("${cram}.flagstat"), emit: cram_stats

    """
    samtools flagstat \
        --threads ${params.bwa_threads} \
        ${cram} > ${cram}.flagstat  \
    """
}

process create_cram_md5sum {
    tag { "${sample_id}.cCMD5" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/align/${project_name}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(cram), path(index), path(flagstat)

    output:
    tuple val(sample_id), path("${cram}.md5"), path("${index}.md5"), path("${flagstat}.md5"), emit: cram_all_md5sum

    """
    md5sum ${cram} > ${cram}.md5
    md5sum ${index} > ${index}.md5
    md5sum ${flagstat} > ${flagstat}.md5
    """
}
