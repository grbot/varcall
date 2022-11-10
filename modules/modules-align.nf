#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref = file(params.ref, type: 'file')
known_indels_1 = file(params.known_indels_1, type: 'file')
known_indels_2 = file(params.known_indels_2, type: 'file')
dbsnp = file(params.dbsnp, type: 'file')
outdir = file(params.outdir, type: 'dir')
outdir.mkdir()

process print_sample_info {
    tag { sample_id }
    echo true

    input:
    tuple val(sample_id), val(gender), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane), path(bam), path(gvcf)

    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1}\tFastqR2: ${fastq_r2}\tFlowcell: ${flowcell}\tLane: ${lane}\n"
    """
}

process log_tool_version_samtools_bwa {
    tag { "tool_ver" }
    label 'bwa_samtools'
    echo true
    publishDir "${outdir}/", mode: 'copy', overwrite: false

    output:
    path("tool.samtools.bwa.version"), emit: tool_version_samtools_bwa

    script:
    """
    samtools --version > tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

process log_tool_version_gatk {
    tag { " tool_ver" }
    label 'gatk'
    echo true
    publishDir "${outdir}/", mode: 'copy', overwrite: false
    
    output:
    path("tool.gatk.version"), emit: tool_version_gatk
    
    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process run_bwa {
    tag { "${params.project_name}.${sample_id}.rBwa" }
    label 'bwa_samtools'
    memory { 64.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), val(gender), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane), path(bam), path(gvcf)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: raw_bam
      
    script:
    readgroup_info="@RG\\tID:${flowcell}.${lane}\\tLB:LIBA\\tSM:${sample_id}\\tPL:Illumina"
      
    if (lane == "0") {
        sample_id = "${sample_id}"
    } else {
  	    sample_id = "${sample_id}-${flowcell}.${lane}"
    }
    
    nr_threads = task.cpus - 1

    """
    bwa mem \
        -R \"${readgroup_info}\" \
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
    tag { "${params.project_name}.${sample_id}.rMD" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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
    tag { "${params.project_name}.${sample_id}.rCRT" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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
    tag { "${params.project_name}.${sample_id}.rRB" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    cpus { 2 }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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
    tag { "${params.project_name}.${sample_id}.btC" }
    label 'bwa_samtools'
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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
    tag { "${params.project_name}.${sample_id}.rCF" }
    label 'bwa_samtools'
    memory { 4.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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
    tag { "${params.project_name}.${sample_id}.cCMD5" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: false

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

workflow.onComplete {
    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
