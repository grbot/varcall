#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fastq_r1_file = file(row['FastqR1'])
            def fastq_r2_file = file(row['FastqR2'])
            return [ sample_id, fastq_r1_file, fastq_r2_file ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_amb = Channel.fromPath(params.ref_seq_amb).toList()
ref_seq_ann = Channel.fromPath(params.ref_seq_ann).toList()
ref_seq_bwt = Channel.fromPath(params.ref_seq_bwt).toList()
ref_seq_pac = Channel.fromPath(params.ref_seq_pac).toList()
ref_seq_sa = Channel.fromPath(params.ref_seq_sa).toList()
known_indels_1 = Channel.fromPath(params.known_indels_1).toList()
known_indels_1_index = Channel.fromPath(params.known_indels_1_index).toList()
known_indels_2 = Channel.fromPath(params.known_indels_2).toList()
known_indels_2_index = Channel.fromPath(params.known_indels_2_index).toList()
dbsnp = Channel.fromPath(params.dbsnp).toList()
dbsnp_index = Channel.fromPath(params.dbsnp_index).toList()

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\n"
    """
}

process run_bwa {
    tag { "${params.project_name}.${sample_id}.rBwa" }
    memory { 64.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_2
    file (ref) from ref_seq
    file (ref_index) from ref_seq_index
    file (ref_amb) from ref_seq_amb
    file (ref_ann) from ref_seq_ann
    file (ref_bwt) from ref_seq_bwt
    file (ref_pac) from ref_seq_pac
    file (ref_sa) from ref_seq_sa

    output:
    set val(sample_id), file("${sample_id}.bam")  into raw_bam

    script:
    nr_threads = task.cpus - 1
    readgroup_info="@RG\\tID:$sample_id.0\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"
    """
    bwa mem \
    -R \"${readgroup_info}\" \
    -t ${nr_threads}  \
    -M \
    ${ref} \
    ${fastq_r1_file} \
    ${fastq_r2_file} | \
    samtools sort \
    -@ ${nr_threads} \
    -m ${params.memory_per_thread} \
    - > ${sample_id}.bam
    """
}

process run_mark_duplicates {
    tag { "${params.project_name}.${sample_id}.rMD" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file) from raw_bam

    output:
    set val(sample_id), file("${sample_id}.md.bam"), file("${sample_id}.md.bai")  into md_bam

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms4g -Xmx${mem}g"  \
    MarkDuplicates \
    --MAX_RECORDS_IN_RAM 5000 \
    --INPUT ${bam_file} \
    --METRICS_FILE ${bam_file}.metrics \
    --TMP_DIR . \
    --ASSUME_SORT_ORDER coordinate \
    --CREATE_INDEX true \
    --OUTPUT ${sample_id}.md.bam
    """
}

process run_create_recalibration_table {
    tag { "${params.project_name}.${sample_id}.rCRT" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file), file(bam_file_index) from md_bam
    file (dbsnp_file) from dbsnp
    file (dbsnp_index_file) from dbsnp_index
    file (known_indels_1_file) from known_indels_1
    file (known_indels_1_index_file) from known_indels_1_index
    file (known_indels_2_file) from known_indels_2
    file (known_indels_2_index_file) from known_indels_2_index

    output:
    set val(sample_id), file("${sample_id}.md.bam"), file("${sample_id}.md.bai"), file("${sample_id}.recal.table")  into recal_table

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    BaseRecalibrator \
    --input ${bam_file} \
    --output ${sample_id}.recal.table \
    --tmp-dir . \
    -R ${params.ref_seq} \
    --known-sites ${dbsnp_file} \
    --known-sites ${known_indels_1_file} \
    --known-sites ${known_indels_2_file}
    """
}

process run_recalibrate_bam {
    tag { "${params.project_name}.${sample_id}.rRB" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file), file(bam_file_index), file(recal_table_file) from recal_table

    output:
    set val(sample_id), file("${sample_id}.md.recal.bam")  into recal_bam
    set val(sample_id), file("${sample_id}.md.recal.bai")  into recal_bam_index

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
     ApplyBQSR \
    --input ${bam_file} \
    --output ${sample_id}.md.recal.bam \
    --tmp-dir . \
    -R ${params.ref_seq} \
    --create-output-bam-index true \
    --bqsr-recal-file ${recal_table_file}
    """
}

process run_samtools_flagstat {
    tag { "${params.project_name}.${sample_id}.rSF" }
    memory { 4.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'move', overwrite: false
    label 'bwa_samtools'

    input:
    set val(sample_id), file(bam_file) from recal_bam

    output:
    set val(sample_id), file("${sample_id}.md.recal.flagstat")  into recal_stats

    """
    samtools flagstat  \
    --threads ${params.bwa_threads} \
    ${bam_file} > ${sample_id}.md.recal.flagstat  \
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
