#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_dir = file(row['BAMdir'])
            return [ sample_id, bam_dir ]
        }.set{samples}

process log_tool_version {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.version") into tool_version

    script:
    """
    samtools --version > tool.version
    """
}

process merge_bams {
    tag { "${params.project_name}.${sample_id}.mB" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_dir) from samples

    output:
    set val(sample_id), file("${bam_dir.baseName}.bam") into bam_file_1

    script:
    """
    samtools merge \
    ${bam_dir.baseName}.bam \
    ${bam_dir}/*
    """
}

process index_bam {
    tag { "${params.project_name}.${sample_id}.iB" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'move', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file) from bam_file_1

    output:
    set val(sample_id), file("${sample_id}.bai") into bam_index


    script:
    """
    samtools index  \
    ${bam_file} ${sample_id}.bai
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
