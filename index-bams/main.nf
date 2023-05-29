#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM/CRAM'])
            return [ sample_id, bam_file ]
        }.set{samples}

process log_tool_version_samtools_bwa {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/bam-index", mode: 'move', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.samtools.bwa.version") into tool_version_samtools_bwa

    script:
    """
    samtools --version > tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

process index_bam {
    tag { "${sample_id}" }
    publishDir "${params.out_dir}/bam-index/${sample_id}", mode: 'move', overwrite: false
    label 'bwa_tools'
    input:
    set val(sample_id), file(bam_file) from samples
    output:
    file("*.bai")
    script:
    """
    samtools index ${bam_file}
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
