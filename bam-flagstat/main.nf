#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM'])
            return [ sample_id, bam_file ]
        }.set{samples}


process run_flagstat {
    tag { "${params.project_name}.${sample_id}.rF" }
    echo true
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    input:
    set val(sample_id), file(bam_file) from samples

    output:
    set val(sample_id), file("${bam_file}.flagstat") into cram_file

    script:
    """
    ${params.samtools_base}/samtools flagstat \
    -@ 1 \
    ${bam_file} > ${bam_file}.flagstat  \
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
