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


process cram_to_fastq {
    tag { "${params.project_name}.${sample_id}.ctF" }
    echo true
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    input:
    set val(sample_id), val(bam_file) from samples

    output:
    set val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") into fastq

    script:
    """
    ${params.samtools_base}/samtools collate -O --reference ${params.ref_seq} ${bam_file} tmp.collate | ${params.samtools_base}/samtools fastq -1 ${sample_id}_R1.fastq.gz -2 ${sample_id}_R2.fastq.gz -0 /dev/null -s /dev/null -N -F 0x900 --reference ${params.ref_seq} -
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
