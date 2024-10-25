#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM/CRAM'])
            return [ sample_id, bam_file ]
        }.set{samples}

ref_seq = Channel.fromPath(params.ref_seq).toList()

process collate {
    tag { "${params.project_name}.${sample_id}.C" }
    memory { 32.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'samtools'
    input:
    set val(sample_id), file(bam_file) from samples
    file (ref) from ref_seq

    output:
    set val(sample_id), file("${sample_id}.collate.bam") into collate

    script:
    """
    samtools collate --reference ${ref} -O ${bam_file}  ${sample_id}.collate
    """
}

process cram_to_fastq {
    tag { "${params.project_name}.${sample_id}.ctF" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'move', overwrite: false
    label 'samtools'
    input:
    set val(sample_id), file(bam_file) from collate
    file (ref) from ref_seq
    
    output:
    set val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") into fastq
    
    script:
    """
    samtools fastq -1 ${sample_id}_R1.fastq.gz -2 ${sample_id}_R2.fastq.gz -0 /dev/null -s /dev/null -N -F 0x900 --reference ${ref} ${bam_file}
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
