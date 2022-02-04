#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fastq_dir = file(row['FastqDir'])
            return [ sample_id, fastq_dir ]
        }.into{samples_1; samples_2}

process run_combine_fastq_r1 {
    tag { "${params.project_name}.${sample_id}.rCfr1" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'combine_fastq'
  
    input:
    set val(sample_id), file(fastq_dir) from samples_1
 
    output:
    set val("$sample_id"), file("${sample_id}_1.fq.gz")  into fastq_r1
      
    script:
    // Barcoded fastqs are sometimes distributed as *_1.fastq.gz or *_1_${sample_id}.fastq.gz or *_1_${sample_id}R.fastq.gz . 
    """
    shopt -s extglob
    w=`pwd`
    cd ${fastq_dir}
    zcat @(*_1.fq.gz|*_1_${sample_id}.fq.gz|*_1_${sample_id}R.fq.gz) | gzip - > \$w/${sample_id}_1.fq.gz
    cd \$w
    """
}

process run_combine_fastq_r2 {
    tag { "${params.project_name}.${sample_id}.rCfr2" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'combine_fastq'
 
    input:
    set val(sample_id), file(fastq_dir) from samples_2
 
    output:
    set val("$sample_id"), file("${sample_id}_2.fq.gz")  into fastq_r2
      
    script:
    // Barcoded fastqs are sometimes distributed as *_2.fastq.gz or *_2_${sample_id}.fastq.gz or *_2_${sample_id}R.fastq.gz. 
    """
    shopt -s extglob
    w=`pwd`
    cd ${fastq_dir}
    zcat @(*_2.fq.gz|*_2_${sample_id}.fq.gz|*_2_${sample_id}R.fq.gz) | gzip - > \$w/${sample_id}_2.fq.gz
    cd \$w
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
