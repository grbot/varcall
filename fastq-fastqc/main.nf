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
        }.into{samples_1; samples_2}

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file)  from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\n"
    """
}
*/

process log_tool_version_fastqc {
    tag { "${params.project_name}.ltVF" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'fastqc'

    output:
    file("tool.fastqc.version") into tool_version_fastq

    script:
    """
    fastqc --version > tool.fastqc.version
    """
}

process log_tool_version_multiqc {
    tag { "${params.project_name}.ltVM" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.multiqc.version") into tool_version_multiqc

    script:
    """
    multiqc --version > tool.multiqc.version
    """
}

process run_fastqc {
    tag { "${params.project_name}.${sample_id}.rFqc" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    cpus { 2 } 
    label 'fastqc'

    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_2

    output:
    file("${sample}_fastqc/*.zip") into fastqc_files
    
    script:
    
    """
    mkdir ${sample_id}_fastqc
    fastqc --outdir ${sample_id}_fastqc \
    -t ${task.cpus} \
    ${fastq_r1_file} \
    ${fastq_r2_file}   
    """
}

process run_multiqc {
    tag { "${params.project_name}.${sample_id}.rMqc" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'multiqc'

    input:
    file('*') from fastqc_files.collect()

    output:
    file('multiqc_report.html') 

    script:
    """
    multiqc .
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
