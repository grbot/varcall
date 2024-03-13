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

/* 
// Unable to fine trimmomatic version, need to look into it another time
process log_tool_version_trimmomatic {
    tag { "${params.project_name}.ltVT" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'trimmomatic'

    output:
    file("tool.trimmomatic.version") into tool_version_trimmomatic

    script:
    """
    trimmomatic --version >tool_version_trimmomatic 2>&1
    """
}
*/

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

process run_trimmomatic {
    tag { "${params.project_name}.${sample_id}.rT" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    cpus { "${params.trimmomatic_threads}" }
    label 'trimmomatic'

    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_2

    output:
    file("*.gz") into fastq_files
    file("${sample_id}_trim_out.log") into trimmomatic_trim_files
    
    script:
    """
    trimmomatic \
    PE -phred33 \
    -threads ${params.trimmomatic_threads} \
    ${fastq_r1_file} \
    ${fastq_r2_file} \
    ${sample_id}_R1.gz \
    ${sample_id}_R1_unpaired.gz \
    ${sample_id}_R2.gz \
    ${sample_id}_R2_unpaired.gz \
    CROP:${params.crop} 2> ${sample_id}_trim_out.log
    """
}

process run_multiqc {
    tag { "${params.project_name}.rMqc" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'multiqc'

    input:
    file('*') from trimmomatic_trim_files.collect()

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
