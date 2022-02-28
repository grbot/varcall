#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM'])
            return [ sample_id, bam_file ]
        }.into{samples_1; samples_2}

process log_tool_version_fastqc {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/bam-qc", mode: 'copy', overwrite: false
    label 'fastqc'

    output:
    file("tool.fastqc.version") into tool_version_fastqc

    script:
    """
    fastqc -v > tool.fastqc.version
    """
}

process log_tool_version_multiqc {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/bam-qc", mode: 'copy', overwrite: false
    label 'multiqc'

    output:
    file("tool.multiqc.version") into tool_version_multiqc

    script:
    """
    multiqc --version > tool.multiqc.version
    """
}

process log_tool_version_samtools {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/bam-qc", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.samtools.version") into tool_version_samtool

    script:
    """
    samtools --version > tool.samtools.version
    """
}

process run_flagstat {
    tag { "${params.project_name}.${sample_id}.rF" }
    echo true
    publishDir "${params.out_dir}/bam-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file) from samples_1

    output:
       file("${bam_file}.flagstat") into flagstat_file

    script:
    """
    samtools flagstat \
    -@ 1 \
    ${bam_file} > ${bam_file}.flagstat
    """
}


process runFastQC{
    tag { "${params.project_name}.${sample_id}.rFQC" }
    publishDir "${params.out_dir}/bam-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'fastqc'

    input:
        set val(sample_id), file(bam_file) from samples_2

    output:
        file("${sample_id}_fastqc/*.zip") into fastqc_file

    """
    mkdir ${sample_id}_fastqc
    fastqc --outdir ${sample_id}_fastqc \
    ${bam_file}
    """
}

combined = flagstat_file.mix(fastqc_file)

process runMultiQC{
    tag { "${params.project_name}.rMQC" }
    publishDir "${params.out_dir}/bam-qc/", mode: 'move', overwrite: false
    label 'multiqc'

    input:
        file('*') from combined.collect()

    output:
        file('multiqc_report.html')

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
