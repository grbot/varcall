#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process print_sample_info {
    tag { sample_id }
    echo true

    input:
    tuple val(sample_id), val(gender), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane), path(bam), path(gvcf)

    script:
    """
    printf "[Sample Info]"
    printf "Sample   : ${sample_id}"
    printf "Fastq R1 : ${fastq_r1}"
    printf "Fastq R2 : ${fastq_r2}"
    printf "Flowcell : ${flowcell}"
    printf "Lane     : ${lane}\n"
    printf "BAM      : ${lane}\n"
    printf "GVCF     : ${lane}\n"
    printf "--------------------------------------------------
    """
}

process log_tool_version_samtools_bwa {
    tag { "tool_ver" }
    label 'bwa_samtools'
    echo true
    publishDir "${outdir}/", mode: 'copy', overwrite: true

    output:
    path("tool.samtools.bwa.version"), emit: tool_version_samtools_bwa

    script:
    """
    samtools --version > tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

process log_tool_version_gatk {
    tag { " tool_ver" }
    label 'gatk'
    echo true
    publishDir "${outdir}/", mode: 'copy', overwrite: true
    
    output:
    path("tool.gatk.version"), emit: tool_version_gatk
    
    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

// process generate_sample_sheet {
//     tag { "${params.project_name}.${sample_id}.btC" }
//     label 'bwa_samtools'
//     memory { 4.GB * task.attempt }
//     publishDir "${outdir}/${sample_id}", mode: 'copy', overwrite: true

//     input:
//     tuple val(sample_id), val(gender), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane), path(bam), path(gvcf), path(bam), path(index)

//     output:
//     tuple val(sample_id), path("${bam.baseName}.tsv"), emit: sample_sheet_bam

//     """
    
//     """
// }
