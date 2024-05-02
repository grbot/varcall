#!/usr/bin/env nextflow
nextflow.enable.dsl=2

outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process print_sample_info {
    tag { sample_id }
    debug true

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
    debug true
    publishDir "${outdir}/${params.worflow}", mode: 'copy', overwrite: true

    output:
    path("tool.samtools.bwa.version"), emit: tool_version_samtools_bwa

    script:
    """
    samtools --version > tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

process log_tool_version_gatk {
    tag { "tool_ver" }
    label 'gatk'
    debug true
    publishDir "${outdir}/${params.workflow}", mode: 'copy', overwrite: true
    
    output:
    path("tool.gatk.version"), emit: tool_version_gatk
    
    script:
    mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

