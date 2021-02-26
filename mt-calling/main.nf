#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM'])
            return [ sample_id, bam_file ]
        }.into{samples; samples_1}

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_dict = Channel.fromPath(params.ref_seq_dict).toList()
dbsnp = Channel.fromPath(params.dbsnp).toList()
dbsnp_index = Channel.fromPath(params.dbsnp_index).toList()

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(bam_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tBAM: ${bam_file}\n"
    """
}
*/

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVG" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process log_tool_version_bcftools {
    tag { "${params.project_name}.ltViB" }
    echo true
    publishDir "${params.out_dir}/", mode: 'move', overwrite: false
    label 'bcftools'

    output:
    file("tool.bcftools.version") into tool_version_bcftools

    script:
    mem = task.memory.toGiga() - 3
    """
    bcftools --version > tool.bcftools.version 2>&1
    """
}

process subset_bam_to_mt {
    tag { "${params.project_name}.${sample_id}.sBM" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'gatk'
    input:
    set val(sample_id), file(bam_file) from samples
    file (ref_seq)
    file (ref_seq_index)
    file (ref_seq_dict)

    output:
    set val(sample_id), file("${sample_id}.MT.bam") into mt_bam_file

    script:
    """
      gatk PrintReads \
      -R ${ref_seq} \
      -L MT \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      -I ${bam_file} \
      -O ${sample_id}.MT.bam
    """
}

process call_mt {
    tag { "${params.project_name}.${sample_id}.cM" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'gatk'
    input:
    set val(sample_id), file (bam_file) from mt_bam_file
    file (ref_seq)
    file (ref_seq_index)
    file (ref_seq_dict)

    output:
    set val(sample_id), file("${sample_id}.MT.vcf.gz"), file("${sample_id}.MT.vcf.gz.tbi"), file("${sample_id}.MT.vcf.gz.stats") into mt_call_file

    script:
    mem = task.memory.toGiga() - 4
    """
      gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" Mutect2 \
      -R ${ref_seq} \
      -I ${bam_file} \
      -O ${sample_id}.MT.vcf.gz \
      --annotation StrandBiasBySample \
      --mitochondria-mode \
      --max-reads-per-alignment-start 75 \
      --max-mnp-distance 0
    """
}


process filter_mt_calls {
    tag { "${params.project_name}.${sample_id}.fMC" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'gatk'
    input:
    set val(sample_id), file (call_file), file (call_file_index), file (call_file_stats) from mt_call_file
    file (ref_seq)
    file (ref_seq_index)
    file (ref_seq_dict)

    output:
    set val(sample_id), file("${sample_id}.MT.filtered.vcf.gz"), file("${sample_id}.MT.filtered.vcf.gz.tbi") into mt_filtered_call_file

    script:
    mem = task.memory.toGiga() - 4
    """
      gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" FilterMutectCalls \
      -R ${ref_seq} \
      -V ${call_file} \
      -O ${sample_id}.MT.filtered.vcf.gz \
      --stats ${call_file_stats} \
      --mitochondria-mode \
    """
}

process select_pass {
    tag { "${params.project_name}.${sample_id}.sP" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bcftools'
    input:
      set val(sample_id), file (call_file), file (call_file_index) from mt_filtered_call_file
    output:
      set val(sample_id), file("${sample_id}.MT.filtered-pass.vcf.gz"), file("${sample_id}.MT.filtered-pass.vcf.gz.tbi") into mt_filtered_call_file_pass

    script:
    """
    bcftools view \
    --include "FILTER='PASS'" \
    -O z \
    -o "${sample_id}.MT.filtered-pass.vcf.gz" \
    ${call_file} 
    bcftools index -t "${sample_id}.MT.filtered-pass.vcf.gz"
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
