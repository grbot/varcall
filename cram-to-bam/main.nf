#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def cram_file = file(row['CRAM'])
            return [ sample_id, cram_file ]
        }.set{samples}

ref_seq = Channel.fromPath(params.ref_seq).toList()

process log_tool_version {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.version") into tool_version

    script:
    """
    samtools --version > tool.version
    """
}

process cram_to_bam {
    tag { "${params.project_name}.${sample_id}.ctB" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(cram_file) from samples
    file (ref) from ref_seq

    output:
    set val(sample_id), file("${cram_file.baseName}.bam") into bam_file

    script:
    """
    samtools view \
    -b \
    --reference ${ref} \
    -o ${cram_file.baseName}.bam  ${cram_file}
    """
}

bam_file.into{ bam_file_1; bam_file_2; bam_file_3 }

process index_bam {
    tag { "${params.project_name}.${sample_id}.iB" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file) from bam_file_1

    output:
    set val(sample_id), file("${bam_file}.bai") into bam_index

    script:
    """
    samtools index  \
    ${bam_file} ${bam_file}.bai
    """
}

process run_flagstat {
    tag { "${params.project_name}.${sample_id}.rF" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file) from bam_file_2

    output:
    set val(sample_id), file("${bam_file}.flagstat") into bam_stats

    script:
    """
    samtools flagstat \
    -@ 1 \
    ${bam_file} > ${bam_file}.flagstat  \
    """
}

bam_file_3.mix(bam_index,bam_stats).groupTuple().set{bam_all}

process create_md5sum {
    tag { "${params.project_name}.${sample_id}.cMD5S" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'move', overwrite: false
    input:
    set val(sample_id), file(bam_file) from bam_all

    output:
    set val(sample_id), file(bam_file), file("${bam_file[0]}.md5"), file("${bam_file[1]}.md5") into bam_all_md5sum

    script:
    """
    md5sum ${bam_file[0]} > ${bam_file[0]}.md5
    md5sum ${bam_file[1]} > ${bam_file[1]}.md5
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
