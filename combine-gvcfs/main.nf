#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and gVCF file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gvcf_file = file(row['gVCF'])
            return [ sample_id, gvcf_file ]
        }
        .tap{samples_1; samples_2}
        .map { sample_id, gvcf_file ->
            return [ gvcf_file ]
        }
        .collect().set { gvcf_files }

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    
    input:
    set val(sample_id), file(gvcf_file) from samples_1
    
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tgVCF: ${gvcf_file}\n"
    """
}

process create_variant_list {
    tag { "${params.project_name}.${params.cohort_id}.cVL" }
    memory { 4.GB * task.attempt }
    cpus { 4 }
    publishDir "${params.out_dir}/${params.cohort_id}/combine-gvcfs", mode: 'copy', overwrite: false

    input:
    val gvcf_file from gvcf_files
    
    output:
    file("gvcf.list") into gvcf_list
    
    script:
    """
    echo "${gvcf_file.join('\n')}" > gvcf.list
    """
}

process run_combine_gvcfs {
    tag { "${params.project_name}.${params.cohort_id}.rCG" }
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

    input:
    file(gvcf_list)

    output:
    file("${params.cohort_id}.g.vcf.gz")  into cohort_calls
    file("${params.cohort_id}.g.vcf.gz.tbi") into cohort_indexes

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_cg_mem}"  \
    CombineGVCFs \
    -R ${params.ref_seq} \
    --variant ${gvcf_list} \
    -O ${params.cohort_id}.g.vcf.gz
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
