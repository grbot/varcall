#!/usr/bin/env nextflow

chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')

process run_remove_samples_combined_gvcf {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rRSCG" }
    label 'bigmem'
    publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

    input:
    file (gvcf) from file(params.combined_gvcf)
    file (gvcf_index) from file(params.combined_gvcf_index)
    each chr from chroms

    output:
    file("${params.cohort_id}.${chr}.g.vcf.gz")  into cohort_chr_calls
    file("${params.cohort_id}.${chr}.g.vcf.gz.tbi") into cohort_chr_indexes

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
    SelectVariants \
    -R ${params.ref_seq} \
    --L ${chr} \
    -V ${gvcf} \
    --exclude-sample-name ${params.exclude_samples} \
    -O ${params.cohort_id}.${chr}.g.vcf.gz
    """
}

cohort_chr_calls.toList().into{ cohort_calls }

process run_concat_combine_gvcf {
     tag { "${params.project_name}.${params.cohort_id}.rCCG" }
     cpus { 20 }
     publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

     input:
     file(gvcf) from cohort_calls

     output:
	   set val(params.cohort_id), file("${params.cohort_id}.g.vcf.gz") into combine_calls
	   set val(params.cohort_id), file("${params.cohort_id}.g.vcf.gz.tbi") into combine_calls_indexes

     script:
     """
     echo "${gvcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.Y\\.g.vcf.gz" >> ${params.cohort_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.MT\\.g\\.vcf\\.gz" >> ${params.cohort_id}.gvcf.list
    
     ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
     GatherVcfs \
     -I ${params.cohort_id}.gvcf.list \
     -O ${params.cohort_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     ${params.tabix_base}/tabix -p vcf ${params.cohort_id}.g.vcf.gz 
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
