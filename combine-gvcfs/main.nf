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

chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')

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
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rCG" }
    label 'bigmem'
    publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

    input:
    file(gvcf_list)
    each chr from chroms

    output:
    file("${params.cohort_id}.${chr}.g.vcf.gz")  into cohort_chr_calls
    file("${params.cohort_id}.${chr}.g.vcf.gz.tbi") into cohort_chr_indexes

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
    CombineGVCFs \
    -R ${params.ref_seq} \
    --L ${chr} \
    --variant ${gvcf_list} \
    -O ${params.cohort_id}.${chr}.g.vcf.gz
    """
}

cohort_chr_calls.groupTuple().set{cohort_calls}

process run_concat_combine_gvcf {
     tag { "${params.project_name}.${params.cohort_id}.rCCG" }
     cpus { 20 }
     publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

     input:
     set val(params.cohort_id), file(gvcf) from cohort_calls

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
