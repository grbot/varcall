#!/usr/bin/env nextflow

if (params.build == "b37") {
  //chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
} else if (params.build == "b38"){
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM".split(',')
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_dict = Channel.fromPath(params.ref_seq_dict).toList()

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVG" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/genotype-refinement", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

if (params.supporting)
  supporting = "-supporting " + params.supporting
else
  supporting = ""

if (params.family)
  family = "-ped " +  params.family
else
  family = ""

Channel.from( file(params.vcf) )
        .set{ vcf_ch }
Channel.from( file(params.vcf_index) )
        .set{ vcf_index_ch }

process run_calculate_genotype_posteriors {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rCGP" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genotype-refinement", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    file vcf from vcf_ch
    file vcf_index from vcf_index_ch
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    each chr from chroms
 
    output:
    set chr, file("${params.cohort_id}.cgp.${chr}.vcf.gz"), file("${params.cohort_id}.cgp.${chr}.vcf.gz.tbi") into cgp_vcf

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    CalculateGenotypePosteriors ${supporting} ${family} \
   -R ${ref_seq} \
   -L ${chr} \
   -V ${vcf} \
   -O "${params.cohort_id}.cgp.${chr}.vcf.gz" \
   """
}

process run_variant_filtration {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rVF" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genotype-refinement", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    set chr, file(vcf), file(vcf_index) from cgp_vcf
    file ref_seq
    file ref_seq_index
    file ref_seq_dict

    output:
    set chr, file("${params.cohort_id}.cgp.vf.${chr}.vcf.gz"), file("${params.cohort_id}.cgp.vf.${chr}.vcf.gz.tbi") into vf_vcf

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantFiltration \
    -R ${params.ref_seq} \
    -L ${chr} \
    --filter-expression "GQ<20" \
    --filter-name "lowGQ" \
    -V ${vcf} \
    -O "${params.cohort_id}.cgp.vf.${chr}.vcf.gz"
    """
}

process run_variant_annotator {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rVA" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genotype-refinement", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    set chr, file(vcf), file(vcf_index) from vf_vcf
    file ref_seq
    file ref_seq_index
    file ref_seq_dict

    output:
    file("${params.cohort_id}.cgp.vf.va.${chr}.vcf.gz") into va_vcf
    file("${params.cohort_id}.cgp.vf.va.${chr}.vcf.gz.tbi") into va_vcf_index

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantAnnotator \
   -R ${params.ref_seq} \
   -L ${chr} \
   -V ${vcf} \
   -A PossibleDeNovo \
   -O "${params.cohort_id}.cgp.vf.va.${chr}.vcf.gz" \
   """
}

va_vcf.toList().set{ concat_ready  }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genotype-refinement", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.cgp.vf.va.vcf.gz"), file("${params.cohort_id}.cgp.vf.va.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.1\\.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.2\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.3\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.4\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.5\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.6\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.7\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.8\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.9\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.10\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.11\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.12\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.13\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.14\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.15\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.16\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.17\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.18\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.19\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.20\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.21\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.22\\.vcf.gz" >> ${params.cohort_id}.vcf.list
 /*      
       echo "${vcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.Y\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.MT\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list 
 */      
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.cgp.vf.va.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.cgp.vf.va.vcf.gz
       """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build38 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genotype-gvcf", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.cgp.vf.va.vcf.gz"), file("${params.cohort_id}.cgp.vf.va.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.chr1\\.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr2\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr3\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr4\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr5\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr6\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr7\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr8\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr9\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr10\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr11\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr12\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr13\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr14\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr15\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr16\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr17\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr18\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr19\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr20\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr21\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr22\\.vcf.gz" >> ${params.cohort_id}.vcf.list
/*
       echo "${vcf.join('\n')}" | grep "\\.chrX\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrY\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrM\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
*/       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.cgp.vf.va.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.cgp.vf.va.vcf.gz
       """
  }
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
