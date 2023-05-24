#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
build                     = params.build
project_name              = params.project_name
cohort_id                 = params.cohort_id
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_combine_gvcfs {
    tag { "${project_name}.${chr}.rCG" }
    label 'gatk'
    memory { 64.GB * task.attempt }
    time '120h'
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true

    input:
    path(gvcf_list)
    each chr

    output:
    tuple path("${cohort_id}.${chr}.g.vcf.gz"), path("${cohort_id}.${chr}.g.vcf.gz.tbi"), emit: cohort_chr_calls

    script:
    mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" CombineGVCFs \
        -R ${ref} \
        --L ${chr} \
        --variant ${gvcf_list} \
        -O ${cohort_id}.${chr}.g.vcf.gz
    """
}

process run_concat_gvcfs {
    tag { "${project_name}.rCV" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    time '120h'
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true

    input:
    tuple path(gvcf), path(index)
  
    output:
  	tuple path("${cohort_id}.g.vcf.gz"), file("${cohort_id}.g.vcf.gz.tbi"), emit: combined_calls
  
    script:
    mem = task.memory.toGiga() - 4
    if (build == 'b37') {
        chr = ''
        mt = 'MT'
    } else if(build == 'b38') {
        chr = 'chr'
        mt = 'chrM'
    }
    
    """
    find . -iname "${cohort_id}.${chr}1.g.vcf.gz" > ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}2.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}3.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}4.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}5.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}6.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}7.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}8.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}9.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}10.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}11.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}12.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}13.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}14.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}15.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}16.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}17.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}18.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}19.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}20.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}21.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}22.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}X.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${chr}Y.g.vcf.gz" >> ${cohort_id}.gvcf.list
    find . -iname "${cohort_id}.${mt}.g.vcf.gz" >> ${cohort_id}.gvcf.list

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GatherVcfs \
        -I ${cohort_id}.gvcf.list \
        -O ${cohort_id}.g.vcf.gz
 
    # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
    tabix -p vcf ${cohort_id}.g.vcf.gz
    """
}
