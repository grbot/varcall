#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
known_indels_1            = file(params.known_indels_1, type: 'file')
known_indels_2            = file(params.known_indels_2, type: 'file')
dbsnp                     = file(params.dbsnp, type: 'file')
hapmap                    = file(params.hapmap, type: 'file')
omni                      = file(params.omni, type: 'file')
phase1_snps               = file(params.phase1_snps, type: 'file')
golden_indels             = file(params.golden_indels, type: 'file')

db                        = file(params.db_path, type: 'dir')

ts_filter_level_snps      = params.ts_filter_level_snps
ts_filter_level_indels    = params.ts_filter_level_indels
max_gaussians_snps        = params.max_gaussians_snps
max_gaussians_indels      = params.max_gaussians_indels
sample_coverage           = params.sample_coverage
project_name              = params.project_name
cohort_id                 = params.cohort_id

outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

// CALL_CONF
if ( sample_coverage == "high" ) {
    call_conf = 30
} else if ( sample_coverage == "low" ) {
    call_conf = 10
} else {
    call_conf = 30
}

process run_genotype_gvcf_on_genome_gvcf {
    tag { "${project_name}.${cohort_id}.${chr}.rGGoG" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true
    
    input:
    tuple path(gvcf), path(index)
    each interval
  
    output:
    tuple path("${cohort_id}.${chr}.vcf.gz"), path("${cohort_id}.${chr}.vcf.gz.tbi"), emit: gg_vcf_set
  
    script:
    mem = task.memory.toGiga() - 4
    if (db_import == "no") {
        variants = "-V ${gvcf}"
    } else if (db_import == "yes") {
        variants = "-V gendb://${db}/${chr}.gdb"
    } else {
        exit 1
    }

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenotypeGVCFs \
        -R ${ref} \
        -L ${chr} \
        ${variants} \
        -stand-call-conf ${call_conf} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        --allow-old-rms-mapping-quality-annotation-data \
        -O "${cohort_id}.${chr}.vcf.gz"
    """
}

process run_concat_vcf {
    tag { "${project_name}.${cohort_id}.rCV" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true

    input:
    tuple path(vcf), path(index)
  
    output:
  	tuple path("${cohort_id}.vcf.gz"), path("${cohort_id}.vcf.gz.tbi"), emit: combined_calls
  
    script:
    mem = task.memory.toGiga() - 4
    if (params.build == 'b37') {
        chr = ''
        mt = 'MT'
    } else if(build == 'b38') {
        chr = 'chr'
        mt = 'chrM'
    }
    
    """
    find . -iname "${cohort_id}.${chr}1.vcf.gz" > ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}2.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}3.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}4.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}5.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}6.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}7.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}8.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}9.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}10.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}11.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}12.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}13.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}14.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}15.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}16.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}17.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}18.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}19.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}20.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}21.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}22.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}X.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${chr}Y.vcf.gz" >> ${cohort_id}.vcf.list
    find . -iname "${cohort_id}.${mt}.vcf.gz" >> ${cohort_id}.vcf.list

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GatherVcfs \
        -I ${cohort_id}.vcf.list \
        -O ${cohort_id}.vcf.gz 
    
    # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
    tabix -p vcf ${cohort_id}.vcf.gz
    """
}

process run_vqsr_on_snps {
    tag { "${project_name}.${cohort_id}.rVoS" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)
 
    output:
    tuple path(vcf), path(vcf_index),
        path("${cohort_id}.vcf.recal-SNP.recal"),
        path("${cohort_id}.vcf.recal-SNP.recal.idx"),
        path("${cohort_id}.vcf.recal-SNP.tranches"), emit: snps_vqsr_recal

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" VariantRecalibrator \
        -R ${ref} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
        -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1_snps} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
        -an DP \
        -an FS \
        -an SOR \
        -an MQ \
        -an MQRankSum \
        -an QD \
        -an ReadPosRankSum \
        -mode SNP \
        --max-gaussians "${max_gaussians_snps}" \
        -V ${vcf} \
        -O "${cohort_id}.vcf.recal-SNP.recal" \
        --tranches-file "${cohort_id}.vcf.recal-SNP.tranches"
    """
}

process run_apply_vqsr_on_snps {
    tag { "${project_name}.${cohort_id}.rAVoS" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index), path(snp_recal), path(snp_recal_index), path(snp_tranches)

    output:
    tuple path("${cohort_id}.recal-SNP.vcf.gz"), path("${cohort_id}.recal-SNP.vcf.gz.tbi"), emit: snps_vqsr_vcf

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" ApplyVQSR \
        -R ${ref} \
        --recal-file ${snp_recal} \
        --tranches-file ${snp_tranches} \
        -mode SNP \
        -ts-filter-level "${ts_filter_level_snps}" \
        -V ${vcf} \
        -O "${cohort_id}.recal-SNP.vcf.gz"
    """
}

process run_vqsr_on_indels {
    tag { "${project_name}.${cohort_id}.rVoI" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path(vcf), path(vcf_index),
    path("${cohort_id}.recal-SNP.vcf.recal-INDEL.recal"),
    path("${cohort_id}.recal-SNP.vcf.recal-INDEL.recal.idx"),
    path("${cohort_id}.recal-SNP.vcf.recal-INDEL.tranches"), emit: indel_vqsr_recal

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" VariantRecalibrator \
        -R ${ref} \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 ${golden_indels} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
        -an DP \
        -an FS \
        -an SOR \
        -an MQ \
        -an MQRankSum \
        -an QD \
        -an ReadPosRankSum \
        -mode INDEL \
        --max-gaussians "${max_gaussians_indels}" \
        -V ${vcf} \
        -O "${cohort_id}.recal-SNP.vcf.recal-INDEL.recal" \
        --tranches-file "${cohort_id}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${project_name}.${cohort_id}.rAVoI" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${cohort_id}", mode: 'copy', overwrite: true

    input:
    tuple path(vcf), path(vcf_index), path(indel_recal), path(indel_recal_index), path(indel_tranches)
 
    output:
    tuple path("${cohort_id}.recal-SNP.recal-INDEL.vcf.gz"), path("${cohort_id}.recal-SNP.recal-INDEL.vcf.gz.tbi"), emit: indel_vqsr_vcf

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" ApplyVQSR \
        -R ${ref} \
        --recal-file ${indel_recal} \
        --tranches-file ${indel_tranches} \
        -mode INDEL \
        -ts-filter-level "${ts_filter_level_indels}" \
        -V ${vcf} \
        -O "${cohort_id}.recal-SNP.recal-INDEL.vcf.gz"
    """
}
