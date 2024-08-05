#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
known_indels_1            = file(params.known_indels_1, type: 'file')
known_indels_2            = file(params.known_indels_2, type: 'file')
dbsnp                     = file(params.dbsnp, type: 'file')
outdir                    = file(params.outdir, type: 'dir')

build                     = params.build
sample_coverage           = params.sample_coverage
project_name              = params.project_name
cohort_id                 = params.cohort_id

db                        = file(params.db_path, type: 'dir')

hapmap                    = file(params.hapmap, type: 'file')
omni                      = file(params.omni, type: 'file')
phase1_snps               = file(params.phase1_snps, type: 'file')
golden_indels             = file(params.golden_indels, type: 'file')

ts_filter_level_snps      = params.ts_filter_level_snps
ts_filter_level_indels    = params.ts_filter_level_indels
max_gaussians_snps        = params.max_gaussians_snps
max_gaussians_indels      = params.max_gaussians_indels

outdir.mkdir()

// CALL_CONF
if ( sample_coverage == "high" ) {
    call_conf = 30
} else if ( sample_coverage == "low" ) {
    call_conf = 10
} else {
    call_conf = 30
}

process run_genotype_gvcf_on_genome_db {
    tag { "${project_name}.${cohort_id}.${interval}.rGGoG" }
    memory { 48.GB * task.attempt }
    publishDir "${outdir}/${params.workflow}/${project_name}/${interval.replaceAll(":","_")}_vcf", mode: 'copy', overwrite: true
    label 'gatk'
    time = 24.h
    
    input:
    each interval
    
    output:
    tuple val("${project_name}"), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: db_gg_vcf_set

    script:
    mem = task.memory.toGiga() - 16 
    
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenotypeGVCFs \
        --reference ${ref} \
        --intervals ${interval} \
        --variant gendb://${db}/${interval.replaceAll(":","_")}.gdb \
        -stand-call-conf ${call_conf} \
        --annotation Coverage \
        --annotation FisherStrand \
        --annotation StrandOddsRatio \
        --annotation MappingQualityRankSumTest \
        --annotation QualByDepth \
        --annotation RMSMappingQuality \
        --annotation ReadPosRankSumTest \
        --allow-old-rms-mapping-quality-annotation-data \
        --output "${project_name}.${interval.replaceAll(":","_")}.vcf.gz"
    """
}

// process run_genotype_gvcf_on_genome_gvcf {
//     tag { "${project_name}.${cohort_id}.${chr}.rGGoG" }
//     label 'gatk'
//     memory { 16.GB * task.attempt }
//     publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true
    
//     input:
//     tuple path(gvcf), path(index)
//     each interval
  
//     output:
//     tuple path("${project_name}.${chr}.vcf.gz"), path("${project_name}.${chr}.vcf.gz.tbi"), emit: gg_vcf_set
  
//     script:
//     mem = task.memory.toGiga() - 4
//     if (db_import == "no") {
//         variants = "-V ${gvcf}"
//     } else if (db_import == "yes") {
//         variants = "-V gendb://${db}/${chr}.gdb"
//     } else {
//         exit 1
//     }

//     """
//     gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenotypeGVCFs \
//         --reference ${ref} \
//         --intervals ${chr} \
//         ${variants} \
//         -stand-call-conf ${call_conf} \
//         --annotation Coverage \
//         --annotation FisherStrand \
//         --annotation StrandOddsRatio \
//         --annotation MappingQualityRankSumTest \
//         --annotation QualByDepth \
//         --annotation RMSMappingQuality \
//         --annotation ReadPosRankSumTest \
//         --allow-old-rms-mapping-quality-annotation-data \
//         -output "${project_name}.${chr}.vcf.gz"
//     """
// }

process run_concat_vcf {
    tag { "${project_name}.${project_name}.rCV" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true

    input:
    tuple val(cohort), path(vcf), path(index)
  
    output:
  	tuple path("${project_name}.vcf.gz"), path("${project_name}.vcf.gz.tbi"), emit: combined_calls

    
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
    find . -iname "${project_name}.${chr}1.vcf.gz" > ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}2.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}3.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}4.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}5.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}6.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}7.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}8.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}9.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}10.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}11.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}12.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}13.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}14.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}15.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}16.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}17.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}18.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}19.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}20.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}21.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}22.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}X.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${chr}Y.vcf.gz" >> ${project_name}.vcf.list
    find . -iname "${project_name}.${mt}.vcf.gz" >> ${project_name}.vcf.list

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GatherVcfs \
        --INPUT ${project_name}.vcf.list \
        --OUTPUT ${project_name}.vcf.gz 
    
    # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
    tabix -p vcf ${project_name}.vcf.gz
    """
}

process run_vqsr_on_snps {
    tag { "${project_name}.${cohort_id}.rVoS" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)
 
    output:
    tuple path(vcf), path(vcf_index),
        path("${project_name}.vcf.recal-SNP.recal"),
        path("${project_name}.vcf.recal-SNP.recal.idx"),
        path("${project_name}.vcf.recal-SNP.tranches"), emit: snps_vqsr_recal

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" VariantRecalibrator \
        --reference ${ref} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1_snps} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
        --use-annotation DP \
        --use-annotation FS \
        --use-annotation SOR \
        --use-annotation MQ \
        --use-annotation MQRankSum \
        --use-annotation QD \
        --use-annotation ReadPosRankSum \
        --mode SNP \
        --max-gaussians "${max_gaussians_snps}" \
        --variant ${vcf} \
        --output "${project_name}.vcf.recal-SNP.recal" \
        --tranches-file "${project_name}.vcf.recal-SNP.tranches"
    """
}
 
process run_apply_vqsr_on_snps {
    tag { "${project_name}.${cohort_id}.rAVoS" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index), path(snp_recal), path(snp_recal_index), path(snp_tranches)

    output:
    tuple path("${project_name}.recal-SNP.vcf.gz"), path("${project_name}.recal-SNP.vcf.gz.tbi"), emit: snps_vqsr_vcf

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" ApplyVQSR \
        --reference ${ref} \
        --recal-file ${snp_recal} \
        --tranches-file ${snp_tranches} \
        --mode SNP \
        --truth-sensitivity-filter-level "${ts_filter_level_snps}" \
        --variant ${vcf} \
        --output "${project_name}.recal-SNP.vcf.gz"
    """
}

process run_vqsr_on_indels {
    tag { "${project_name}.${cohort_id}.rVoI" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path(vcf), path(vcf_index),
    path("${project_name}.recal-SNP.vcf.recal-INDEL.recal"),
    path("${project_name}.recal-SNP.vcf.recal-INDEL.recal.idx"),
    path("${project_name}.recal-SNP.vcf.recal-INDEL.tranches"), emit: indel_vqsr_recal

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" VariantRecalibrator \
        --reference ${ref} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${golden_indels} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
        --use-annotation DP \
        --use-annotation FS \
        --use-annotation SOR \
        --use-annotation MQ \
        --use-annotation MQRankSum \
        --use-annotation QD \
        --use-annotation ReadPosRankSum \
        --mode INDEL \
        --max-gaussians "${max_gaussians_indels}" \
        --variant ${vcf} \
        --output "${project_name}.recal-SNP.vcf.recal-INDEL.recal" \
        --tranches-file "${project_name}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${project_name}.${cohort_id}.rAVoI" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true

    input:
    tuple path(vcf), path(vcf_index), path(indel_recal), path(indel_recal_index), path(indel_tranches)
 
    output:
    tuple path("${project_name}.recal-SNP.recal-INDEL.vcf.gz"), path("${project_name}.recal-SNP.recal-INDEL.vcf.gz.tbi"), emit: indel_vqsr_vcf

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" ApplyVQSR \
        --reference ${ref} \
        --recal-file ${indel_recal} \
        --tranches-file ${indel_tranches} \
        --mode INDEL \
        --truth-sensitivity-filter-level "${ts_filter_level_indels}" \
        --variant ${vcf} \
        --output "${project_name}.recal-SNP.recal-INDEL.vcf.gz"
    """
}

