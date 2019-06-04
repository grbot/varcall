#!/usr/bin/env nextflow

chroms = params.chromosomes.split(',')

Channel.from( file(params.gvcf_file) )
        .set{ gvcf_file_ch }

process run_genotype_gvcf_on_genome {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rGGoG" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:		
    val (gvcf_file) from gvcf_file_ch
    each chr from chroms

    output:
    set chr, file("${params.cohort_id}.${chr}.vcf.gz"), file("${params.cohort_id}.${chr}.vcf.gz.tbi") into gg_vcf

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GenotypeGVCFs \
    -R ${params.ref_seq} \
    -L $chr \
    -V ${gvcf_file} \
    -stand-call-conf ${call_conf} \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -O "${params.cohort_id}.${chr}.vcf.gz" 
    """
}

process run_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rVoS" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:		
    set val (chr), file(vcf), file(vcf_index) from gg_vcf

    output:
    set chr, file(vcf), file(vcf_index), file("${params.cohort_id}.${chr}.vcf.recal-SNP.recal"), file("${params.cohort_id}.${chr}.vcf.recal-SNP.recal.idx"), file("${params.cohort_id}.${chr}.vcf.recal-SNP.tranches") into snps_vqsr_recal

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    VariantRecalibrator \
   -R ${params.ref_seq} \
   -L $chr \
   -resource hapmap,known=false,training=true,truth=true,prior=15.0:${params.hapmap} \
   -resource omni,known=false,training=true,truth=true,prior=12.0:${params.omni} \
   -resource 1000G,known=false,training=true,truth=false,prior=10.0:${params.phase1_snps} \
   -resource dbsnp,known=true,training=false,truth=false,prior=2.0:${params.dbsnp} \
   -an DP \
   -an FS \
   -an SOR \
   -an MQ \
   -an MQRankSum \
   -an QD \
   -an ReadPosRankSum \
   -mode SNP \
    --max-gaussians "${params.max_gaussians_snps}" \
   -V ${vcf} \
   -O "${params.cohort_id}.${chr}.vcf.recal-SNP.recal" \
   --tranches-file "${params.cohort_id}.${chr}.vcf.recal-SNP.tranches"
   """
}

process run_apply_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rAVoS" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:		
    set val (chr), file(vcf), file(vcf_index), file(snp_recal), file(snp_recal_index), file(snp_tranches) from snps_vqsr_recal

    output:
    set chr, file("${params.cohort_id}.${chr}.recal-SNP.vcf.gz"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.gz.tbi") into snps_vqsr_vcf

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    ApplyVQSR \
    -R ${params.ref_seq} \
    --recal-file ${snp_recal} \
    --tranches-file ${snp_tranches} \
    -mode SNP \
    -ts-filter-level "${params.ts_filter_level_snps}" \
    -V ${vcf} \
    -O "${params.cohort_id}.${chr}.recal-SNP.vcf.gz"
    """
}

process run_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rVoI" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:
    set val (chr), file(vcf), file(vcf_index) from snps_vqsr_vcf

    output:
    set chr, file(vcf), file(vcf_index), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal.idx"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.tranches") into indel_vqsr_recal

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    VariantRecalibrator \
   -R ${params.ref_seq} \
   -L $chr \
   -resource mills,known=false,training=true,truth=true,prior=12.0:${params.golden_indels} \
   -resource dbsnp,known=true,training=false,truth=false,prior=2.0:${params.dbsnp} \
   -an DP \
   -an FS \
   -an SOR \
   -an MQ \
   -an MQRankSum \
   -an QD \
   -an ReadPosRankSum \
   -mode INDEL \
    --max-gaussians "${params.max_gaussians_indels}" \
   -V ${vcf} \
   -O "${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal" \
   --tranches-file "${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rAVoI" }
    label "bigmem"
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    input:		
    set val (chr), file(vcf), file(vcf_index), file(indel_recal), file(indel_recal_index), file(indel_tranches) from indel_vqsr_recal

    output:
    set chr, file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz.tbi") into indel_vqsr_vcf
    file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz") into vcf_concat_ready

    script:
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    ApplyVQSR \
    -R ${params.ref_seq} \
    --recal-file ${indel_recal} \
    --tranches-file ${indel_tranches} \
    -mode INDEL \
    -ts-filter-level "${params.ts_filter_level_indels}" \
    -V ${vcf} \
    -O "${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz"
    """
}

vcf_concat_ready.toList().set{ concat_ready  }

process run_concat_vcf {
     tag { "${params.project_name}.${params.cohort_id}.rCV" }
     label "bigmem"
     publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
 
     input:
     file(vcf) from concat_ready

     output:
	   set val("${params.cohort_id}"), file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz.tbi") into combine_calls

     script:
     """
     echo "${vcf.join('\n')}" | grep "\\.1\\.*recal-SNP.recal-INDEL.vcf.gz" > ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.2\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.3\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.4\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.5\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.6\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.7\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.8\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.9\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.10\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.11\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.12\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.13\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.14\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.15\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.16\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.17\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.18\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.19\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.20\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.21\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
     echo "${vcf.join('\n')}" | grep "\\.22\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
    
     ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
     GatherVcfs \
     -I ${params.cohort_id}.vcf.list \
     -O ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     ${params.tabix_base}/tabix -p vcf ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz 
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
