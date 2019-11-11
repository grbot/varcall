#!/usr/bin/env nextflow

Channel.from( file(params.gvcf_file) )
        .set{ gvcf_file_ch }

if (params.build == "b37") {
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')
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
dbsnp = Channel.fromPath(params.dbsnp).toList()
dbsnp_index = Channel.fromPath(params.dbsnp_index).toList()
hapmap = Channel.fromPath(params.hapmap).toList()
hapmap_index = Channel.fromPath(params.hapmap_index).toList()
omni = Channel.fromPath(params.omni).toList()
omni_index = Channel.fromPath(params.omni_index).toList()
phase1_snps = Channel.fromPath(params.phase1_snps).toList()
phase1_snps_index = Channel.fromPath(params.phase1_snps_index).toList()
golden_indels = Channel.fromPath(params.golden_indels).toList()
golden_indels_index = Channel.fromPath(params.golden_indels_index).toList()

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVG" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process run_genotype_gvcf_on_genome {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rGGoG" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
  
    input:
    val (gvcf_file) from gvcf_file_ch
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    each chr from chroms

    output:
    set chr, file("${params.cohort_id}.${chr}.vcf.gz"), file("${params.cohort_id}.${chr}.vcf.gz.tbi") into gg_vcf_set
    file("${params.cohort_id}.${chr}.vcf.gz") into gg_vcf

    script:
      mem = task.memory.toGiga() - 4
      call_conf = 30 // set default
      if ( params.sample_coverage == "high" )
        call_conf = 30
      else if ( params.sample_coverage == "low" )
        call_conf = 10
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    GenotypeGVCFs \
    -R ${ref_seq} \
    -L $chr \
    -V ${gvcf_file} \
    -stand-call-conf ${call_conf} \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -O "${params.cohort_id}.${chr}.vcf.gz"
    """
}

gg_vcf.toList().set{ concat_ready  }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.vcf.gz"), file("${params.cohort_id}.vcf.gz.tbi") into combined_calls
  
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
       echo "${vcf.join('\n')}" | grep "\\.X\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.Y\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.MT\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.vcf.gz
       """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.vcf.gz"), file("${params.cohort_id}.vcf.gz.tbi") into combined_calls
  
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
       echo "${vcf.join('\n')}" | grep "\\.chrX\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrY\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrM\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.vcf.gz
       """
  }
}

process run_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.rVoS" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    set file(vcf), file(vcf_index) from combined_calls
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    file hapmap
    file hapmap_index
    file omni
    file omni_index
    file phase1_snps
    file phase1_snps_index
    file dbsnp
    file dbsnp_index
 
    output:
    set file(vcf), file(vcf_index), file("${params.cohort_id}.vcf.recal-SNP.recal"), file("${params.cohort_id}.vcf.recal-SNP.recal.idx"), file("${params.cohort_id}.vcf.recal-SNP.tranches") into snps_vqsr_recal

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantRecalibrator \
   -R ${ref_seq} \
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
    --max-gaussians "${params.max_gaussians_snps}" \
   -V ${vcf} \
   -O "${params.cohort_id}.vcf.recal-SNP.recal" \
   --tranches-file "${params.cohort_id}.vcf.recal-SNP.tranches"
   """
}

process run_apply_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.rAVoS" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    set file(vcf), file(vcf_index), file(snp_recal), file(snp_recal_index), file(snp_tranches) from snps_vqsr_recal
    file ref_seq
    file ref_seq_index
    file ref_seq_dict

    output:
    set file("${params.cohort_id}.recal-SNP.vcf.gz"), file("${params.cohort_id}.recal-SNP.vcf.gz.tbi") into snps_vqsr_vcf

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    ApplyVQSR \
    -R ${params.ref_seq} \
    --recal-file ${snp_recal} \
    --tranches-file ${snp_tranches} \
    -mode SNP \
    -ts-filter-level "${params.ts_filter_level_snps}" \
    -V ${vcf} \
    -O "${params.cohort_id}.recal-SNP.vcf.gz"
    """
}

process run_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.rVoI" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    set file(vcf), file(vcf_index) from snps_vqsr_vcf
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    file golden_indels
    file golden_indels_index
    file dbsnp
    file dbsnp_index

    output:
    set file(vcf), file(vcf_index), file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal"), file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal.idx"), file("${params.cohort_id}.recal-SNP.vcf.recal-INDEL.tranches") into indel_vqsr_recal

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantRecalibrator \
   -R ${params.ref_seq} \
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
    --max-gaussians "${params.max_gaussians_indels}" \
   -V ${vcf} \
   -O "${params.cohort_id}.recal-SNP.vcf.recal-INDEL.recal" \
   --tranches-file "${params.cohort_id}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.rAVoI" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'

    input:
    set file(vcf), file(vcf_index), file(indel_recal), file(indel_recal_index), file(indel_tranches) from indel_vqsr_recal
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
 
    output:
    set file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz.tbi") into indel_vqsr_vcf

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    ApplyVQSR \
    -R ${ref_seq} \
    --recal-file ${indel_recal} \
    --tranches-file ${indel_tranches} \
    -mode INDEL \
    -ts-filter-level "${params.ts_filter_level_indels}" \
    -V ${vcf} \
    -O "${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz"
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
