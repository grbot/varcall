#!/usr/bin/env nextflow

if (params.build == "b37") {
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')
} else if (params.build == "b38"){
    // chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM".split(',')
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22".split(',')
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

Channel.fromPath( file(params.intervals) )
        .splitText()
        .map{it -> it.trim()}
        .set {intervals}

intervals_file = Channel.fromPath( file(params.intervals) ).toList()

db = file(params.db_path)
db_import = params.db_import

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

process log_tool_version_bcftools {
    tag { "${params.project_name}.ltViB" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'move', overwrite: false
    label 'bcftools'

    output:
    file("tool.bcftools.version") into tool_version_bcftools

    script:
    mem = task.memory.toGiga() - 3
    """
    bcftools --version > tool.bcftools.version 2>&1
    """
}

process run_genotype_gvcf_on_genome_db {
    tag { "${params.project_name}.${params.cohort_id}.${interval}.rGGoG" }
    memory { 48.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    time = 24.h
 
    input:
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    file db
    each interval from intervals  
    
    output:
    file("*") into gg_vcf

    script:
      mem = task.memory.toGiga() - 16 
      call_conf = 30 // set default
      if ( params.sample_coverage == "high" )
        call_conf = 30
      else if ( params.sample_coverage == "low" )
        call_conf = 10
    """
    a=$interval
    b=\${a/:/_}
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    GenotypeGVCFs \
    -R ${ref_seq} \
    -L $interval \
    -V gendb://${db}/\$b.gdb \
    -stand-call-conf ${call_conf} \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    --allow-old-rms-mapping-quality-annotation-data \
    -O "${params.cohort_id}.\$b.vcf.gz"
    """
}

gg_vcf.flatten().toList().set{ concat_ready }

process create_vcf_list {
    tag { "${params.project_name}.${params.cohort_id}.rCVL" }
    memory { 4.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
  
    input:
    file intervals_file  
    
    output:
    file("${params.cohort_id}.vcf.list") into vcf_list

    script:
    """
    while read line; do
    a=\$line
    b=\${a/:/_}
    echo "${params.cohort_id}."\$b".vcf.gz"
    done < $intervals_file > ${params.cohort_id}.vcf.list;
    """
}

process run_concat_vcf{
     tag { "${params.project_name}.${params.cohort_id}.rCV" }
     memory { 16.GB * task.attempt }  
     publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
     label 'gatk'
     time = 336.h

     input:
     //set file(vcf), file(vcf_index) from concat_ready
     file(vcf) from concat_ready
     file vcf_list

     output:
	   file("${params.cohort_id}.vcf.gz") into vcf_concat
	   file("${params.cohort_id}.vcf.gz.tbi") into vcf_concat_index

     script:
       mem = task.memory.toGiga() - 4
     """
     gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
     GatherVcfs \
     -I ${vcf_list} \
     -O ${params.cohort_id}.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     tabix -p vcf ${params.cohort_id}.vcf.gz
     """
}

process run_select_chromosomes {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rSC" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'bcftools'
    time = 336.h
    input:
    file(vcf) from vcf_concat
    file(vcf_index) from vcf_concat_index
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    each chr from chroms

    output:
    set chr, file("${params.cohort_id}.${chr}.vcf.gz"), file("${params.cohort_id}.${chr}.vcf.gz.tbi") into chr_vcf

    script:
    """
    bcftools \
    view \
    -O z \
    --regions $chr \
    -o "${params.cohort_id}.${chr}.vcf.gz" \
    ${vcf}
    bcftools index -t ${params.cohort_id}.${chr}.vcf.gz
    """
}

process run_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rVoS" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    time = 336.h    
    input:
    set chr, file(vcf), file(vcf_index) from chr_vcf
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
    set chr, file(vcf), file(vcf_index), file("${params.cohort_id}.${chr}.vcf.recal-SNP.recal"), file("${params.cohort_id}.${chr}.vcf.recal-SNP.recal.idx"), file("${params.cohort_id}.${chr}.vcf.recal-SNP.tranches") into snps_vqsr_recal

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantRecalibrator \
   -R ${ref_seq} \
   -L ${chr} \
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
   -O "${params.cohort_id}.${chr}.vcf.recal-SNP.recal" \
   --tranches-file "${params.cohort_id}.${chr}.vcf.recal-SNP.tranches"
   """
}

process run_apply_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rAVoS" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    time = 336.h    
   
    input:
    set chr, file(vcf), file(vcf_index), file(snp_recal), file(snp_recal_index), file(snp_tranches) from snps_vqsr_recal
    file ref_seq
    file ref_seq_index
    file ref_seq_dict

    output:
    set chr, file("${params.cohort_id}.${chr}.recal-SNP.vcf.gz"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.gz.tbi") into snps_vqsr_vcf

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    ApplyVQSR \
    -R ${params.ref_seq} \
    -L ${chr} \
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
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    time = 336.h    
    
    input:
    set chr, file(vcf), file(vcf_index) from snps_vqsr_vcf
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    file golden_indels
    file golden_indels_index
    file dbsnp
    file dbsnp_index

    output:
    set chr, file(vcf), file(vcf_index), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal.idx"), file("${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.tranches") into indel_vqsr_recal

    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    VariantRecalibrator \
   -R ${params.ref_seq} \
   -L $chr \
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
   -O "${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.recal" \
   --tranches-file "${params.cohort_id}.${chr}.recal-SNP.vcf.recal-INDEL.tranches"
   """
}

process run_apply_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rAVoI" }
    memory { 16.GB * task.attempt }  
    publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'
    time = 336.h    

    input:
    set chr, file(vcf), file(vcf_index), file(indel_recal), file(indel_recal_index), file(indel_tranches) from indel_vqsr_recal
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
 
    output:
    set chr, file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz.tbi") into indel_vqsr_vcf_set
    file("${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz") into indel_vqsr_vcf


    script:
       mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    ApplyVQSR \
    -R ${ref_seq} \
    -L ${chr} \
    --recal-file ${indel_recal} \
    --tranches-file ${indel_tranches} \
    -mode INDEL \
    -ts-filter-level "${params.ts_filter_level_indels}" \
    -V ${vcf} \
    -O "${params.cohort_id}.${chr}.recal-SNP.recal-INDEL.vcf.gz"
    """
}

indel_vqsr_vcf.toList().set{ concat_ready  }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
       label 'gatk'
       time = 336.h    
  
       input:
       file(vcf) from concat_ready
  
       output:
           set file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz.tbi") into combine_calls

       script:
         mem = task.memory.toGiga() - 4
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
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz 
      """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build38 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/genome-calling", mode: 'copy', overwrite: false
       label 'gatk'
       time = 336.h    
       
       input:
       file(vcf) from concat_ready
  
       output:
           set file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz"), file("${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz.tbi") into combine_calls

       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.chr1\\.*recal-SNP.recal-INDEL.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr2\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr3\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr4\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr5\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr6\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr7\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr8\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr9\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr10\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr11\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr12\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr13\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr14\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr15\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr16\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr17\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr18\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr19\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr20\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr21\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr22\\.*recal-SNP.recal-INDEL.vcf.gz" >> ${params.cohort_id}.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.recal-SNP.recal-INDEL.vcf.gz 
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
