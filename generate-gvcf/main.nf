#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gender = row['Gender']
            def bam_file = file(row['BAM'])
            return [ sample_id, gender, bam_file ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')

process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), val(gender), file(bam_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tGender: ${gender}\tBAM: ${bam_file}\n"
    """
}

process run_haplotype_caller_on_autosomes {
    tag { "${params.project_name}.${sample_id}.${chr}.rHCoA" }
    memory { 4.GB * task.attempt }
    cpus { 1 }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), val(bam_file) from samples_2
	  each chr from autosomes

    output:
	  set val(sample_id), file("${sample_id}.${chr}.g.vcf.gz")  into autosome_calls
	  set val(sample_id), file("${sample_id}.${chr}.g.vcf.gz.tbi") into autosome_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
    HaplotypeCaller \
    -R ${params.ref_seq} \
    -I $bam_file \
    --emit-ref-confidence GVCF \
    --dbsnp ${params.dbsnp_sites} \
    --L $chr \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.${chr}.g.vcf.gz
    """
}

// Now do X and Y calling
samples_3.filter{it[1] == 'M'}.into{samples_male_1; samples_male_2; samples_male_3; samples_male_4; samples_male_5; samples_male_6}
samples_4.filter{it[1] == 'F'}.set{samples_female}

// Males
process run_haplotype_caller_on_x_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP1M" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_1

     output:
	   set val(sample_id), file("${sample_id}.X_PAR1.g.vcf.gz") into x_par1_calls
	   set val(sample_id), file("${sample_id}.X_PAR1.g.vcf.gz.tbi") into x_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L X:60001-2699520 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR1.g.vcf.gz
     """
}

process run_haplotype_caller_on_x_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP2M" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_2

     output:
	   set val(sample_id), file("${sample_id}.X_PAR2.g.vcf.gz") into x_par2_calls
	   set val(sample_id), file("${sample_id}.X_PAR2.g.vcf.gz.tbi") into x_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L X:154931044-155260560 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR2.g.vcf.gz
     """
}

process run_haplotype_caller_on_x_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoXNPM" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_3

     output:
	   set val(sample_id), file("${sample_id}.X_nonPAR.g.vcf.gz") into x_nonpar_calls
	   set val(sample_id), file("${sample_id}.X_nonPAR.g.vcf.gz.tbi") into x_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L X -XL X:60001-2699520 -XL X:154931044-155260560 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.X_nonPAR.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP1M" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_4

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR1.g.vcf.gz") into y_par1_calls
	   set val(sample_id), file("${sample_id}.Y_PAR1.g.vcf.gz.tbi") into y_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L Y:10001-2649520 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR1.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP2M" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_5

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR2.g.vcf.gz") into y_par2_calls
	   set val(sample_id), file("${sample_id}.Y_PAR2.g.vcf.gz.tbi") into y_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L Y:59034050-59363566 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR2.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoYNPM" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_6

     output:
	   set val(sample_id), file("${sample_id}.Y_nonPAR.g.vcf.gz") into y_nonpar_calls
	   set val(sample_id), file("${sample_id}.Y_nonPAR.g.vcf.gz.tbi") into y_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L Y -XL Y:10001-2649520 -XL Y:59034050-59363566 \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.Y_nonPAR.g.vcf.gz
    """
}

// Females
process run_haplotype_caller_on_x_female {
     tag { "${params.project_name}.${sample_id}.rHCoXF" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_female

     output:
	   set val(sample_id), file("${sample_id}.X.g.vcf.gz") into x_calls
	   set val(sample_id), file("${sample_id}.X.g.vcf.gz.tbi") into x_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
     HaplotypeCaller \
     -R ${params.ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${params.dbsnp_sites} \
     --L X \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X.g.vcf.gz
     """

}

process run_haplotype_caller_on_mt {
    tag { "${params.project_name}.${sample_id}.rHCoMT" }

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), file(bam_file) from samples_5

    output:
	  set val(sample_id), file("${sample_id}.MT.g.vcf.gz") into mt_calls
	  set val(sample_id), file("${sample_id}.MT.g.vcf.gz.tbi") into mt_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    """
    ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_hc_mem}"  \
    HaplotypeCaller \
    -R ${params.ref_seq} \
    -I $bam_file \
    --emit-ref-confidence GVCF \
    --dbsnp ${params.dbsnp_sites} \
    --L MT \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.MT.g.vcf.gz
    """
}

autosome_calls.mix(mt_calls,x_par1_calls,x_nonpar_calls,x_par2_calls,x_calls,y_par1_calls,y_nonpar_calls,y_par2_calls).groupTuple().set{all_calls}

process combine_gVCFs {
     tag { "${params.project_name}.${sample_id}.cCgVCF" }

     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), file(gvcf) from all_calls

     output:
	   set val(sample_id), file("${sample_id}.g.vcf.gz") into combine_calls
	   set val(sample_id), file("${sample_id}.g.vcf.gz.tbi") into combine_calls_indexes

     script:
     if (gvcf.size() == 29) // working with a male sample
     """
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_sv_mem}"  \
     SortVcf \
     -I ${sample_id}.X_PAR1.g.vcf.gz \
     -I ${sample_id}.X_PAR2.g.vcf.gz \
     -I ${sample_id}.X_nonPAR.g.vcf.gz \
     -O ${sample_id}.X.g.vcf.gz

     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_sv_mem}"  \
     SortVcf \
     -I ${sample_id}.Y_PAR1.g.vcf.gz \
     -I ${sample_id}.Y_PAR2.g.vcf.gz \
     -I ${sample_id}.Y_nonPAR.g.vcf.gz \
     -O ${sample_id}.Y.g.vcf.gz

     echo "${gvcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo ${sample_id}.X.g.vcf.gz >> ${sample_id}.gvcf.list
     echo ${sample_id}.Y.g.vcf.gz >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.MT\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list
    
     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_gv_mem}"  \
     GatherVcfs \
     -I ${sample_id}.gvcf.list \
     -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.

     ${params.tabix_base}/tabix -p vcf ${sample_id}.g.vcf.gz 
     """
     else if (gvcf.size() == 24) // working with a female  sample
     """
     echo "${gvcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${sample_id}.gvcf.list
     echo "${gvcf.join('\n')}" | grep "\\.MT\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list

     ${params.gatk_base}/gatk --java-options "-Xmx${params.gatk_gv_mem}"  \
     GatherVcfs \
     -I ${sample_id}.gvcf.list \
     -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
     
     ${params.tabix_base}/tabix -p vcf ${sample_id}.g.vcf.gz
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
