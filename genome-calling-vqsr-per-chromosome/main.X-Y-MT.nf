#!/usr/bin/env nextflow

chromosomes = "X,Y,MT"

chroms = chromosomes.split(',')

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
