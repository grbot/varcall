#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { select_samples; select_samples as select_samples_2 } from './modules/compare/main.nf'
include { get_stats } from './modules/compare/main.nf'
include { get_unique_sites_and_overlaps; get_unique_sites_and_overlap_index } from './modules/compare/main.nf'
include { get_novel_stats; combine_stats_novel } from './modules/compare/main.nf'
include { parse_stats } from './modules/compare/main.nf'

vcf_1 = params.vcf_1    
vcf_2 = params.vcf_2

dbsnp = params.dbsnp

sample_list = file(params.sample_list)

vcf_1_ch = Channel.fromFilePairs("${vcf_1}*", checkIfExists:true) { file -> file.getSimpleName() }
vcf_2_ch = Channel.fromFilePairs("${vcf_2}*", checkIfExists:true) { file -> file.getSimpleName() }
dbsnp_ch = Channel.fromFilePairs("${dbsnp}*", checkIfExists:true) { file -> file.getSimpleName() }

workflow {
 
    vcf_1_ch_ss = select_samples(vcf_1_ch,sample_list)
    vcf_2_ch_ss = select_samples_2(vcf_2_ch,sample_list)

    stats_main = get_stats(vcf_1_ch_ss, vcf_2_ch_ss, sample_list)

    unique_sites_and_intersections = get_unique_sites_and_overlaps(vcf_1_ch_ss,vcf_2_ch_ss,sample_list)  
    unique_sites_and_intersections_indexed = unique_sites_and_intersections | flatten | get_unique_sites_and_overlap_index
    
    (stats_novel,vcf) = get_novel_stats(dbsnp_ch, unique_sites_and_intersections_indexed)

    stats_novel_combined = combine_stats_novel(stats_novel.collect(),vcf)
   
    parse_stats(stats_main, stats_novel_combined)
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
