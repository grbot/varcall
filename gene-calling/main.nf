#!/usr/bin/env nextflow

"""
Author: Mamana M.
Affiliation: University of Cape Town

Latest modification: Gerrit divided gene and genome calling into two separate nextflow workflows.
"""

chromosomes = params.chromosomes.split(',')

Channel.from( file(params.gvcf_file) )
        .set{ gvcf_file_cha }

process split_bed_to_chr {
    tag { "${params.cohort_id}.${chr}.sBtC" }
    label "bigmem"
    input:
        set file(gvcf_file), file(gene_region_bed) from gvcf_file_cha.combine([file(params.gene_region_bed)])
        each chr from chromosomes
    output:
        set chr, file(gvcf_file), file(outBED) into split_bed_to_chr
    script:
        inBED = gene_region_bed
        outBED = "${gene_region_bed.baseName}_chr${chr}.bed"
        template "split_bed_to_chr.py"
}


process run_genotype_gvcf_on_genes {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rGGoG" }
    label "bigmem"
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
        set chr, file(gvcf_file), file(bed_file) from split_bed_to_chr
    output:
        set chr, file(vcf_out), file(vcf_index_out) into vcf
    script:
         call_conf = 30 // set default
         if ( params.sample_coverage == "high" )
           call_conf = 30
         else if ( params.sample_coverage == "low" )
           call_conf = 10
        base = file(file(gvcf_file.baseName).baseName).baseName
        vcf_out = "${base}_chr${chr}_genes.vcf.gz"
        vcf_index_out = "${base}_chr${chr}_genes.vcf.gz.tbi"
        """
        ${params.tabix_base}/tabix ${gvcf_file}
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${bed_file} \
            -V ${gvcf_file} \
            -stand-call-conf ${call_conf} \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O ${vcf_out}
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
