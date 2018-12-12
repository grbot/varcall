#!/usr/bin/env nextflow

"""
Author: Mamana M.
Affiliation: University of Cape Town
Aim: A simple.
Date: Mon June 
Run:

Latest modification:
- TODO
"""
chromosomes = params.chromosomes.split(',')

println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Chromosomes used: ${chromosomes.join(',')}"


Channel.from( file(params.gVCF_file) )
        .set{ gVCF_file_cha }
gVCF_file_cha.into{ gVCF_file_cha;  gVCF_file_cha1 }

"""
Step 1
"""
process GenotypeGVCF {
    tag { "GenotypeGVCF_${params.project_name}_${chrm}" }
    label "bigmem"
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
        file(gVCF_file) from gVCF_file_cha1
        each chrm from chromosomes
    output:
        set chrm, file(vcf_out), file(vcf_index_out)  into GenotypeGVCF
    script:
         call_conf = 30 // set default
         if ( params.sample_coverage == "high" )
           call_conf = 30
         else if ( params.sample_coverage == "low" )
           call_conf = 10
        base = file(file(gVCF_file.baseName).baseName).baseName
        vcf_out = "${base}_chr${chrm}.vcf.gz"
        vcf_index_out = "${base}_chr${chrm}.vcf.gz.tbi"
        """
        ${params.tabix_base}/tabix ${gVCF_file}
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${chrm} \
            -V ${gVCF_file} \
            -stand-call-conf ${call_conf} \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O ${vcf_out}
        """
}


process split_bed_to_chr {
    tag { "split_bed_to_chr_${params.project_name}_${chrm}" }
    label "bigmem"
    input:
        set file(gVCF_file), file(gene_region_bed) from gVCF_file_cha.combine([file(params.gene_region_bed)])
        each chrm from chromosomes
    output:
        set chrm, file(gVCF_file), file(outBED) into split_bed_to_chr
    script:
        inBED = gene_region_bed
        outBED = "${gene_region_bed.baseName}_chr${chrm}.bed"
        template "split_bed_to_chrm.py"
}


process GenotypeGVCF_genes {
    tag { "GenotypeGVCF_gene_${params.project_name}_${chrm}" }
    label "bigmem"
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    input:
        set chrm, file(gVCF_file), file(bed_file) from split_bed_to_chr
    output:
        set chrm, file(vcf_out), file(vcf_index_out) into GenotypeGVCF_genes
    script:
         call_conf = 30 // set default
         if ( params.sample_coverage == "high" )
           call_conf = 30
         else if ( params.sample_coverage == "low" )
           call_conf = 10
        base = file(file(gVCF_file.baseName).baseName).baseName
        vcf_out = "${base}_chr${chrm}_genes.vcf.gz"
        vcf_index_out = "${base}_chr${chrm}_genes.vcf.gz.tbi"
        """
        ${params.tabix_base}/tabix ${gVCF_file}
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${bed_file} \
            -V ${gVCF_file} \
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
