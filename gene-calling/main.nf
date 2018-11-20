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
        set chrm, file(vcf_out) into GenotypeGVCF
    script:
        base = file(file(gVCF_file.baseName).baseName).baseName
        vcf_out = "${base}_chr${chrm}.vcf.gz"
        """
        tabix ${gVCF_file}
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${chrm} \
            -V ${gVCF_file} \
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
        set chrm, file(vcf_out) into GenotypeGVCF_genes
    script:
        base = file(file(gVCF_file.baseName).baseName).baseName
        vcf_out = "${base}_chr${chrm}_genes.vcf.gz"
        """
        tabix ${gVCF_file}
        ${params.gatk_base}/gatk \
            GenotypeGVCFs \
            -R ${params.ref_seq} \
            -L ${params.gene_region_bed} \
            -V ${gVCF_file} \
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
