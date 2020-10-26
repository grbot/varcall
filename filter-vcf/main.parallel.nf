#!/usr/bin/env nextflow

in_files = params.in_files
out_dir = file(params.out_dir)

// Problems with VQSR X, Y and MT only, only doing on chr 1 -> 22

if (params.build == "b37") {
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
} else if (params.build == "b38"){
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22".split(',')
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

Channel.fromFilePairs(in_files)
        { file ->
          b = file.baseName
          m = b =~ /(.*)\.vcf.*/
          return m[0][1]
        }.into { vcfs1 ; vcfs2 }

process log_tool_version_bcftools {
    tag { "${params.project_name}.ltViB" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/filter-vcf", mode: 'move', overwrite: false
    label 'bcftools'

    output:
    file("tool.bcftools.version") into tool_version_bcftools

    script:
    mem = task.memory.toGiga() - 3
    """
    bcftools --version > tool.bcftools.version 2>&1
    """
}

process filter_short_snps_indels {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.fP" }
    publishDir "${out_dir}/${params.cohort_id}/filter-vcf", mode: 'move', overwrite: false
    label 'bcftools'
    time = 336.h    
   
    input:
      set val (file_name), file (vcf) from vcfs1
      each chr from chroms

    output:
    set file("${filebase}.${chr}.filter-pass.vcf.gz"), file("${filebase}.${chr}.filter-pass.vcf.gz.tbi") into vcf_pass_out

    script:
    filebase = (file(vcf[0].baseName)).baseName
    """
    bcftools view \
    --include "FILTER='PASS'" \
    -r ${chr} \
    -O z \
    -o "${filebase}.${chr}.filter-pass.vcf.gz" \
    ${vcf[0]} 
    bcftools index \
    -t \
    "${filebase}.${chr}.filter-pass.vcf.gz"
    """
}

process filter_other {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.fO" }
    publishDir "${out_dir}/${params.cohort_id}/filter-vcf", mode: 'move', overwrite: false
    label 'bcftools'
    time = 336.h

    input:
      set val (file_name), file (vcf) from vcfs2
      each chr from chroms
    output:
    set file("${filebase}.${chr}.filter-other.vcf.gz"), file("${filebase}.${chr}.filter-other.vcf.gz.tbi") into vcf_other_out

    script:
    filebase = (file(vcf[0].baseName)).baseName
    """
    bcftools view \
    --include "FILTER='.'" \
    -r ${chr} \
    -O z \
    -o "${filebase}.${chr}.filter-other.vcf.gz" \
    ${vcf[0]}
    bcftools index \
    -t \
    "${filebase}.${chr}.filter-other.vcf.gz"
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
