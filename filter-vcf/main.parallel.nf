#!/usr/bin/env nextflow

in_files = params.in_files
out_dir = file(params.out_dir)

// Problems with VQSR X, Y and MT only, only doing on chr 1 -> 22

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
    publishDir "${out_dir}/${params.cohort_id}/filter-vcf", mode: 'copy', overwrite: false
    label 'bcftools'
    time = 336.h    
   
    input:
      set val (file_name), file (vcf) from vcfs1
      each chr from chroms

    output:
    file("${filebase}.${chr}.filter-pass.vcf.gz") into vcf_pass_out
    file("${filebase}.${chr}.filter-pass.vcf.gz.tbi") into vcf_pass_out_index

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
    publishDir "${out_dir}/${params.cohort_id}/filter-vcf", mode: 'copy', overwrite: false
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

vcf_pass_out.toList().set{ concat_ready  }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/filter-vcf", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.filter-pass.vcf.gz"), file("${params.cohort_id}.filter-pass.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.1\\.filter-pass.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.2\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.3\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.4\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.5\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.6\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.7\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.8\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.9\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.10\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.11\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.12\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.13\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.14\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.15\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.16\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.17\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.18\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.19\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.20\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.21\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.22\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.X\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.Y\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.MT\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.filter-pass.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.filter-pass.vcf.gz
       """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build38 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/filter-vcf", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.filter-pass.vcf.gz"), file("${params.cohort_id}.filter-pass.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.chr1\\.filter-pass.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr2\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr3\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr4\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr5\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr6\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr7\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr8\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr9\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr10\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr11\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr12\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr13\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr14\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr15\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr16\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr17\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr18\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr19\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr20\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr21\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr22\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrX\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrY\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrM\\.filter-pass.vcf.gz" >> ${params.cohort_id}.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.filter-pass.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.filter-pass.vcf.gz
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
