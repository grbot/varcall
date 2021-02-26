#!/usr/bin/env nextflow

vcfs = Channel.fromPath( file(params.vcfs_path) )

process log_tool_version_gatk {
    tag { "ltVV" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process validate_vcf {
     tag { "${vcf}.vV" }
     memory { 16.GB * task.attempt }
     label 'gatk'
     publishDir "${params.out_dir}/validate-vcf/", mode: 'copy', overwrite: false

     input:
     file(vcf) from vcfs 

     output:
          file("${vcf}.validatevariants") into validatevariants_male_file

     script:
     mem = task.memory.toGiga() - 4
     """
     tabix -p vcf $vcf
     gatk --java-options "-Xmx${mem}g" \
     ValidateVariants \
     --validation-type-to-exclude ALL \
     -V $vcf \
     > ${vcf}.validatevariants 2>&1
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
