#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and gVCF file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gvcf_file = file(row['gVCF'])
            return [ sample_id, gvcf_file ]
        }
        .into{samples_1; samples_2}

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    
    input:
    set val(sample_id), file(gvcf_file) from samples_1
    
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tgVCF: ${gvcf_file}\n"
    """
}
*/

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVG" }
    echo true
    publishDir "${params.out_dir}/genome-calling", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process index_vcf {
    tag { "${sample_id}.iV" }
    publishDir "${params.out_dir}/index-vcf/${sample_id}", mode: 'move', overwrite: false
    label 'gatk'
    input:
    set val(sample_id), file(vcf_file) from samples_2
    output:
    file("*.tbi")
    script:
    """
    tabix -p vcf ${vcf_file}
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
