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
        .tap{samples_1; samples_2}
        .map { sample_id, gvcf_file ->
            return [ gvcf_file ]
        }
        .collect().set { gvcf_files }


// Build is in the intervals file
Channel.fromPath( file(params.intervals) )
        .splitText()
        .map{it -> it.trim()}
        .set {intervals}

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
    publishDir "${params.out_dir}/genomics-db-import", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    mem = task.memory.toGiga() - 3
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

process create_variant_list {
    tag { "${params.project_name}.cVL" }
    publishDir "${params.out_dir}/genomics-db-import", mode: 'symlink', overwrite: false

    input:
    val gvcf_file from gvcf_files

    output:
    file("gvcf.list") into gvcf_list

    script:
    """
    echo "${gvcf_file.join('\n')}" > gvcf.list
    """
}

process run_genomics_db_import_new {
    tag { "${params.project_name}.${interval}.rGDIN" }
    memory { 10.GB }  
    cpus { 5 }
    publishDir "${params.out_dir}/genomics-db-import", mode: 'copy', overwrite: false
    label 'gatk'
  
    input:
    file(gvcf_list)
    each interval from intervals 
  
    output:
    file("*.gdb")  into cohort_chr_calls
  
    script:
    mem = task.memory.toGiga() - 2
    // We delete the database first if it was created on a failed run. E.g. when memory was to low on previous nextflow process.
    """
    a=$interval
    b=\${a/:/_}
    rm -rf ${interval}.gdb && \
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g"  \
    GenomicsDBImport \
    --L ${interval} \
    --variant ${gvcf_list} \
    --batch-size 50 \
    --reader-threads 5 \
    --genomicsdb-workspace-path \$b.gdb
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
