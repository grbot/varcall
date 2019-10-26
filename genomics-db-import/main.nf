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

db = file(params.db_path)

if ( (db_update == "yes") && db.exists()   ){
  println "Db:${db} exists. Making a copy before update."

  process makeCopy {
    publishDir "${baseDir}", mode: 'move', overwrite: false
    input:
    file db

    output:
    file("${db.getName()}.backup") into db_backup

    script:
    """
    cp -L ${db} ${db.getName()}.backup
    """
  }
} else if ( (db_update == "yes") && !db.exists() ){
  println "Did not find existing db:${db} please specify the correct path."
  exit 1
} else if ( db_update == "no" ){
  println "Creating new db:${db}."
}

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

process create_variant_list {
    tag { "${params.project_name}.${params.cohort_id}.cVL" }
    publishDir "${params.out_dir}/${params.cohort_id}/genomics-db-import", mode: 'copy', overwrite: false

    input:
    val gvcf_file from gvcf_files

    output:
    file("gvcf.list") into gvcf_list

    script:
    """
    echo "${gvcf_file.join('\n')}" > gvcf.list
    """
}

process run_genomics_db_import {
    tag { "${params.project_name}.${params.cohort_id}.${chr}.rGDI" }
    label 'bigmem'
    publishDir "${params.out_dir}/${params.cohort_id}", mode: 'copy', overwrite: false

    input:
    file(gvcf_list)
    each chr from chroms
    file db
    val db_update

    output:
    file("${params.cohort_id}.${chr}.gdb")  into cohort_chr_calls

    script:
    if (db_update == "no"){
      """
      ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
      GenomicsDBImport \
      --L ${chr} \
      --variant ${gvcf_list} \
      --genomicsdb-workspace-path ${db}/${chr}.gdb
      """
    }else if (db_update == "yes"){
    ${params.gatk_base}/gatk --java-options "-Xmx${task.memory.toGiga()}g"  \
    GenomicsDBImport \
    --L ${chr} \
    --variant ${gvcf_list} \
    --genomicsdb-update-workspace-path ${db}/${chr}.gdb
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
