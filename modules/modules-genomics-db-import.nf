#!/usr/bin/env nextflow

db                        = file(params.db_path, type: 'dir')
db_update                 = params.db_update
project_name              = params.project_name
cohort_id                 = params.cohort_id
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_genomics_db_import_new {
    tag { "${project_name}.${chr}.rGDIN" }
    label 'gatk'
    memory { 230.GB * task.attempt }  
    publishDir "${outdir}/genomics-db-import", mode: 'symlink', overwrite: false
  
    input:
    path(gvcf_list)
    each chr
  
    script:
    mem = task.memory.toGiga() - 32

    """
    mkdir -p ${db}
    rm -rf ${db}/${chr}.gdb

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenomicsDBImport \
        --L ${chr} \
        --variant ${gvcf_list} \
        --genomicsdb-workspace-path ${db}/${chr}.gdb
    """
}

process run_genomics_db_import_update {
    tag { "${params.project_name}.${chr}.rGDIU" }
    label 'gatk'
    memory { 24.GB * task.attempt }  
    publishDir "${outdir}/genomics-db-import", mode: 'copy', overwrite: false
  
    input:
    path(gvcf_list)
    each chr
  
    script:
    mem = task.memory.toGiga() - 10
    """
    ## backup=`date +"%Y_%m_%d-%H_%M_%S"`
    cp -rL ${db} ${db}.backup
    
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenomicsDBImport \
        --L ${chr} \
        --variant ${gvcf_list} \
        --genomicsdb-update-workspace-path ${db}/${chr}.gdb
    """
}
