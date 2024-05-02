#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref                       = file(params.ref, type: 'file')
db                        = file(params.db_path, type: 'dir')
db_update                 = params.db_update
project_name              = params.project_name
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_genomics_db_import_new {
    tag { "${project_name}.${chr}.rGDIN" }
    label 'gatk'
    memory { 250.GB * task.attempt }
    time '240h'
    cpus 2
    errorStrategy 'finish'
    publishDir "${db}", mode: 'copy', overwrite: true
    
    input:
    path(gvcf_list)
    each chr

    output:
    path("${chr.replaceAll(":","_")}.gdb"), emit: interval_db
    
    script:
    mem = task.memory.toGiga() - ( task.memory.toGiga() * 1/5 )
    rthreads = task.cpus - 1
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        GenomicsDBImport \
        --reference ${ref} \
        --intervals ${chr} \
        --sample-name-map ${gvcf_list} \
        --batch-size 50 \
        --bypass-feature-reader \
        --consolidate \
        --genomicsdb-shared-posixfs-optimizations \
        --genomicsdb-workspace-path ${chr.replaceAll(":","_")}.gdb
    """
}

// --reader-threads ${rthreads} \

// I'VE COMMENTED THE ACTUAL DB BACKUP COMMAND - DB T00 LARGE
process run_backup_genomic_db {
    tag { "DB_Backup" }
    output:
    path("backup_status.txt"), emit: backup_status
    
    """
    backup=`date +"%Y_%m_%d-%H_%M_%S"`
    ## rsync -avhP ${db} ${db}.\$backup
    echo -e "Done backing up ${db}." | tee backup_status.txt
    """
}

process run_genomics_db_import_update {
    tag { "${params.project_name}.${chr}.rGDIU" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    time '168h'
    publishDir "${db}", mode: 'copy', overwrite: true
  
    input:
    path(gvcf_list)
    each interval
    path(backup_status)

    output:
    path("${chr.replaceAll(":","_")}.gdb"), emit: interval_db
    
    script:
    mem = task.memory.toGiga() - 2
    """
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenomicsDBImport \
        --L ${chr} \
        --variant ${gvcf_list} \
        --batch-size 50 \
        --reader-threads 5 \
        --genomicsdb-update-workspace-path ${chr.replaceAll(":","_")}.gdb
    """
}