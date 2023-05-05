#!/usr/bin/env nextflow

db                        = file(params.db_path, type: 'dir')
db_update                 = params.db_update
project_name              = params.project_name
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

// THIS IS ONLY FOR RUNNING ON THE WITS CLUSTER =======================================================
if (params.build == "b37") {
    chroms = (1..22).toList() +  ["X","Y","MT"]
} else if (params.build == "b38") {
    chroms = (1..22).toList().collect { 'chr' + "${it}" } + ["chrX", "chrY","chrM"]
}
machines                  = ([2,5,10,19,24,25,26,27]+(29..45)).collect { it.intValue() }
def chromAlloc = [:];
chroms.eachWithIndex { chrom, i -> chromAlloc[chrom] = String.format("-w n%02d",machines[i]) }
println chromAlloc
// ====================================================================================================

process run_genomics_db_import_new {
    tag { "${project_name}.${chr}.rGDIN" }
    label 'gatk'
    memory { 60.GB * task.attempt }
    time '120h'
    clusterOptions { chromAlloc["${chr}"] }
    // publishDir "${outdir}/genomics-db-import", mode: 'symlink', overwrite: false
  
    input:
    path(gvcf_list)
    each chr
    
    script:
    mem = task.memory.toGiga() - 10

    """
    mkdir -p ${db}
    rm -rf ${db}/${chr}.gdb

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenomicsDBImport \
        --batch-size 50 \
        --L ${chr} \
        --variant ${gvcf_list} \
        --genomicsdb-workspace-path ${db}/${chr}.gdb
    """
}

process run_backup_genomic_db {
    tag { "DB_Backup" }
    output:
    path("backup_status.txt"), emit: backup_status
    
    """
    backup=`date +"%Y_%m_%d-%H_%M_%S"`
    rsync -avhP ${db} ${db}.\$backup
    echo -e "Done backing up ${db}." | tee backup_status.txt
    """
}

process run_genomics_db_import_update {
    tag { "${params.project_name}.${chr}.rGDIU" }
    label 'gatk'
    memory { 24.GB * task.attempt }
    time '120h'
    clusterOptions { chromAlloc["${chr}"] }
    // publishDir "${outdir}/genomics-db-import", mode: 'copy', overwrite: false
  
    input:
    path(gvcf_list)
    each chr
    path(backup_status)
  
    script:
    mem = task.memory.toGiga() - 10
    """
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenomicsDBImport \
        --batch-size 50 \
        --L ${chr} \
        --variant ${gvcf_list} \
        --genomicsdb-update-workspace-path ${db}/${chr}.gdb
    """
}
