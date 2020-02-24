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
db_update = params.db_update

if ( (db_update == "yes") && db.exists() ){
  println "Db:${db} exists. Making a copy before update."

  process make_copy {
    publishDir "${params.out_dir}/genomics-db-import", mode: 'move', overwrite: false
    input:
    file db

    output:
    file("${db.getName()}.backup") into db_backup

    script:
    """
    cp -rL ${db} ${db.getName()}.backup
    """
  }
} else if ( (db_update == "yes") && !db.exists() ){
    println "Did not find existing db:${db} please specify the correct path."
    exit 1
} else if ( ( db_update == "no" ) && !db.exists() ){
    println "Creating new db:${db}"
    result = db.mkdir()
    println result ? "OK" : "Cannot create directory: $db"
} else if ( ( db_update == "no" ) && db.exists() ){
    println "Db exists:${db} . Please check if you want to create new. If so please delete first."
    exit 1
}


if (params.build == "b37") {
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
} else if (params.build == "b38"){
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

chrom_list = Arrays.asList(chroms)

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

process create_interval_list {
    tag { "${params.project_name}.cIL" }
    publishDir "${params.out_dir}/genomics-db-import", mode: 'symlink', overwrite: false

    input:
    chrom_list

    output:
    file("interval.list") into intervals

    script:
    """
      echo "${chrom_list.join('\n')}" > interval.list
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

if (db_update == "no"){
  process run_genomics_db_import_new {
      tag { "${params.project_name}.rGDIN" }
      memory { 232.GB * task.attempt }  
      publishDir "${params.out_dir}/genomics-db-import", mode: 'symlink', overwrite: false
      label 'gatk'
  
      input:
      file(gvcf_list)
      file db
      file intervals
  
      output:
      file("${db}.gdb")  into cohort_chr_calls
  
      script:
      mem = task.memory.toGiga() - 32
        // We delete the database first if it was created on a failed run. E.g. when memory was to low on previous nextflow process.
      """
      mkdir tmp && \
      rm -rf ${db}.gdb && \
      gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g"  \
      GenomicsDBImport \
      --intervals ${intervals} \
      --variant ${gvcf_list} \
      --genomicsdb-workspace-path ${db}.gdb \
      --batch-size 50 \
      --reader-threads 5 \
     --tmp-dir tmp
      """
  }
}

if (db_update == "yes"){
  process run_genomics_db_import_update {
      tag { "${params.project_name}.rGDIU" }
      memory { 232.GB * task.attempt }  
      publishDir "${params.out_dir}/genomics-db-import", mode: 'symlink', overwrite: false
      label 'gatk'
  
      input:
      file(gvcf_list)
      file db
      file db_backup // this staging forces the process to wait for the backup to be made before continuing
      file intervals 
  
      output:
      file("${db}.gdb")  into cohort_chr_calls
  
      script:
      mem = task.memory.toGiga() - 10
      """
      mkdir tmp && \
      gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g"  \
      GenomicsDBImport \
      --intervals ${intervals} \
      --variant ${gvcf_list} \
      --genomicsdb-update-workspace-path ${db}.gdb \
      --batch-size 50 \
      --reader-threads 5 \
      --tmp-dir tmp
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
