#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fastq_r1_file = file(row['FastqR1'])
            def fastq_r2_file = file(row['FastqR2'])
            def flowcell = row['Flowcell']
            def lane = row['Lane']
            return [ sample_id, fastq_r1_file, fastq_r2_file, flowcell, lane ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_dict = Channel.fromPath(params.ref_seq_dict).toList()
ref_seq_amb = Channel.fromPath(params.ref_seq_amb).toList()
ref_seq_ann = Channel.fromPath(params.ref_seq_ann).toList()
ref_seq_bwt = Channel.fromPath(params.ref_seq_bwt).toList()
ref_seq_pac = Channel.fromPath(params.ref_seq_pac).toList()
ref_seq_sa = Channel.fromPath(params.ref_seq_sa).toList()
if (params.build == "b38") {
  ref_seq_alt = Channel.fromPath(params.ref_seq_alt).toList()
}
known_indels_1 = Channel.fromPath(params.known_indels_1).toList()
known_indels_1_index = Channel.fromPath(params.known_indels_1_index).toList()
known_indels_2 = Channel.fromPath(params.known_indels_2).toList()
known_indels_2_index = Channel.fromPath(params.known_indels_2_index).toList()
dbsnp = Channel.fromPath(params.dbsnp).toList()
dbsnp_index = Channel.fromPath(params.dbsnp_index).toList()

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file), val(flowcell), val(lane) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\tFlowcell: ${flowcell}\tLane: ${lane}\n"
    """
}
*/

process log_tool_version_samtools_bwa {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.samtools.bwa.version") into tool_version_samtools_bwa

    script:
    """
    samtools --version > tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVG" }
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

if (params.build == "b37") {
  process run_bwa_build37 {
      tag { "${params.project_name}.${sample_id}.rBwa" }
      memory { 64.GB * task.attempt }
      cpus { "${params.bwa_threads}" }
      publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
      label 'bwa_samtools'
  
      input:
      set val(sample_id), file(fastq_r1_file), file(fastq_r2_file), val(flowcell), val(lane) from samples_2
      file (ref) from ref_seq
      file (ref_index) from ref_seq_index
      file (ref_amb) from ref_seq_amb
      file (ref_ann) from ref_seq_ann
      file (ref_bwt) from ref_seq_bwt
      file (ref_pac) from ref_seq_pac
      file (ref_sa) from ref_seq_sa
  
      output:
      set val("$sample_id"), file("${sample_id}.bam")  into raw_bam
      
      script:
     
      readgroup_info="@RG\\tID:$flowcell.$lane\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"
      
      if(lane == "0") {
          sample_id = "$sample_id"
      } else {
  	sample_id = "$sample_id-${flowcell}.${lane}"
      }
     
      nr_threads = task.cpus - 1
      
      """
      bwa mem \
      -R \"${readgroup_info}\" \
      -t ${nr_threads}  \
      -K 100000000 \
      -Y \
      ${ref} \
      ${fastq_r1_file} \
      ${fastq_r2_file} | \
      samtools sort \
      -@ ${nr_threads} \
      -m ${params.memory_per_thread} \
      - > ${sample_id}.bam
      """
  }
}

if (params.build == "b38") {
  process run_bwa_build38 {
      tag { "${params.project_name}.${sample_id}.rBwa" }
      memory { 64.GB * task.attempt }
      cpus { "${params.bwa_threads}" }
      publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
      label 'bwa_samtools'
  
      input:
      set val(sample_id), file(fastq_r1_file), file(fastq_r2_file), val(flowcell), val(lane) from samples_2
      file (ref) from ref_seq
      file (ref_index) from ref_seq_index
      file (ref_amb) from ref_seq_amb
      file (ref_ann) from ref_seq_ann
      file (ref_bwt) from ref_seq_bwt
      file (ref_pac) from ref_seq_pac
      file (ref_sa) from ref_seq_sa
      file (ref_alt) from ref_seq_alt
  
      output:
      set val("$sample_id"), file("${sample_id}.bam")  into raw_bam
      
      script:
     
      readgroup_info="@RG\\tID:$flowcell.$lane\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"
      
      if(lane == "0") {
          sample_id = "$sample_id"
      } else {
  	sample_id = "$sample_id-${flowcell}.${lane}"
      }
     
      nr_threads = task.cpus - 1
      
      """
      bwa mem \
      -R \"${readgroup_info}\" \
      -t ${nr_threads}  \
      -K 100000000 \
      -Y \
      ${ref} \
      ${fastq_r1_file} \
      ${fastq_r2_file} | \
      samtools sort \
      -@ ${nr_threads} \
      -m ${params.memory_per_thread} \
      - > ${sample_id}.bam
      """
  }
}
process run_mark_duplicates {
    tag { "${params.project_name}.${sample_id}.rMD" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file) from raw_bam

    output:
    set val(sample_id), file("${sample_id}.md.bam"), file("${sample_id}.md.bai")  into md_bam

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms4g -Xmx${mem}g"  \
    MarkDuplicates \
    --MAX_RECORDS_IN_RAM 5000 \
    --INPUT ${bam_file} \
    --METRICS_FILE ${bam_file}.metrics \
    --TMP_DIR . \
    --ASSUME_SORT_ORDER coordinate \
    --CREATE_INDEX true \
    --OUTPUT ${sample_id}.md.bam
    """
}

process run_create_recalibration_table {
    tag { "${params.project_name}.${sample_id}.rCRT" }
    memory { 16.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file), file(bam_file_index) from md_bam
    file (ref) from ref_seq
    file (ref_index) from ref_seq_index
    file (ref_dict) from ref_seq_dict
    file (dbsnp_file) from dbsnp
    file (dbsnp_index_file) from dbsnp_index
    file (known_indels_1_file) from known_indels_1
    file (known_indels_1_index_file) from known_indels_1_index
    file (known_indels_2_file) from known_indels_2
    file (known_indels_2_index_file) from known_indels_2_index

    output:
    set val(sample_id), file("${sample_id}.md.bam"), file("${sample_id}.md.bai"), file("${sample_id}.recal.table")  into recal_table

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    BaseRecalibrator \
    --input ${bam_file} \
    --output ${sample_id}.recal.table \
    --tmp-dir . \
    -R ${ref} \
    --known-sites ${dbsnp_file} \
    --known-sites ${known_indels_1_file} \
    --known-sites ${known_indels_2_file}
    """
}

process run_recalibrate_bam {
    tag { "${params.project_name}.${sample_id}.rRB" }
    memory { 16.GB * task.attempt }
    cpus { 2 }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'gatk'

    input:
    set val(sample_id), file(bam_file), file(bam_file_index), file(recal_table_file) from recal_table
    file (ref) from ref_seq
    file (ref_index) from ref_seq_index
    file (ref_dict) from ref_seq_dict

    output:
    set val(sample_id), file("${sample_id}.md.recal.bam")  into recal_bam
    set val(sample_id), file("${sample_id}.md.recal.bai")  into recal_bam_index

    script:
      mem = task.memory.toGiga() - 4
    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
     ApplyBQSR \
    --input ${bam_file} \
    --output ${sample_id}.md.recal.bam \
    --tmp-dir . \
    -R ${ref} \
    --create-output-bam-index true \
    --bqsr-recal-file ${recal_table_file}
    """
}

process bam_to_cram {
    tag { "${params.project_name}.${sample_id}.btC" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file) from recal_bam
    file (ref) from ref_seq

    output:
    set val(sample_id), file("${bam_file.baseName}.cram") into cram_file

    script:
    """
    samtools view \
    --reference ${ref} \
    --output-fmt cram,version=3.0 \
    -o ${bam_file.baseName}.cram  ${bam_file}
    """
}

cram_file.into{ cram_file_1; cram_file_2; cram_file_3}

process index_cram {
    tag { "${params.project_name}.${sample_id}.iC" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(cram_file) from cram_file_1

    output:
    set val(sample_id), file("${cram_file}.crai") into cram_index

    script:
    """
    samtools index  \
    ${cram_file} ${cram_file}.crai
    """
}

process run_cram_flagstat {
    tag { "${params.project_name}.${sample_id}.rCF" }
    memory { 4.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(cram_file) from cram_file_2

    output:
    set val(sample_id), file("${cram_file}.flagstat") into cram_stats

    script:
    """
    samtools flagstat \
    --threads ${params.bwa_threads} \
    ${cram_file} > ${cram_file}.flagstat  \
    """
}

cram_file_3.mix(cram_index,cram_stats).groupTuple().set{cram_all}

process create_cram_md5sum {
    tag { "${params.project_name}.${sample_id}.cCMD5" }
    memory { 4.GB * task.attempt }
    publishDir "${params.out_dir}/${sample_id}", mode: 'move', overwrite: false
    input:
    set val(sample_id), file(cram_file) from cram_all

    output:
    set val(sample_id), file(cram_file), file("${cram_file[0]}.md5"), file("${cram_file[1]}.md5") into cram_all_md5sum

    script:
    """
    md5sum ${cram_file[0]} > ${cram_file[0]}.md5
    md5sum ${cram_file[1]} > ${cram_file[1]}.md5
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
