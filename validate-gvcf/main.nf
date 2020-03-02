#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gender = row['Gender']
            def gvcf_file = file(row['gVCF'])
            return [ sample_id, gender, gvcf_file ]
        }.into{samples_1; samples_2; samples_3}

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

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_dict = Channel.fromPath(params.ref_seq_dict).toList()

chrom_list = Arrays.asList(chroms)

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), val(gender), file(gvcf_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tGender: ${gender}\tgVCF: ${gvcf_file}\n"
    """
}
*/

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

process create_interval_list_male {
    tag { "${params.project_name}.cILM" }
    publishDir "${params.out_dir}/validate-gvcf", mode: 'symlink', overwrite: false

    input:
    chrom_list

    output:
    file("interval.list") into intervals_male

    script:
    """
      echo "${chrom_list.join('\n')}" > interval.list
    """
}

process create_interval_list_female {
    tag { "${params.project_name}.cILF" }
    publishDir "${params.out_dir}/validate-gvcf", mode: 'symlink', overwrite: false

    input:
    chrom_list

    output:
    file("interval.list") into intervals_female

    script:
    """
      echo "${chrom_list.join('\n')}" > tmp.list
      grep -v "Y" tmp.list > interval.list
    """
}


// Now do Male and Female validation
samples_2.filter{it[1] == 'M'}.set{samples_male}
samples_3.filter{it[1] == 'F'}.set{samples_female}

// Males
process validate_male {
     tag { "${params.project_name}.${sample_id}.vM" }
     memory { 16.GB * task.attempt }
     label 'gatk'
     publishDir "${params.out_dir}/validate-gvcf/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(gvcf_file) from samples_male
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (intervals_male)

     output:
          file("${sample_id}.validatevariants") into validatevariants_male_file

     script:
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-Xmx${mem}g" \
     ValidateVariants \
     -R $ref_seq \
     -L $intervals_male \
     "--validate-GVCF" \
     -V $gvcf_file \
     > ${sample_id}.validatevariants 2>&1
     """
}

process validate_female {
     tag { "${params.project_name}.${sample_id}.vM" }
     memory { 16.GB * task.attempt }
     label 'gatk'
     publishDir "${params.out_dir}/validate-gvcf/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(gvcf_file) from samples_female
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (intervals_female)

     output:
          file("${sample_id}.validatevariants") into validatevariants_female_file

     script:
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-Xmx${mem}g" \
     ValidateVariants \
     -R $ref_seq \
     -L $intervals_female \
     "--validate-GVCF" \
     -V $gvcf_file \
     > ${sample_id}.validatevariants 2>&1
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
