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

process create_variant_list {
    tag { "${params.project_name}.cVL" }
    publishDir "${params.out_dir}/combine-gvcfs", mode: 'copy', overwrite: false

    input:
    val gvcf_file from gvcf_files
    
    output:
    file("gvcf.list") into gvcf_list
    
    script:
    """
    echo "${gvcf_file.join('\n')}" > gvcf.list
    """
}

process run_combine_gvcfs {
    tag { "${params.project_name}.${chr}.rCG" }
    memory { 16.GB * task.attempt } 
    publishDir "${params.out_dir}/combine-gvcfs", mode: 'copy', overwrite: false
    label 'gatk'
    
    input:
    file(gvcf_list)
    file ref_seq
    file ref_seq_index
    file ref_seq_dict
    each chr from chroms

    output:
    file("${params.cohort_id}.${chr}.g.vcf.gz")  into cohort_chr_calls
    file("${params.cohort_id}.${chr}.g.vcf.gz.tbi") into cohort_chr_indexes

    script:
    mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g"  \
    CombineGVCFs \
    -R ${params.ref_seq} \
    --L ${chr} \
    --variant ${gvcf_list} \
    -O ${params.cohort_id}.${chr}.g.vcf.gz
    """
}

cohort_chr_calls.toList().set{ cohort_calls }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/combine-gvcfs", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from cohort_calls 
  
       output:
  	   set file("${params.cohort_id}.g.vcf.gz"), file("${params.cohort_id}.g.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.Y\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.MT\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.g.vcf.gz # GatherVCF does not index the gVCF. The gVCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.g.vcf.gz
       """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build38 {
       tag { "${params.project_name}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/combine-gvcfs", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from cohort_calls
  
       output:
  	   set file("${params.cohort_id}.g.vcf.gz"), file("${params.cohort_id}.g.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.chr1\\.g.vcf.gz" > ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr2\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr3\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr4\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr5\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr6\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr7\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr8\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr9\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr10\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr11\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr12\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr13\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr14\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr15\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr16\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr17\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr18\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr19\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr20\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr21\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr22\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrX\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrY\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chrM\\.g.vcf.gz" >> ${params.cohort_id}.g.vcf.list
       
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.g.vcf.list \
       -O ${params.cohort_id}.g.vcf.gz # GatherVCF does not index the gVCF. The gVCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.g.vcf.gz
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
