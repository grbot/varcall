#!/usr/bin/env nextflow

// Adapted from https://github.com/mamanambiya/phasing_vcf_eagle/ , recognision to Mamana Mbiyavanga!

if (params.build == "b37") {
  //chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT".split(',')
  chroms = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
} else if (params.build == "b38"){
    //chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM".split(',')
    chroms = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22".split(',')
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

eagle_genetic_map = Channel.fromPath(params.eagle_genetic_map).toList()

phasing_method = params.phasing_method
minRatio = params.minRatio

Channel.from( file(params.vcf) )
        .set{ vcf_ch }
Channel.from( file(params.vcf_index) )
        .set{ vcf_index_ch }
/*
Channel
  .from('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')
  .map { chr -> tuple("chr$chr",file("/scratch3/users/gerrit/projects/adme/1000GP_Phase3_b38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz"), file("/scratch3/users/gerrit/projects/adme/1000GP_Phase3_b38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz.tbi")) }
  .set { ref_panel }
*/

process log_tool_version_tabix {
    tag { "${params.project_name}.ltVT" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    label 'gatk'
   
    output:
    file("tool.tabix.version") into tool_version_tabix

    script:
    """
    tabix 2>&1 | grep Version  > tool.tabix.version
    """
}

process log_tool_version_bcftools {
    tag { "${params.project_name}.ltVB" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    label 'impute'

    output:
    file("tool.bcftools.version") into tool_version_bcftools

    script:
    """
    bcftools -v > tool.bcftools.version 2>&1
    """
}

process log_tool_version_eagle2 {
    tag { "${params.project_name}.ltVE" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    label 'impute'

    output:
    file("tool.eagle2.version") into tool_version_eagle2

    script:
    """
    eagle -h > tool.eagle2.version 2>&1
    """
}

process split_vcf_to_chrm {
    tag { "${params.project_name}.${chr}.svtC" }
    echo true
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    memory { 16.GB * task.attempt }
    cpus {16}
    label 'impute'

    input:
    file vcf from vcf_ch
    file vcf_index from vcf_index_ch
    each chr from chroms
   
    output:
    set chr, file("${params.cohort_id}.phasing-ready.${chr}.vcf.gz"), file("${params.cohort_id}.phasing-ready.${chr}.vcf.gz.tbi") into phasing_ready_vcf
    
    script:
    """
    bcftools norm \
        --regions ${chr} \
        -m -both \
        --threads ${task.cpus} \
        -Oz \
        -o ${params.cohort_id}.phasing-ready.${chr}.vcf.gz \
        ${vcf}
    bcftools index --threads ${task.cpus} -t ${params.cohort_id}.phasing-ready.${chr}.vcf.gz
    """
}

// ref_panel.mix(phasing_ready_vcf).groupTuple().set{cvcfs}

/*
process phase_vcf_with_reference {
    tag { "${params.project_name}.${chr}.pV" }
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    memory { 64.GB * task.attempt }    
    label 'impute'
    
    input:
    set chr, vcf, vcf_index from phasing_ready_vcf
    file ref_index from "/scratch3/users/gerrit/projects/adme/1000GP_Phase3_b38/ALL.*_GRCh38.genotypes.20170504.bcf.csi"
    file eagle_genetic_map from eagle_genetic_map
    
    output:
    file("${params.cohort_id}.phased.${chr}.vcf.gz") into phased_vcf
    file("${params.cohort_id}.phased.${chr}.vcf.gz.tbi") into phased_vcf_index
    
    script:
        """
        eagle \
            --numThreads=${task.cpus} \
            --vcfTarget=${vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfRef=/scratch3/users/gerrit/projects/adme/1000GP_Phase3_b38/ALL.${chr}_GRCh38.genotypes.20170504.bcf \
            --vcfOutFormat=z \
            --chrom=${chr} \
            --outPrefix=${params.cohort_id}.phased.${chr} 2>&1 | tee ${params.cohort_id}.phased.${chr}.vcf.gz.log
        bcftools index -t ${params.cohort_id}.phased.${chr}.vcf.gz
        """
}
*/

process phase_vcf_with_out_reference {
    tag { "${params.project_name}.${chr}.pV" }
    publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
    memory { 64.GB * task.attempt }    
    cpus {16}  
    label 'impute'
    
    input:
    set chr, vcf, vcf_index from phasing_ready_vcf
    file eagle_genetic_map from eagle_genetic_map
    
    output:
    file("${params.cohort_id}.phased.${chr}.vcf.gz") into phased_vcf
    file("${params.cohort_id}.phased.${chr}.vcf.gz.tbi") into phased_vcf_index
    
    script:
        """
        eagle \
            --numThreads=${task.cpus} \
            --vcf=${vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfOutFormat=z \
            --chrom=${chr} \
            --outPrefix=${params.cohort_id}.phased.${chr} 2>&1 | tee ${params.cohort_id}.phased.${chr}.vcf.gz.log
        bcftools index -t ${params.cohort_id}.phased.${chr}.vcf.gz
        """
}

phased_vcf.toList().set{ concat_ready }

if (params.build == "b37") {
  process run_concat_vcf_build37 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.phased.vcf.gz"), file("${params.cohort_id}.phased.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.1\\.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.2\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.3\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.4\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.5\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.6\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.7\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.8\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.9\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.10\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.11\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.12\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.13\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.14\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.15\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.16\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.17\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.18\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.19\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.20\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.21\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.22\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.Y\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.MT\\.g.vcf.gz" >> ${params.cohort_id}.vcf.list 
       
      gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.phased.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.phased.vf.va.vcf.gz
       """
  }
}

if (params.build == "b38") {
  process run_concat_vcf_build38 {
       tag { "${params.project_name}.${params.cohort_id}.rCV" }
       memory { 16.GB * task.attempt }  
       publishDir "${params.out_dir}/${params.cohort_id}/phasing", mode: 'copy', overwrite: false
       label 'gatk'
  
       input:
       file(vcf) from concat_ready
  
       output:
  	   set file("${params.cohort_id}.phased.vcf.gz"), file("${params.cohort_id}.phased.vcf.gz.tbi") into combined_calls
  
       script:
         mem = task.memory.toGiga() - 4
       """
       echo "${vcf.join('\n')}" | grep "\\.chr1\\.vcf.gz" > ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr2\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr3\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr4\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr5\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr6\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr7\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr8\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr9\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr10\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr11\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr12\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr13\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr14\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr15\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr16\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr17\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr18\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr19\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr20\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr21\\.vcf.gz" >> ${params.cohort_id}.vcf.list
       echo "${vcf.join('\n')}" | grep "\\.chr22\\.vcf.gz" >> ${params.cohort_id}.vcf.list
  
       gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
       GatherVcfs \
       -I ${params.cohort_id}.vcf.list \
       -O ${params.cohort_id}.phased.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
       tabix -p vcf ${params.cohort_id}.phased.vcf.gz
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
