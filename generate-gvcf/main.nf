#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def gender = row['Gender']
            def bam_file = file(row['BAM'])
            return [ sample_id, gender, bam_file ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}

if (params.build == "b37") {
  autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
  x = "X"
  y = "Y"
  // Coordinates are from here: https://www.ncbi.nlm.nih.gov/grc/human
  x_par1 = "60001-2699520"
  x_par2 = "154931044-155260560"
  y_par1 = "10001-2649520"
  y_par2 = "59034050-59363566"
  mt = "MT"
} else if (params.build == "b38"){
  autosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22".split(',')
  x = "chrX"
  y = "chrY"
  // No autosomal region changes from build 37 to 38
  x_par1 = "10001-2781479"
  x_par2 = "155701383-156030895"
  y_par1 = "10001-2781479"
  y_par2 = "56887903-57217415"
  mt = "chrM"
} else {
  println "\n============================================================================================="
  println "Please specify a genome build (b37 or b38)!"
  println "=============================================================================================\n"
  exit 1
}

ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_seq_index = Channel.fromPath(params.ref_seq_index).toList()
ref_seq_dict = Channel.fromPath(params.ref_seq_dict).toList()
dbsnp = Channel.fromPath(params.dbsnp).toList()
dbsnp_index = Channel.fromPath(params.dbsnp_index).toList()

/*
process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), val(gender), file(bam_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tGender: ${gender}\tBAM: ${bam_file}\n"
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

process run_haplotype_caller_on_autosomes {
    tag { "${params.project_name}.${sample_id}.${chr}.rHCoA" }
    memory { 8.GB * task.attempt }
    cpus { 2 }
    label 'gatk'
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), val(bam_file) from samples_2
    file (ref_seq)
    file (ref_seq_index)
    file (ref_seq_dict)
    file (dbsnp)
    file (dbsnp_index)
	  each chr from autosomes

    output:
	  set val(sample_id), file("${sample_id}.${chr}.g.vcf.gz")  into autosome_calls
	  set val(sample_id), file("${sample_id}.${chr}.g.vcf.gz.tbi") into autosome_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
    HaplotypeCaller \
    -R ${ref_seq} \
    -I $bam_file \
    --emit-ref-confidence GVCF \
    --dbsnp ${dbsnp} \
    --L $chr \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.${chr}.g.vcf.gz
    """
}

// Now do X and Y calling
samples_3.filter{it[1] == 'M'}.into{samples_male_1; samples_male_2; samples_male_3; samples_male_4; samples_male_5; samples_male_6}
samples_4.filter{it[1] == 'F'}.set{samples_female}

// Males
process run_haplotype_caller_on_x_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP1M" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_1
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.X_PAR1.g.vcf.gz") into x_par1_calls
	   set val(sample_id), file("${sample_id}.X_PAR1.g.vcf.gz.tbi") into x_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${x}:${x_par1} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR1.g.vcf.gz
     """
}

process run_haplotype_caller_on_x_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoXP2M" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_2
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.X_PAR2.g.vcf.gz") into x_par2_calls
	   set val(sample_id), file("${sample_id}.X_PAR2.g.vcf.gz.tbi") into x_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${x}:${x_par2} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X_PAR2.g.vcf.gz
     """
}

process run_haplotype_caller_on_x_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoXNPM" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_3
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.X_nonPAR.g.vcf.gz") into x_nonpar_calls
	   set val(sample_id), file("${sample_id}.X_nonPAR.g.vcf.gz.tbi") into x_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${x} -XL ${x}:${x_par1} -XL ${x}:${x_par2} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.X_nonPAR.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_par1_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP1M" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_4
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR1.g.vcf.gz") into y_par1_calls
	   set val(sample_id), file("${sample_id}.Y_PAR1.g.vcf.gz.tbi") into y_par1_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${y}:${y_par1} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR1.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_par2_male {
     tag { "${params.project_name}.${sample_id}.rHCoYP2M" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_5
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.Y_PAR2.g.vcf.gz") into y_par2_calls
	   set val(sample_id), file("${sample_id}.Y_PAR2.g.vcf.gz.tbi") into y_par2_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${y}:${y_par2} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.Y_PAR2.g.vcf.gz
     """
}

process run_haplotype_caller_on_y_nonpar_male {
     tag { "${params.project_name}.${sample_id}.rHCoYNPM" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_male_6
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.Y_nonPAR.g.vcf.gz") into y_nonpar_calls
	   set val(sample_id), file("${sample_id}.Y_nonPAR.g.vcf.gz.tbi") into y_nonpar_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${y} -XL ${y}:${y_par1} -XL ${y}:${y_par2} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 1 \
     -O ${sample_id}.Y_nonPAR.g.vcf.gz
    """
}

// Females
process run_haplotype_caller_on_x_female {
     tag { "${params.project_name}.${sample_id}.rHCoXF" }
     memory { 8.GB * task.attempt }
     cpus { 2 }
     label 'gatk'
     publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

     input:
     set val(sample_id), val(gender), val(bam_file) from samples_female
     file (ref_seq)
     file (ref_seq_index)
     file (ref_seq_dict)
     file (dbsnp)
     file (dbsnp_index)

     output:
	   set val(sample_id), file("${sample_id}.X.g.vcf.gz") into x_calls
	   set val(sample_id), file("${sample_id}.X.g.vcf.gz.tbi") into x_calls_indexes

     script:
     call_conf = 30 // set default
     if ( params.sample_coverage == "high" )
       call_conf = 30
     else if ( params.sample_coverage == "low" )
       call_conf = 10
     mem = task.memory.toGiga() - 4
     """
     gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
     HaplotypeCaller \
     -R ${ref_seq} \
     -I $bam_file \
     --emit-ref-confidence GVCF \
     --dbsnp ${dbsnp} \
     --L ${x} \
     --genotyping-mode DISCOVERY \
     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
     -stand-call-conf ${call_conf} \
     --sample-ploidy 2 \
     -O ${sample_id}.X.g.vcf.gz
     """

}

process run_haplotype_caller_on_mt {
    tag { "${params.project_name}.${sample_id}.rHCoMT" }
    memory { 8.GB * task.attempt }
    cpus { 2 }
    label 'gatk'
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false

    input:
    set val(sample_id), val(gender), file(bam_file) from samples_5
    file (ref_seq)
    file (ref_seq_index)
    file (ref_seq_dict)
    file (dbsnp)
    file (dbsnp_index)

    output:
	  set val(sample_id), file("${sample_id}.MT.g.vcf.gz") into mt_calls
	  set val(sample_id), file("${sample_id}.MT.g.vcf.gz.tbi") into mt_calls_indexes

    script:
    call_conf = 30 // set default
    if ( params.sample_coverage == "high" )
      call_conf = 30
    else if ( params.sample_coverage == "low" )
      call_conf = 10
    mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g"  \
    HaplotypeCaller \
    -R ${ref_seq} \
    -I $bam_file \
    --emit-ref-confidence GVCF \
    --dbsnp ${dbsnp} \
    --L ${mt} \
    --genotyping-mode DISCOVERY \
    -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    -stand-call-conf ${call_conf} \
    --sample-ploidy 2 \
    -O ${sample_id}.MT.g.vcf.gz
    """
}

autosome_calls.mix(mt_calls,x_par1_calls,x_nonpar_calls,x_par2_calls,x_calls,y_par1_calls,y_nonpar_calls,y_par2_calls).groupTuple().set{all_calls}

if (params.build == "b37"){
  process combine_gVCFs_b37 {
       tag { "${params.project_name}.${sample_id}.cCgVCFb37" }
       memory { 8.GB * task.attempt }
       label 'gatk'
       publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
  
       input:
       set val(sample_id), file(gvcf) from all_calls
  
       output:
  	   set val(sample_id), file("${sample_id}.g.vcf.gz") into combine_calls
  	   set val(sample_id), file("${sample_id}.g.vcf.gz.tbi") into combine_calls_indexes
  
       script:
       mem = task.memory.toGiga() - 4
       if (gvcf.size() == 29) // working with a male sample
       """
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       SortVcf \
       -I ${sample_id}.X_PAR1.g.vcf.gz \
       -I ${sample_id}.X_PAR2.g.vcf.gz \
       -I ${sample_id}.X_nonPAR.g.vcf.gz \
       -O ${sample_id}.X.g.vcf.gz
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       SortVcf \
       -I ${sample_id}.Y_PAR1.g.vcf.gz \
       -I ${sample_id}.Y_PAR2.g.vcf.gz \
       -I ${sample_id}.Y_nonPAR.g.vcf.gz \
       -O ${sample_id}.Y.g.vcf.gz
  
       echo "${gvcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo ${sample_id}.X.g.vcf.gz >> ${sample_id}.gvcf.list
       echo ${sample_id}.Y.g.vcf.gz >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.MT\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       GatherVcfs \
       -I ${sample_id}.gvcf.list \
       -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
  
       tabix -p vcf ${sample_id}.g.vcf.gz
       """
       else if (gvcf.size() == 24) // working with a female  sample
       """
       echo "${gvcf.join('\n')}" | grep "\\.1\\.g.vcf.gz" > ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.X\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.MT\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       GatherVcfs \
       -I ${sample_id}.gvcf.list \
       -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
  
       tabix -p vcf ${sample_id}.g.vcf.gz
       """
  }
}

if (params.build == "b38"){
  process combine_gVCFs_b38 {
       tag { "${params.project_name}.${sample_id}.cCgVCFb38" }
       memory { 8.GB * task.attempt }
       label 'gatk'
       publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
  
       input:
       set val(sample_id), file(gvcf) from all_calls
  
       output:
  	   set val(sample_id), file("${sample_id}.g.vcf.gz") into combine_calls
  	   set val(sample_id), file("${sample_id}.g.vcf.gz.tbi") into combine_calls_indexes
  
       script:
       mem = task.memory.toGiga() - 4
       if (gvcf.size() == 29) // working with a male sample
       """
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       SortVcf \
       -I ${sample_id}.X_PAR1.g.vcf.gz \
       -I ${sample_id}.X_PAR2.g.vcf.gz \
       -I ${sample_id}.X_nonPAR.g.vcf.gz \
       -O ${sample_id}.chrX.g.vcf.gz
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       SortVcf \
       -I ${sample_id}.Y_PAR1.g.vcf.gz \
       -I ${sample_id}.Y_PAR2.g.vcf.gz \
       -I ${sample_id}.Y_nonPAR.g.vcf.gz \
       -O ${sample_id}.chrY.g.vcf.gz
  
       echo "${gvcf.join('\n')}" | grep "\\.chr1\\.g.vcf.gz" > ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo ${sample_id}.chrX.g.vcf.gz >> ${sample_id}.gvcf.list
       echo ${sample_id}.chrY.g.vcf.gz >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chrM\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       GatherVcfs \
       -I ${sample_id}.gvcf.list \
       -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
  
       tabix -p vcf ${sample_id}.g.vcf.gz
       """
       else if (gvcf.size() == 24) // working with a female  sample
       """
       mv ${sample_id}.X.g.vcf.gz ${sample_id}.chrX.g.vcf.gz
       echo "${gvcf.join('\n')}" | grep "\\.chr1\\.g.vcf.gz" > ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr2\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr3\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr4\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr5\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr6\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr7\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr8\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr9\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr10\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr11\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr12\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr13\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr14\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr15\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr16\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr17\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr18\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr19\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr20\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr21\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chr22\\.g.vcf.gz" >> ${sample_id}.gvcf.list
       echo ${sample_id}.chrX.g.vcf.gz >> ${sample_id}.gvcf.list
       echo "${gvcf.join('\n')}" | grep "\\.chrMT\\.g\\.vcf\\.gz" >> ${sample_id}.gvcf.list
  
       gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g"  \
       GatherVcfs \
       -I ${sample_id}.gvcf.list \
       -O ${sample_id}.g.vcf.gz # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
  
       tabix -p vcf ${sample_id}.g.vcf.gz
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
