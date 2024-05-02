#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
known_indels_1            = file(params.known_indels_1, type: 'file')
known_indels_2            = file(params.known_indels_2, type: 'file')
dbsnp                     = file(params.dbsnp, type: 'file')
outdir                    = file(params.outdir, type: 'dir')
target_regions            = file(params.target_regions, type: 'file')
build                     = params.build
type                      = params.type
sample_coverage           = params.sample_coverage
project_name              = params.project_name
cohort_id                 = params.cohort_id
outdir.mkdir()

// 
if (build == "b37") {
    x = "X"
    y = "Y"
    // Coordinates are from here: https://www.ncbi.nlm.nih.gov/grc/human
    x_par1 = "60001-2699520"
    x_par1_target_default = "60001-60001"
    x_par2 = "154931044-155260560"
    x_par2_target_default = "154931044-154931044"
    y_par1 = "10001-2649520"
    y_par1_target_default = "10001-10001"
    y_par2 = "59034050-59363566"
    y_par2_target_default = "59034050-59034050"
    mt = "MT"
} else if (build == "b38"){
    x = "chrX"
    y = "chrY"
    // No autosomal region changes from build 37 to 38
    x_par1 = "10001-2781479"
    x_par1_target_default = "10001-10001"
    x_par2 = "155701383-156030895"
    x_par2_target_default = "155701383-155701383"
    y_par1 = "10001-2781479"
    y_par1_target_default = "10001-10001"
    y_par2 = "56887903-57217415"
    y_par2_target_default = "56887903-56887903"
    mt = "chrM"
} else {
    println "\n============================================================================================="
    println "Please specify a genome build (b37 or b38)!"
    println "=============================================================================================\n"
    exit 1
}

// CALL_CONF
if ( sample_coverage == "high" ) {
    call_conf = 30
} else if ( sample_coverage == "low" ) {
    call_conf = 10
} else {
    call_conf = 30
}

// RUN HAPLOTYPE CALLER ON AUTOSOMES
process run_haplotype_caller_auto {
    tag { "${sample_id}.${autosome}.rHCoA" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }

    input:
    tuple val(sample_id), val(gender), path(bam)
	each autosome
    
    output:
	tuple val(sample_id), path("${sample_id}.${autosome}.g.vcf.gz"), path("${sample_id}.${autosome}.g.vcf.gz.tbi"), emit: auto_calls

    script:
    mem = task.memory.toGiga() - 4
    if (type == "wes") {
    } else {
        intervals = "--intervals ${autosome}"
    }
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        HaplotypeCaller \
        --reference ${ref} \
        --input ${bam} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp} \
        ${intervals} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        -stand-call-conf ${call_conf} \
        --sample-ploidy 2 \
        -O ${sample_id}.${autosome}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON MALES
process run_haplotype_caller_males {
    tag { "${sample_id}:${nonautosome}.rHCoNAmales" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }

    input:
    tuple val(sample_id), val(gender), path(bam)
    each nonautosome
    
    output:
	tuple val(sample_id), path("${sample_id}.${base}.g.vcf.gz"), path("${sample_id}.${base}.g.vcf.gz.tbi"), emit: male_calls

    script:
    mem = task.memory.toGiga() - 4
    switch (nonautosome) {
        case['x_par1_male']:
            base = "X_PAR1"
            ploidy = "2"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals  ${x}:${x_par1}"
            }
            break
            // =====
        case['x_par2_male']:
            base = "X_PAR2"
            ploidy = "2"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals ${x}:${x_par2}"
            }
            break
            // =====
        case['x_nonpar_male']:
            base = "X_nonPAR"
            ploidy = "1"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals ${x} --exclude-intervals ${x}:${x_par1} --exclude-intervals ${x}:${x_par2}"
            }
            break
            // =====
        case['y_par1_male']:
            base = "Y_PAR1"
            ploidy = "2"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals ${y}:${y_par1}"
            }
            break
            // =====
        case['y_par2_male']:
            base = "Y_PAR2"
            ploidy = "2"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals ${y}:${y_par2}"
            }                    
            break
            // =====
        case['y_nonpar_male']:
            base = "Y_nonPAR"
            ploidy = "1"
            if (type == "wes") {
                target_regions = file(target_regions)
            } else {
                intervals = "--intervals ${y} --exclude-intervals ${y}:${y_par1} --exclude-intervals ${y}:${y_par2}"
            }
            break
            // =====
    }

    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        HaplotypeCaller \
        --reference ${ref} \
        --input ${bam} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp} \
        ${intervals} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        -stand-call-conf ${call_conf} \
        --sample-ploidy ${ploidy} \
        -O ${sample_id}.${base}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON FEMALES
process run_haplotype_caller_females {
    tag { "${sample_id}:${x}.rHCoNAfemales" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }

    input:
    tuple val(sample_id), val(gender), path(bam)
    
    output:
	tuple val(sample_id), path("${sample_id}.${x}.g.vcf.gz"), path("${sample_id}.${x}.g.vcf.gz.tbi"), emit: female_calls

    script:
    mem = task.memory.toGiga() - 4
    if (type == "wes") {
        target_regions = file(target_regions)
    } else {
        intervals = "--intervals ${x}"
    }
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        HaplotypeCaller \
        --reference ${ref} \
        --input ${bam} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp} \
        ${intervals} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        -stand-call-conf ${call_conf} \
        --sample-ploidy 2 \
        -O ${sample_id}.${x}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON MT
process run_haplotype_caller_mt {
    tag { "${sample_id}:${mt}.rHCoNAmt" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }

    input:
    tuple val(sample_id), val(gender), path(bam)
    
    output:
	tuple val(sample_id), path("${sample_id}.${mt}.g.vcf.gz"), path("${sample_id}.${mt}.g.vcf.gz.tbi"), emit: mt_calls

    script:
    mem = task.memory.toGiga() - 4
    if (type == "wes") {
        target_regions = file(target_regions)
    } else {
        intervals = "--intervals ${mt}"
    }
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        HaplotypeCaller \
        --reference ${ref} \
        --input ${bam} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp} \
        ${intervals} \
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        -stand-call-conf ${call_conf} \
        --sample-ploidy 2 \
        -O ${sample_id}.${mt}.g.vcf.gz
    """
}

process run_sort_male_gvcfs {
    tag { "${sample_id}.sMgVCF" }
    label 'gatk'
    memory { 8.GB * task.attempt }

    input:
    tuple val(sample_id), path(vcf), path(index) 
    
    output:
    tuple val(sample_id), path("${sample_id}.X.g.vcf.gz"), path("${sample_id}.Y.g.vcf.gz"), path("${sample_id}.X.g.vcf.gz.tbi"), path("${sample_id}.Y.g.vcf.gz.tbi"), emit: male_calls_combined
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g" SortVcf \
        --INPUT ${sample_id}.X_PAR1.g.vcf.gz \
        --INPUT ${sample_id}.X_PAR2.g.vcf.gz \
        --INPUT ${sample_id}.X_nonPAR.g.vcf.gz \
        --OUTPUT ${sample_id}.X.g.vcf.gz
  
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g" SortVcf \
        --INPUT ${sample_id}.Y_PAR1.g.vcf.gz \
        --INPUT ${sample_id}.Y_PAR2.g.vcf.gz \
        --INPUT ${sample_id}.Y_nonPAR.g.vcf.gz \
        --OUTPUT ${sample_id}.Y.g.vcf.gz
    """
}

process run_combine_sample_gvcfs {
    tag { "${sample_id}.cCgVCF" }
    label 'gatk'
    memory { 8.GB * task.attempt }
    publishDir "${outdir}/${params.workflow}/${project_name}/${sample_id}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), path(gvcfs), path(indexes)
  
    output:
	tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: combined_calls
    tuple val(sample_id), val("${outdir}/${params.workflow}/${project_name}/${sample_id}/${sample_id}.g.vcf.gz"), emit: gvcf
  
    script:
    mem = task.memory.toGiga() - 4
    if (build == 'b37') {
        chr = ''
        mt = 'MT'
    } else if(build == 'b38') {
        chr = 'chr'
        mt = 'chrM'
    }

    """
    find . -iname "${sample_id}.${chr}1.g.vcf.gz" > ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}2.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}3.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}4.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}5.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}6.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}7.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}8.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}9.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}10.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}11.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}12.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}13.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}14.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}15.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}16.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}17.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}18.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}19.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}20.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}21.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}22.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}X.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${chr}Y.g.vcf.gz" >> ${sample_id}.gvcf.list
    find . -iname "${sample_id}.${mt}.g.vcf.gz" >> ${sample_id}.gvcf.list

    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem}g" GatherVcfs \
        --INPUT ${sample_id}.gvcf.list \
        --OUTPUT ${sample_id}.g.vcf.gz
 
    # GatherVCF does not index the VCF. The VCF will be indexed in the next tabix operation.
    tabix -p vcf ${sample_id}.g.vcf.gz
    """
}

process run_create_gvcf_md5sum {
    tag { "${sample_id}.cGMD5" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/${params.workflow}/${project_name}/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(sample_id), path(gvcf), path(index)

    output:
    tuple val(sample_id), path("${gvcf}.md5"), path("${index}.md5"), emit: gvcf_md5sum

    """
    md5sum ${gvcf} > ${gvcf}.md5
    md5sum ${index} > ${index}.md5
    """
}

process run_create_gvcf_samplesheet {
    tag { "write_samplesheet" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/${params.workflow}/${project_name}/", mode: 'copy', overwrite: false

    input:
    path(samplesheet)

    output:
    path("samplesheet_${project_name}_gvcfs.tsv"), emit: bams_samplesheet

    """
    echo -e "SampleID\tGender\tFastqR1\tFastqR2\tFlowcell\tLane\tBAM\tgVCF" > samplesheet_${project_name}_gvcfs.tsv
    cat ${samplesheet} >> samplesheet_${project_name}_gvcfs.tsv
    """
}
