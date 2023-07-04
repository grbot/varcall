process select_samples {
    tag { "${params.project_name}.${params.protocol}.select_samples.${name}" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bcftools'

    input:
    tuple val(name), file(vcf)
    file(sample_list)
  
    output:
    tuple val("${name}.select_samples"), file("${name}.select_samples.vcf.gz*") 

    script:
    """
    bcftools view --threads ${task.cpus} -c 1 -O z -S ${sample_list} --force-samples -o ${name}.select_samples.vcf.gz ${vcf[0]}
    bcftools index --threads ${task.cpus} -t ${name}.select_samples.vcf.gz
    """ 
}

process get_stats {
    tag { "${params.project_name}.${params.protocol}.get_stats" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bcftools'

    input:
    tuple val(name_1), file(vcf_1)
    tuple val(name_2), file(vcf_2)
    file(sample_list)
  
    output:
    file("${name_1}_vs_${name_2}.bcftools.stats") 

    script:
    """
    bcftools stats -v -S ${sample_list} ${vcf_1[0]} ${vcf_2[0]} > ${name_1}_vs_${name_2}.bcftools.stats
    """ 
}

process get_unique_sites_and_overlaps {
    tag { "${params.project_name}.${params.protocol}.get_unique_sites_and_overlaps" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bcftools'

    input:
    tuple val(name_1), file(vcf_1)
    tuple val(name_2), file(vcf_2)
    file(sample_list)
  
    output:
    file("${name_1}_vs_${name_2}.bcftools.isec/*.vcf") 

    script:
    """
    bcftools isec -p ${name_1}_vs_${name_2}.bcftools.isec ${vcf_1[0]} ${vcf_2[0]}
    """ 
}

process get_unique_sites_and_overlap_index {
    tag { "${params.project_name}.${params.protocol}.get_unique_sites_and_overlaps_index" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bcftools'

    input:
    file isec
  
    output:
    file("${isec.getSimpleName()}.sorted.vcf.gz*") 

    script:
    """
    tmp="/scratch3/users/gerrit/tmp/\$RANDOM"
    mkdir -p \$tmp
    bcftools view -O z -o ${isec.getSimpleName()}.vcf.gz ${isec}
    bcftools index -t ${isec.getSimpleName()}.vcf.gz
    bcftools sort -O z -o ${isec.getSimpleName()}.sorted.vcf.gz ${isec}.gz -T \$tmp
    bcftools index -t ${isec.getSimpleName()}.sorted.vcf.gz
    rm -rf \$tmp
    """ 
}

process get_novel_stats {
    tag { "${params.project_name}.${params.protocol}.get_novel_stats" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bcftools'

    input:
    tuple val(name), file(dbsnp)
    each file(isec)
  
    output:
    file("${isec[0].getName()}_vs_${name}.bcftools.stats")
    file("${isec[0].getName()}_vs_${name}.bcftools.isec/*.vcf") // will not be used just including for work to be done later

    script:
    """
    bcftools stats -v ${isec[0]} ${dbsnp[0]} > ${isec[0].getName()}_vs_${name}.bcftools.stats
    bcftools isec -p ${isec[0].getName()}_vs_${name}.bcftools.isec ${isec[0]} ${dbsnp[0]}
    """ 
}

process combine_stats_novel {
    tag { "${params.project_name}.${params.protocol}.combine_stats_novel" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false

    input:
    path bcftools_stats
    path vcf // will not be used just including for work to be done later

    output:
    file("final.stats_novel.txt")

    script:
    """
    for i in `ls -1 *.bcftools.stats`; do
      echo -en "Records in \$i that are not novel\t" >> final.stats_novel.txt
      cat \$i | grep ^SN | awk '{if(\$2=='2')print \$0}'  | grep "number of records" | cut -f 4 >> final.stats_novel.txt
    done  
    """
}

process parse_stats {
    tag { "${params.project_name}.${params.protocol}.parse_stats" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false

    input:
    path stats_main 
    path stats_novel_combined

    output:
    file("final.stats.txt")

    script:
    """
    echo -en "NRD SNPs (%)\t" >> final.stats.txt
    grep "NRD" ${stats_main} | grep -A2 "SNPs" | grep "^NRDs" | cut -f3 >> final.stats.txt
    echo -en "NRD INDELs (%)\t" >> final.stats.txt
    grep "NRD" ${stats_main} |  grep -A2 "indels" | grep "^NRDi" | cut -f3 >> final.stats.txt
    echo -en "Records unique to 1st VCF\t" >> final.stats.txt
    cat ${stats_main} | grep ^SN | awk '{if(\$2=='0')print \$0}'  | grep "number of records" | cut -f 4 >> final.stats.txt
    echo -en "Records unique to 2nd VCF\t" >> final.stats.txt
    cat ${stats_main} | grep ^SN | awk '{if(\$2=='1')print \$0}'  | grep "number of records" | cut -f 4 >> final.stats.txt
    echo -en "Records that overlap\t" >> final.stats.txt
    cat ${stats_main} | grep ^SN | awk '{if(\$2=='2')print \$0}'  | grep "number of records" | cut -f 4 >> final.stats.txt
    cat ${stats_novel_combined} >> final.stats.txt
    """
}

