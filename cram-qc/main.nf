#!/usr/bin/env nextflow
// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, Gender and BAM file location (pointing to CRAM)
/// Retrieving the indexing is very specific to the way we name our files
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def bam_file = file(row['BAM/CRAM'])
	    def bam_index = ''
            if(bam_file.getName().contains('.bam')){
              bam_index = file(bam_file + ".bai")
	      if(!bam_index.exists()){
                bam_index = file((bam_file.getParent()).toString() + "/" + (bam_file.getBaseName()).toString() + ".bai")
              }
	    }
            if(bam_file.getName().contains('.cram')){
              bam_index = file(bam_file + ".crai")
	    }
            return [ sample_id, bam_file, bam_index ]
        }.into{samples_1; samples_2; samples_3; samples_4; samples_5}


ref_seq = Channel.fromPath(params.ref_seq).toList()
ref_dict = Channel.fromPath(params.ref_dict).toList()
high_cov_cutoff = (params.high_cov_cutoff)
medium_cov_cutoff = (params.medium_cov_cutoff)
dup_cutoff = (params.dup_cutoff)
map_cutoff = (params.map_cutoff)
freemix_cutoff = (params.freemix_cutoff)

process log_tool_version_multiqc {
    tag { "${params.project_name}.ltVm" }
    publishDir "${params.out_dir}/cram-qc", mode: 'copy', overwrite: false
    label 'multiqc'

    output:
    file("tool.multiqc.version") into tool_version_multiqc

    script:
    """
    multiqc --version > tool.multiqc.version
    """
}

process log_tool_version_samtools {
    tag { "${params.project_name}.ltVs" }
    publishDir "${params.out_dir}/cram-qc", mode: 'copy', overwrite: false
    label 'samtools'

    output:
    file("tool.samtools.version") into tool_version_samtool

    script:
    """
    samtools 2>&1 | grep Version > tool.samtools.version
    """
}

process log_tool_version_gatk {
    tag { "${params.project_name}.ltVg" }
    publishDir "${params.out_dir}/cram-qc", mode: 'copy', overwrite: false
    label 'gatk'

    output:
    file("tool.gatk.version") into tool_version_gatk

    script:
    """
    gatk --version > tool.gatk.version
    """
}

process log_tool_version_mosdepth {
    tag { "${params.project_name}.ltVm" }
    publishDir "${params.out_dir}/cram-qc", mode: 'copy', overwrite: false
    label 'mosdepth'

    output:
    file("tool.mosdepth.version") into tool_version_mosdepth

    script:
    """
    mosdepth --version > tool.mosdepth.version
    """
}

process log_tool_version_verifybamid2 {
    tag { "${params.project_name}.ltVv" }
    publishDir "${params.out_dir}/cram-qc", mode: 'copy', overwrite: false
    label 'verifybamid2'

    output:
    file("tool.verifybamid2.version") into tool_version_verifybamid2

    script:
    """
    verifybamid2 2>&1 | grep Version | sed 's/ //' > tool.verifybamid2.version
    """
}

process run_flagstat {
    tag { "${params.project_name}.${sample_id}.rF" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'samtools'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_1

    output:
       file("${bam_file}.flagstat") into flagstat_file
       set val(sample_id), file("${bam_file}.flagstat") into flagstat

    script:
    """
    samtools flagstat \
    -@ 1 \
    ${bam_file} > ${bam_file}.flagstat
    """
}

process run_stats {
    tag { "${params.project_name}.${sample_id}.rS" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'samtools'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_2

    output:
       file("${bam_file}.stats") into stats_file

    script:
    """
    samtools stats \
    -@ 1 \
    ${bam_file} > ${bam_file}.stats
    """
}

process run_mosdepth {
    tag { "${params.project_name}.${sample_id}.riM" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'mosdepth'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_3
    file (ref) from ref_seq

    output:
       file("${sample_id}.mosdepth.*") into mosdepth_file
       set val(sample_id), file("${sample_id}.mosdepth.summary.txt") into mosdepth

    script:
    """
    mosdepth \
    -f $ref \
    $sample_id \
    ${bam_file}
    """
}

if (params.build == "b37") {
  process run_verifybamid2_b37 {
      tag { "${params.project_name}.${sample_id}.riV" }
      publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
      label 'verifybamid2'
      input:
      set val(sample_id), file(bam_file), file(bam_index) from samples_4
      file (ref) from ref_seq
  
      output:
         file("${sample_id}.result.*") into verifybamid2_file
         set val(sample_id), file("${sample_id}.result.selfSM") into verifybamid2
  
      script:
      """
      verifybamid2 \
      --SVDPrefix  /usr/local/share/verifybamid2-2.0.1-7/resource/1000g.phase3.100k.b37.vcf.gz.dat \
      --Reference ${ref} \
      --BamFile ${bam_file}
      mv result.Ancestry ${sample_id}.result.Ancestry
      mv result.selfSM ${sample_id}.result.selfSM
      """
  }
}

if (params.build == "b38") {
  process run_verifybamid2_b38 {
      tag { "${params.project_name}.${sample_id}.riV" }
      publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
      label 'verifybamid2'
      input:
      set val(sample_id), file(bam_file), file(bam_index) from samples_4
      file (ref) from ref_seq
        
      output:
         file("${sample_id}.result.*") into verifybamid2_file
         set val(sample_id), file("${sample_id}.result.selfSM") into verifybamid2
      
      script:
      """
      verifybamid2 \
      --SVDPrefix /usr/local/share/verifybamid2-2.0.1-7/resource/1000g.phase3.100k.b38.vcf.gz.dat \
      --Reference ${ref} \
      --BamFile ${bam_file}
      mv result.Ancestry ${sample_id}.result.Ancestry
      mv result.selfSM ${sample_id}.result.selfSM
      """
  }
}

// Not doing GATK CollectAlignmentSummaryMetrics now. Seems like it requires SAM/BAM input.
/*
process run_gatk {
    tag { "${params.project_name}.${sample_id}.riV" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'gatk'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_5
    file (ref) from ref_seq
    file (dict) from ref_dict

    output:
       file("*.gatk.*") into gatk_file

    script:
    """
    gatk CollectAlignmentSummaryMetrics \
    -I ${bam_file} \
    -R ${ref} \
    -H ${sample_id}.gatk.pdf \
    -O ${sample_id}.gatk.metrics 
    """
}
*/

combined = flagstat_file.mix(stats_file, mosdepth_file, verifybamid2_file)
//combined = flagstat_file.mix(stats_file, mosdepth_file, verifybamid2_file, gatk_file)

process runMultiQC{
    tag { "${params.project_name}.rMQC" }
    publishDir "${params.out_dir}/cram-qc/", mode: 'move', overwrite: false
    label 'multiqc'

    input:
        file('*') from combined.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

process check_dup_cov {
    tag { "${params.project_name}.${sample_id}.cDC" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    input:
    set val(sample_id), file(flagstat_file) from flagstat

    output:
       file("${sample_id}.dup_cov") into dup_cov

    script:
    """
    total=`cat $flagstat_file |  grep "total" | cut -d ' ' -f1`
    mapped=`cat $flagstat_file |  grep "mapped (" | cut -d ' ' -f1`
    dup=`cat $flagstat_file |  grep "duplicates" | cut -d ' ' -f1`
    perc_mapped=\$((mapped*100/total))
    perc_dup=\$((dup*100/total))

     if [ \$perc_mapped -gt ${map_cutoff} ]
     then
       echo -e "${sample_id}\tMapping cutoff\tPASS" > ${sample_id}.dup_cov 
     else
       echo -e "${sample_id}\tMapping cutoff\tFAIL" > ${sample_id}.dup_cov
     fi
	
     if [ \$perc_dup -lt ${dup_cutoff} ]
     then
       echo -e "${sample_id}\tDuplication cutoff\tPASS" >> ${sample_id}.dup_cov
     else
       echo -e "${sample_id}\tDuplication cutoff\tFAIL" >> ${sample_id}.dup_cov
     fi

    """
}

process check_freemix {
    tag { "${params.project_name}.${sample_id}.cFM" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    input:
    set val(sample_id), file(verifybamid2_file) from verifybamid2

    output:
       file("${sample_id}.freemix") into freemix

    script:
    """
    freemix=`cat $verifybamid2_file | tail -1 | cut -f 7`
    freemix_mod=`echo \$freemix | sed 's/e/*10^/'`
  
 
    if [ \$(echo "\$freemix_mod < $freemix_cutoff"|bc) -eq 1 ]
     then
       echo -e "${sample_id}\tFreemix cutof\tPASS" > ${sample_id}.freemix 
     else
       echo -e "${sample_id}\tFreemix cutoff\tFAIL" > ${sample_id}.freemix
     fi

    """
}

process check_depth {
    tag { "${params.project_name}.${sample_id}.cD" }
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    input:
      set val(sample_id), file(mosdepth_file) from mosdepth
    
    output:
       file("${sample_id}.depth")into depth

    script:
    """
    depth=`cat $mosdepth_file | tail -1 | cut -f 4`
    depth=`printf "%.0f\n" \$depth`
 
    if [ \$depth -ge $high_cov_cutoff ]
    then
      echo -e "${sample_id}\tDepth\tHigh coverage" > ${sample_id}.depth 
    elif [ \$depth -ge $medium_cov_cutoff ] && [ \$depth -lt $high_cov_cutoff ]
    then
      echo -e "${sample_id}\tDepth\tMedium coverage" > ${sample_id}.depth
    else
      echo -e "${sample_id}\tDepth\tLow coverage" > ${sample_id}.depth
    fi

    """
}

checks = dup_cov.mix(freemix, depth)

process combine_checks {
    tag { "${params.project_name}.cC" }
    publishDir "${params.out_dir}/cram-qc/", mode: 'copy', overwrite: false
    input:
      file('*') from checks.collect()

    output:
       file ("checks.tsv")

    script:
    """
      cat * | sort -k1 > checks.tsv
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
