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

process log_tool_version_multiqc {
    tag { "${params.project_name}.ltVm" }
    echo true
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
    echo true
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
    echo true
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
    echo true
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
    echo true
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
    echo true
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'samtools'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_1

    output:
       file("${bam_file}.flagstat") into flagstat_file

    script:
    """
    samtools flagstat \
    -@ 1 \
    ${bam_file} > ${bam_file}.flagstat
    """
}

process run_stats {
    tag { "${params.project_name}.${sample_id}.rS" }
    echo true
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
    echo true
    publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
    label 'mosdepth'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples_3
    file (ref) from ref_seq

    output:
       file("${sample_id}.mosdepth.*") into mosdepth_file

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
      echo true
      publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
      label 'verifybamid2'
      input:
      set val(sample_id), file(bam_file), file(bam_index) from samples_4
      file (ref) from ref_seq
  
      output:
         file("${sample_id}.result.*") into verifybamid2_file
  
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
      echo true
      publishDir "${params.out_dir}/cram-qc/${sample_id}", mode: 'copy', overwrite: false
      label 'verifybamid2'
      input:
      set val(sample_id), file(bam_file), file(bam_index) from samples_4
      file (ref) from ref_seq
  
      output:
         file("${sample_id}.result.*") into verifybamid2_file 
      
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
    echo true
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
