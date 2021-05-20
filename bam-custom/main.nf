#!/usr/bin/env nextflow

region = params.region

/// Get the bam info from the sheet
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
        }.set{samples}

process log_tool_version_samtools {
    tag { "${params.project_name}.ltV" }
    echo true
    publishDir "${params.out_dir}/custom-bam", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.samtools.version") into tool_version_samtool

    script:
    """
    samtools --version > tool.samtools.version
    """
}

process run_custom_bam {
    tag { "${params.project_name}.${sample_id}.rCB" }
    echo true
    publishDir "${params.out_dir}/custom-bam/${sample_id}", mode: 'move', overwrite: false
    label 'bwa_samtools'
    input:
    set val(sample_id), file(bam_file), file(bam_index) from samples

    output:
    set val(sample_id), file("${bam_file.getBaseName()}.*.bam") into custom_bam_file

    script:
    """
    a=$region
    b=\${a/:/_}
    samtools view \
    -b \
    ${bam_file} ${region} > ${bam_file.getBaseName()}.\$b.bam \
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
