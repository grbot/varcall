#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//========== PARAMETERS & INPUT ==========//
// ALIGNMENT
params.build              = "b37"
params.dbsnp              = "/home/phelelani/data/gatk-b37-bundle/dbsnp_138.b37.vcf"
params.known_indels_1     = "/home/phelelani/data/gatk-b37-bundle/1000G_phase1.indels.b37.vcf"
params.known_indels_2     = "/home/phelelani/data/gatk-b37-bundle/Mills_and_1000G_gold_standard.indels.b37.vcf"
params.hapmap             = "/home/phelelani/data/gatk-b37-bundle/hapmap_3.3.b37.vcf"
params.omni               = "/home/phelelani/data/gatk-b37-bundle/1000G_omni2.5.b37.vcf"
params.phase1_snps        = "/home/phelelani/data/gatk-b37-bundle/1000G_phase1.snps.high_confidence.b37.vcf"
params.golden_indels      = "/home/phelelani/data/gatk-b37-bundle/Mills_and_1000G_gold_standard.indels.b37.vcf"

params.memory_per_thread  = "500M"
params.outdir             = "/home/phelelani/nf-workflows/varcall/WF-GENERATE_GVCFS"
params.project_name       = "NA12878"

params.ref                = "/home/phelelani/data/gatk-b37-bundle/human_g1k_v37_decoy.fasta" 
params.sample_sheet       = "/home/phelelani/nf-workflows/varcall/data/NA12878.samplesheet.tsv"
params.workflow           = null
params.bwa_threads        = "32"
params.type               = "wgs"
params.sample_coverage    = "high"
params.target_regions     = ""

autosomes = (1..22).toList()
println autosomes + ["X","Y","MT"]

// SOME SMALL FUNCTIONS

// CHANNELS
// Get sample info from sample sheet - Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
samples = Channel.fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> [ "${row.SampleID}", "${row.Gender}", "${row.FastqR1}", "${row.FastqR2}", "${row.Flowcell}", "${row.Lane}", "${row.BAM}", "${row.gVCF}" ] }

// GET INPUT FOR ALIGN WORKFLOW: SAMPLE_ID, FORWARD_READ, REVERSE_READ, FLOWCELL, LANE
samples.map { [ it[0], it[2], it[3], it[4], it[5]] }
    .set { samples_aling }

// GET INPUT FOR GENERATE-GVCFS WORKFLOW: SAMPLE_ID, GENDER, BAM
samples.map { [ it[0], it[1], it[6] ] }
    .set { samples_generate_gvcfs }

// GET INPUT (GENDER BASED) FOR GENERATE-GVCFS WORKFLOW: SAMPLE_ID, GENDER, BAM
samples_generate_gvcfs.filter { it[1] == 'M' }.set { samples_generate_gvcfs_males }
samples_generate_gvcfs.filter { it[1] == 'F' }.set { samples_generate_gvcfs_females }
samples_generate_gvcfs.filter { !(it[1] =~ 'F|M') }.set { samples_generate_gvcfs_nosex } 

// WORKFLOWS PARAMETER
workflow = params.workflow

// GENERATE GVCFS WORKFLOW
if (params.build == "b37") {
    autosomes = (1..22).toList()
} else if (params.build == "b38"){
    autosomes = (1..22).toList().collect { 'chr' + "${it}" }
}
nonautosomes = ["x_par1_male", "x_par2_male", "x_nonpar_male", "y_par1_male", "y_par2_male", "y_nonpar_male" ]

// INCLUDE MODULES
include { log_tool_version_samtools_bwa; log_tool_version_gatk;
         run_bwa; run_mark_duplicates; run_create_recalibration_table;
         run_recalibrate_bam; bam_to_cram; run_cram_flagstat; create_cram_md5sum } from './modules/modules-align.nf'
include { run_haplotype_caller_on_autosomes as run_haplotype_caller_on_autosomes_nosex;
         run_haplotype_caller_on_autosomes as run_haplotype_caller_on_autosomes_males;
         run_haplotype_caller_on_autosomes as run_haplotype_caller_on_autosomes_females;
         run_haplotype_caller_on_xy;
         run_haplotype_caller_on_xx;
         run_haplotype_caller_on_mt as run_haplotype_caller_on_mt_nosex;
         run_haplotype_caller_on_mt as run_haplotype_caller_on_mt_males;
         run_haplotype_caller_on_mt as run_haplotype_caller_on_mt_females;
         run_sort_xy_gVCFs;
         run_combine_gVCFs as run_combine_gVCFs_nosex;
         run_combine_gVCFs as run_combine_gVCFs_males;
         run_combine_gVCFs as run_combine_gVCFs_females;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_nosex;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_males;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_females
} from './modules/modules-generate-gvcf.nf'
// include {} from './modules/modules-qc.nf'
// include {} from './modules/modules-variant-calling.nf'

// WORKFLOWS
workflow ALIGN {
    take:
    samples

    main:
    // log_tool_version_samtools_bwa()
    // log_tool_version_gatk()
    run_bwa(samples)
    run_mark_duplicates(run_bwa.out.raw_bam)
    run_create_recalibration_table(run_mark_duplicates.out.md_bam)
    run_recalibrate_bam(run_create_recalibration_table.out.recal_table)
    bam_to_cram(run_recalibrate_bam.out.recal_bam)
    run_cram_flagstat(bam_to_cram.out.cram_file)
    create_cram_md5sum(run_cram_flagstat.out.cram_stats)
}

workflow GENERATE_GVCFS {
    take:
    samples_nosex
    samples_males
    samples_females
    autosomes
    nonautosomes
    
    main:
    // log_tool_version_gatk()    

    // GENERATE_GVCFS: NO SEX
    run_haplotype_caller_on_autosomes_nosex(samples_nosex, autosomes)
    run_haplotype_caller_on_mt_nosex(samples_nosex)
    run_haplotype_caller_on_autosomes_nosex.out.autosome_calls
        .groupTuple()
        .join(run_haplotype_caller_on_mt_nosex.out.nonautosome_calls_mt)
        .flatten()
        .collect()
        .map { it -> [ it[0], it.findAll { it =~ 'g.vcf.gz$' }, it.findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .view()
        .set { nosex }
    run_combine_gVCFs_nosex(nosex).view()
    run_create_gvcf_md5sum_nosex(run_combine_gVCFs_nosex.out.combined_calls)

    // GENERATE_GVCFS: MALES
    run_haplotype_caller_on_autosomes_males(samples_males, autosomes)
    run_haplotype_caller_on_mt_males(samples_males)
    run_haplotype_caller_on_xy(samples_males, nonautosomes)
    run_sort_xy_gVCFs(run_haplotype_caller_on_xy.out.nonautosome_calls_males.groupTuple())
    run_haplotype_caller_on_autosomes_males.out.autosome_calls
        .groupTuple()
        .join(run_haplotype_caller_on_mt_males.out.nonautosome_calls_mt)
        .join(run_sort_xy_gVCFs.out.nonautosome_calls_males_combined)
        .flatten()
        .collect()
        .map { it -> [ it[0], it.findAll { it =~ 'g.vcf.gz$' }, it.findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { males }
    run_combine_gVCFs_males(males)
    run_create_gvcf_md5sum_males(run_combine_gVCFs_males.out.combined_calls)

    // GENERATE_GVCFS: FEMALES
    run_haplotype_caller_on_autosomes_females(samples_females, autosomes)
    run_haplotype_caller_on_mt_females(samples_females)
    run_haplotype_caller_on_xx(samples_females)
    run_haplotype_caller_on_autosomes_females.out.autosome_calls
        .groupTuple()
        .join(run_haplotype_caller_on_mt_females.out.nonautosome_calls_mt)
        .join(run_haplotype_caller_on_xx.out.nonautosome_calls_females)
        .flatten()
        .collect()
        .map { it -> [ it[0], it.findAll { it =~ 'g.vcf.gz$' }, it.findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { females }
    run_combine_gVCFs_females(females)
    run_create_gvcf_md5sum_females(run_combine_gVCFs_females.out.combined_calls)
}

// PICK & CHOOSE WORKFLOW
workflow {
    switch (workflow) {
            // ===== MAIN WORKFLOWS
        case['align']:
            ALIGN(samples_align)
            break
            // =====
        case['generate-gvcf']:
            GENERATE_GVCFS(samples_generate_gvcfs_nosex, samples_generate_gvcfs_males, samples_generate_gvcfs_females, autosomes, nonautosomes)
            break
            // =====
        case['combine-gvcfs']:
            break
            // =====
        case['genomics-db-import']:
            break
            // =====
        case['genome-calling']:
            break
            // =====            
        case['filter-vcf ']:
            break
            // ===== SEPARATE WORKFLOWS
        case['bam-to-cram']:
            break
            // =====            
        case['cram-to-bam']:
            break
            // =====            
        case['cram-to-fastq']:
            break
            // =====            
        case['index-bams']:
            break
            // =====            
        case['index-vcf']:
            break
            // =====            
        case['bam-flagstat']:
            break
            // =====            
        case['mt-calling']:
            break
            // =====            
        case['combine-lanes']:
            break
            // =====            
        case['validate-gvcf']:
            break
            // =====            
        case['genotype-refinement']:
            break
            // =====            
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
            // =====
    }
}
