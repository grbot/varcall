#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//========== PARAMETERS & INPUT ==========//
// ALIGNMENT
params.sample_sheet       = "/home/phelelani/nf-workflows/varcall/data/samplesheets/samplesheet_scd.tsv"
params.outdir             = "/external/diskC/phelelani/varcall_out"
params.project_name       = "scd"
params.cohort_id          = "scd"
params.sample_coverage    = "low"

params.workflow           = null
params.type               = "wgs"
params.target_regions     = "/"

//params.db_path            = "/external/diskC/phelelani/data/b38db"
params.db_path            = "/local/phelelani/b38db"
params.db_update          = "no"

params.gvcf               = "/home/phelelani/nf-workflows/varcall/WF-GENERATE_GVCFS/combine-gvcfs/TEST_COHORT.g.vcf.gz"

// // BUILD 37
// params.build              = "b37"
// params.ref                = "/external/diskC/phelelani/data/gatk/bundle/b37/human_g1k_v37_decoy.fasta" 
// params.dbsnp              = "/external/diskC/phelelani/data/gatk/bundle/b37/dbsnp_138.b37.vcf"
// params.known_indels_1     = "/external/diskC/phelelani/data/gatk/bundle/b37/1000G_phase1.indels.b37.vcf"
// params.known_indels_2     = "/external/diskC/phelelani/data/gatk/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
// params.hapmap             = "/external/diskC/phelelani/data/gatk/bundle/b37/hapmap_3.3.b37.vcf"
// params.omni               = "/external/diskC/phelelani/data/gatk/bundle/b37/1000G_omni2.5.b37.vcf"
// params.phase1_snps        = "/external/diskC/phelelani/data/gatk/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf"
// params.golden_indels      = "/external/diskC/phelelani/data/gatk/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"

// BUILD 38
params.build              = "b38"
params.ref                = "/external/diskC/phelelani/data/gatk/bundle/b38/Homo_sapiens_assembly38.fasta" 
params.dbsnp              = "/external/diskC/phelelani/data/gatk/bundle/b38/dbsnp_146.hg38.vcf.gz"
params.known_indels_1     = "/external/diskC/phelelani/data/gatk/bundle/b38/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.known_indels_2     = "/external/diskC/phelelani/data/gatk/bundle/b38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.hapmap             = "/external/diskC/phelelani/data/gatk/bundle/b38/hapmap_3.3.hg38.vcf.gz"
params.omni               = "/external/diskC/phelelani/data/gatk/bundle/b38/1000G_omni2.5.hg38.vcf.gz"
params.phase1_snps        = "/external/diskC/phelelani/data/gatk/bundle/b38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.golden_indels      = "/external/diskC/phelelani/data/gatk/bundle/b38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

params.memory_per_thread  = "500M"
params.bwa_threads        = "32"

params.ts_filter_level_snps  = 99.5    
params.ts_filter_level_indels = 99.0    
params.max_gaussians_snps = 4
params.max_gaussians_indels = 8

//----- COMMON  WORKFLOWS PARAMETER
workflow                  = params.workflow
db                        = file(params.db_path, type: 'dir')
db_update                 = params.db_update

// GENERATE GVCFS WORKFLOW
if (params.build == "b37") {
    if ( params.workflow == 'generate-gvcfs' ) {
        chroms_auto = (1..22).toList()
        chroms_par = ["x_par1_male", "x_par2_male", "x_nonpar_male", "y_par1_male", "y_par2_male", "y_nonpar_male" ]
    } else if ( params.workflow == 'genome-calling' || params.workflow == 'combine-gvcfs' || params.workflow == 'genomics-db-import' ) {
        chroms_all = (1..22).toList() +  ["X","Y","MT"]
    }
} else if (params.build == "b38") {
    if ( params.workflow == 'generate-gvcfs' ) {
        chroms_auto = (1..22).toList().collect { 'chr' + "${it}" }
        chroms_par = ["x_par1_male", "x_par2_male", "x_nonpar_male", "y_par1_male", "y_par2_male", "y_nonpar_male" ]
    } else if ( params.workflow == 'genome-calling' || params.workflow == 'combine-gvcfs' || params.workflow == 'genomics-db-import') {
        chroms_all = (1..22).toList().collect { 'chr' + "${it}" } + ["chrX", "chrY","chrM"]
    }
} else {
}

//----- CHANNELS
// GET SAMPLE INFO FROM SAMPLE SHEET - MINIMUM INFORMATION THAT ARE NEEDED IN THE SAMPLE SHEET ARE SAMPLEID, FORWARD AND REVERSE READ FILE LOCATION
samples = Channel.fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> [ "${row.SampleID}",
                   "${row.Gender}",
                   "${row.FastqR1}",
                   "${row.FastqR2}",
                   "${row.Flowcell}",
                   "${row.Lane}",
                   "${row.BAM}",
                   "${row.gVCF}" ] }

// GET INPUT FOR ALIGN WORKFLOW: SAMPLE_ID, FORWARD_READ, REVERSE_READ, FLOWCELL, LANE
samples.map { [ it[0], it[2], it[3], it[4], it[5]] }
    .set { samples_align }

// GET INPUT FOR GENERATE-GVCFS WORKFLOW: SAMPLE_ID, GENDER, BAM
samples.map { [ it[0], it[1], it[6]] }
    .set { samples_generate_gvcfs }

// GET INPUT (GENDER BASED) FOR GENERATE-GVCFS WORKFLOW: SAMPLE_ID, GENDER, BAM
samples_generate_gvcfs.filter { it[1] == 'M' }.set { samples_generate_gvcfs_males }
samples_generate_gvcfs.filter { it[1] == 'F' }.set { samples_generate_gvcfs_females }
samples_generate_gvcfs.filter { !(it[1] =~ 'F|M') }.set { samples_generate_gvcfs_nosex } 

// GET INPUT FOR COMBINE-GVCFS WORKFLOW: SAMPLE_ID, GVCF
samples.map { [ it[0], it[7] ] }
    .set { samples_gvcfs_list }

// GVCF FILE ---CHECK<<<<<<<
gvcf = Channel.fromPath([params.gvcf, "${params.gvcf}.tbi"]).toList()

//----- INCLUDE MODULES
// GENERAL
include { print_sample_info; log_tool_version_samtools_bwa; log_tool_version_gatk } from './modules/modules-general.nf'
// ALIGNMENT
include { run_bwa; run_mark_duplicates; run_create_recalibration_table;
         run_recalibrate_bam; bam_to_cram; run_cram_flagstat; create_cram_md5sum } from './modules/modules-align.nf'
// GENERATE GENOME GVCFS
include { run_haplotype_caller_auto as run_haplotype_caller_auto_nosex;
         run_haplotype_caller_auto as run_haplotype_caller_auto_males;
         run_haplotype_caller_auto as run_haplotype_caller_auto_females;
         run_haplotype_caller_males;
         run_haplotype_caller_females;
         run_haplotype_caller_mt as run_haplotype_caller_mt_nosex;
         run_haplotype_caller_mt as run_haplotype_caller_mt_males;
         run_haplotype_caller_mt as run_haplotype_caller_mt_females;
         run_sort_male_gvcfs;
         run_combine_sample_gvcfs as run_combine_sample_gvcfs_nosex;
         run_combine_sample_gvcfs as run_combine_sample_gvcfs_males;
         run_combine_sample_gvcfs as run_combine_sample_gvcfs_females;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_nosex;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_males;
         run_create_gvcf_md5sum as run_create_gvcf_md5sum_females } from './modules/modules-generate-gvcf.nf'
// COMBINE GVCFS
include { run_combine_gvcfs; run_concat_gvcfs} from './modules/modules-combine-gvcfs.nf'
// GENOMICS DB IMPORT
include { run_genomics_db_import_new; run_backup_genomic_db; run_genomics_db_import_update } from './modules/modules-genomics-db-import.nf'
// GENOME CALLING
include { run_genotype_gvcf_on_genome_gvcf; run_concat_vcf; run_vqsr_on_snps; run_apply_vqsr_on_snps;
         run_vqsr_on_indels; run_apply_vqsr_on_indels } from './modules/modules-genome-calling.nf'

// ALIGN WORKFLOW - NOT TESTED THOROUGHLY (NEED ORIGINAL FASTQ FILES)
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

// GENERATE GVCFS WORKFLOW - TESTED THOROUGHLY AND WORKS FINE
workflow GENERATE_GVCFS {
    take:
    samples_nosex
    samples_males
    samples_females
    chroms_auto
    chroms_par
    
    main:
    // GENERATE_GVCFS: NO SEX
    run_haplotype_caller_auto_nosex(samples_nosex, chroms_auto)
    run_haplotype_caller_mt_nosex(samples_nosex)
    run_haplotype_caller_auto_nosex.out.auto_calls.groupTuple()
        .join(run_haplotype_caller_mt_nosex.out.mt_calls.groupTuple())
        .map { it -> [ it[0], it.flatten().findAll { it =~ 'g.vcf.gz$' }, it.flatten().findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { nosex }
    run_combine_sample_gvcfs_nosex(nosex)
    run_create_gvcf_md5sum_nosex(run_combine_sample_gvcfs_nosex.out.combined_calls)
    // GENERATE_GVCFS: MALES
    run_haplotype_caller_auto_males(samples_males, chroms_auto)
    run_haplotype_caller_mt_males(samples_males)
    run_haplotype_caller_males(samples_males, chroms_par)
    run_sort_male_gvcfs(run_haplotype_caller_males.out.male_calls.groupTuple())
    run_haplotype_caller_auto_males.out.auto_calls.groupTuple()
        .join(run_haplotype_caller_mt_males.out.mt_calls.groupTuple())
        .join(run_sort_male_gvcfs.out.male_calls_combined.groupTuple())
        .map { it -> [ it[0], it.flatten().findAll { it =~ 'g.vcf.gz$' }, it.flatten().findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { males }
    run_combine_sample_gvcfs_males(males)
    run_create_gvcf_md5sum_males(run_combine_sample_gvcfs_males.out.combined_calls)
    // GENERATE_GVCFS: FEMALES
    run_haplotype_caller_auto_females(samples_females, chroms_auto)
    run_haplotype_caller_mt_females(samples_females)
    run_haplotype_caller_females(samples_females)
    run_haplotype_caller_auto_females.out.auto_calls.groupTuple()
        .join(run_haplotype_caller_mt_females.out.mt_calls.groupTuple())
        .join(run_haplotype_caller_females.out.female_calls.groupTuple())
        .map { it -> [ it[0], it.flatten().findAll { it =~ 'g.vcf.gz$' }, it.flatten().findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { females }
    run_combine_sample_gvcfs_females(females)
    run_create_gvcf_md5sum_females(run_combine_sample_gvcfs_females.out.combined_calls)
}

// COMBINE GVCFS - SKIPPING TESTING 
workflow COMBINE_GVCFS {
    take:
    samples_gvcfs_list
    chrom_all
    
    main:
    // log_tool_version_gatk()
    samples_gvcfs_list
        .collectFile() { item -> [ 'gvcf.list', "${item.get(1)}" + '\n' ] }
        .set { gvcf_list }
    run_combine_gvcfs( gvcf_list, chrom_all)
    run_combine_gvcfs.out.cohort_chr_calls
        .flatten()
        .collect()
        .map { it -> [ it.findAll { it =~ 'g.vcf.gz$' }, it.findAll { it =~ 'g.vcf.gz.tbi$' } ] }
        .set { cohort_calls }
    run_concat_gvcfs(cohort_calls)
}

// GENOMICS DB IMPORT
workflow GENOMICS_DB_IMPORT {
    take:
    samples_gvcfs_list
    chrom_all

    main:
    switch (db_update) {
        case['yes']:
            if ( db.exists() ) {
                println "DB:${db} exists. Making a backup copy before update."
                run_backup_genomic_db()
                samples_gvcfs_list
                    .collectFile() { item -> [ 'gvcf.list', "${item.get(1)}" + '\n' ] }
                    .set { gvcf_list }
                run_genomics_db_import_update(gvcf_list, chroms_all, run_backup_genomic_db.out.backup_status)
            } else if ( !db.exists() ) {
                println "Could not find existing DB:${db} please specify the correct DB path."
                exit 1
            }
            break
            // =====
        case['no']:
            if ( db.exists() ) {
                println "DB exists:${db}. Please check if you want to create new or update existing DB. If so please delete or backup first."
                exit 1
            } else if ( !db.exists() ) {
                println "Creating new DB:${db}"
                samples_gvcfs_list
                    .collectFile() { item -> [ 'gvcf.list', "${item.get(1)}" + '\n' ] }
                    .set { gvcf_list }
                run_genomics_db_import_new(gvcf_list, chroms_all)
            }
            break
            // =====   
    }
}

// GENOME CALLING
workflow GENOME_CALLING {
    take:
    gvcf
    chroms_all
    
    main:
    run_genotype_gvcf_on_genome_gvcf(gvcf, chroms_all)
    run_genotype_gvcf_on_genome_gvcf.out.gg_vcf_set
        .flatten()
        .collect()
        .map { it -> [ it.findAll { it =~ '.vcf.gz$' }, it.findAll { it =~ '.vcf.gz.tbi$' } ] }
        .set { concat_ready }
    run_concat_vcf(concat_ready)
    run_vqsr_on_snps(run_concat_vcf.out.combined_calls)
    run_apply_vqsr_on_snps(run_vqsr_on_snps.out.snps_vqsr_recal)
    run_vqsr_on_indels(run_concat_vcf.out.combined_calls)
    run_apply_vqsr_on_indels(run_vqsr_on_indels.out.indel_vqsr_recal)
}

// PICK & CHOOSE WORKFLOW
workflow {
    switch (workflow) {
            // ===== MAIN WORKFLOWS
        case['align']:
            ALIGN(samples_align)
            break
            // =====
        case['generate-gvcfs']:
            GENERATE_GVCFS(samples_generate_gvcfs_nosex, samples_generate_gvcfs_males, samples_generate_gvcfs_females, chroms_auto, chroms_par)
            break
            // =====
        case['combine-gvcfs']:
            COMBINE_GVCFS(samples_gvcfs_list, chroms_all)
            break
            // =====
        case['genomics-db-import']:
            GENOMICS_DB_IMPORT(samples_gvcfs_list, chroms_all)
            break
            // =====
        case['genome-calling']:
            GENOME_CALLING(gvcf, chroms_all)
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
    }
}
