# INTRODUCTION
**_Work in progress_**

## WORKFLOW PARAMETERS:
```
## INPUTS AND OUTPUTS
params.sample_sheet       = "/home/phelelani/nf-workflows/varcall/data/NA12878.samplesheet.tsv"
params.project_name       = "NA12878"
params.outdir             = "/home/phelelani/nf-workflows/varcall/WF-GENERATE_GVCFS"
params.sample_coverage    = "high"

## REFERENCE GENOME BUILD
params.build              = "b37"
params.ref                = "/home/phelelani/data/gatk-b37-bundle/human_g1k_v37_decoy.fasta" 

## WORKFLOW & ANALYSIS TYPE
params.workflow           = null
params.type               = "wgs"

## RESOURCES
params.dbsnp              = "/home/phelelani/data/gatk-b37-bundle/dbsnp_138.b37.vcf"
params.known_indels_1     = "/home/phelelani/data/gatk-b37-bundle/1000G_phase1.indels.b37.vcf"
params.known_indels_2     = "/home/phelelani/data/gatk-b37-bundle/Mills_and_1000G_gold_standard.indels.b37.vcf"
params.hapmap             = "/home/phelelani/data/gatk-b37-bundle/hapmap_3.3.b37.vcf"
params.omni               = "/home/phelelani/data/gatk-b37-bundle/1000G_omni2.5.b37.vcf"
params.phase1_snps        = "/home/phelelani/data/gatk-b37-bundle/1000G_phase1.snps.high_confidence.b37.vcf"
params.golden_indels      = "/home/phelelani/data/gatk-b37-bundle/Mills_and_1000G_gold_standard.indels.b37.vcf"

## COMPUTATION RESOURCES
params.memory_per_thread  = "500M"
params.bwa_threads        = "32"

## FOR WHOLE EXOME SEQUENCING
params.target_regions     = ""
```

## NB: SAMPLE SHEET!
| SampleID | Gender | FastqR1 | FastqR2 | Flowcell | Lane | BAM | gVCF |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| NA12878 |	. | NA12878\_R1.fastq.gz | NA12878\_R2.fastq.gz | HHTN2BBXX	 | 6 | NA12878.md.recal.cram | . |
| NA12878\_M | M | NA12878\_M\_R1.fastq.gz | NA12878\_M\_R2.fastq.gz | HHTN2BBXX | 6 | NA12878\_M.md.recal.cram | . |
| NA12878\_F | F | NA12878\_F\_R1.fastq.gz | NA12878\_F\_R2.fastq.gz | HHTN2BBXX | 6 | NA12878\_F.md.recal.cram | . |

## VARIANT CALLING WORKFLOWS
### 1. Alignment to Reference Genome
```
nextflow run main.nf -profile wits --workflow align --sample_sheet <samplesheet.tsv>
```

### 2. Generate GVCFs
```
nextflow run main.nf -profile wits --workflow generate-gvcfs --sample_sheet <samplesheet.tsv>
```

### 3. Combine GVCFs
```
nextflow run main.nf -profile wits --workflow combine-gvcfs --sample_sheet <samplesheet.tsv>
```

### 4. Genomics DB Import
```
nextflow run main.nf -profile wits --workflow genomics-db-import --sample_sheet <samplesheet.tsv>
```

<!-- The repos contains individual wrokflows to process cohorts of human genome samples through alignment, calling, joint calling and variant quality score recalibration. Some additional workflows for pre and post BAM/gVCF manipulation/checking are also included. -->

<!-- ## To run -->

<!-- Each folder contains a workflow which is part of the main workflow or a separate workflow to check/prepare the data for a specific part of the main workflow. -->


<!-- ### Main workflow -->
<!-- * **align** - Need to provide a sample sheet containing paths to forward and reverse reads. BWA, Picard MarkDuplicates and GATK BaseRecalibrator are run on the samples. -->
<!-- * **generate-gvcf** - Need to provide a sample sheet with paths to the BAMs. GATKs HaplotypeCaller is run in gVCF mode on the samples. -->
<!-- * **combine-gvcfs** - Need to provide a sample sheet with paths to the gVCFs. GATKs CombineGVCFs is run on all samples. -->
<!-- * **genomics-db-import** - Need to provide a sample sheet with paths to the gVCFs. GATKs GenomicsDBImport is run on all samples. A new database can be created or an existing database can be updated. -->
<!-- * **genome-calling** - Need to provide the path to the combined gVCF. GenotypeGVCFs are run jointly on all samples on full genome. VQSR is also applied to chromosome 1 to 22, X, Y and MT. GenotypeGVCFs are run on chromosome level where SNP and INDEL VQSR are run on genome level (**according to GATKs best practices**). -->
<!-- * **filter-vcf** - Filter final VCFs based on the PASS flag in the FILTER column. Needed after VQSR. -->

<!-- **Note**: Can use **combine-gvcfs** or **genomics-db-import** for combining the gVCFs. **genomics-db-import** is however the latest method and allows for updating of an existing database / combined set. -->

<!-- ### Separate workflows -->
<!-- * **bam-to-cram** - Need to provide a sample sheet with paths to the BAMs. BAMs are converted to CRAM (v3), indexed, stats are calculated and md5sums are generated. -->
<!-- * **cram-to-bam** - Need to provide a sample sheet with paths to the CRAMs. CRAMs are converted to BAM, indexed, stats are calculated and md5sums are generated. -->
<!-- * **cram-to-fastq** - Need to provide a sample sheet with paths to the CRAMs. CRAMs are converted to Fastq (forward and reverse pair). -->
<!-- * **index-bams** - Need to provide a sample sheet with paths to the BAMs. BAMs are indexed. -->
<!-- * **index-vcf** - Need to provide a sample sheet with paths to the VCFs/gVCFs. VCFs/gVCFs are indexed. -->
<!-- * **bam-flagstat** - Need to provide a sample sheet with paths to the BAMs/CRAMs. Flagstat reporsts are created. -->
<!-- * **mt-calling** - Need to provide a samples sheet with paths to BAMs/CRAMs. Mitochondrial variant calling are doen with the Mutect2 pipeline. -->
<!-- * **combine-lanes** - Need to provide a samplesheet to directories with the multi-lane BAMs. BAMs are merged and indexed. -->
<!-- * **validate-gvcf** - Need to provide a samplesheet with path to gVCFs and gender info. gVCF validation are done on chromosome level using GATK's ValidateVariants.  -->
<!-- * **genotype-refinement** - Need to provide a path to the VQSRed VCF. Genotype Refinment are done on chromosome level using GATK's CalculateGenotypePosteriors, VariantFiltration and VariantAnnotator.  -->

<!-- ## Note -->
<!-- * The main workflow can be followed for a single sample or joint calling sample. For a single sample the **combine-gvcfs** or **genomics-db-import** can be skipped and the **genome-calling** configs can point to single sample `.g.vcf`. -->

<!-- ## Example -->
<!-- The GiaB dataset can be downloaded and used for testing -->
<!-- ``` -->
<!-- wget -O NA12878_R1.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz -->
<!-- wget -O NA12878_R2.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz -->
<!-- ``` -->
