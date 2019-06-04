# Intro

The Nextflow scripts in this repository were mostly used for a specific H3Africa project. The project contained datasets at different stages e.g Fastq, BAM and gVCF so it made sense just to process the set at stages where necessary.

## To run

Each folder contains part of the pipeline or separate scripts to prepare the data for a specific part of the pipeline or get it ready for archival purposes.


### Main pipeline
* **align** - Need to provide a sample sheet containing paths to forward and reverse reads. BWA, Picard MarkDuplicates and GATK BaseRecalibrator are run on the samples.
* **generate-gvcf** - Need to provide a sample sheet with paths to the BAMs. GATKs HaplotypeCaller was run in gVCF mode on the samples.
* **combine-gvcfs** - Need to provide a sample sheet with paths to the gVCFs. GATKs CombineGVCFs were run on all samples.
* **gene-calling** - Need to provide the path to the combined gVCF. GenotypeGVCFs are run jointly on all samples on gene regions only.
* **genome-calling** - Need to provide the path to the combined gVCF. GenotypeGVCFs are run jointly on all samples on full genome. VQSR is also applied to chromosome 1 to 22, X, Y and MT. GenotypeGVCFs are run on chromosome level where SNP and INDEL VQSR are run on genome level (**according to GATKs best practices**).
* **genome-calling-vqsr-per-chromosome** Need to provide the path to the combined gVCF. GenotypeGVCFs are run jointly on all samples on full genome. VQSR is also applied to chromosome 1 to 22. Calling is only done for X, Y and MT in a separate Nextflow script. GenotypeGVCFs are run on chromosome level where SNP and INDEL VQSR are run on chromosome level.


### Separate pipelines
* **bam-to-cram** - Need to provide a sample sheet with paths to the BAMs. BAMs are converted to CRAM (v3), indexed, stats are calculated and md5sums are generated.
* **cram-to-fastq** - Need to provide a sample sheet with paths to the CRAMs. CRAMs are converted to Fastq (forward and reverse pair).`
* **index-bams** - Need to provide a sample sheet with paths to the BAMs. BAMs are indexed.`
* **select-from-combined-gvcf** - Need to provide the path to the combined GVCF and list of sample ids to be removed. Sample ids are removed from the combined gVCF and a new combined gVCF are generated.`

## Note
* The main pipeline can be followed for a single sample or joint calling sample. For a single sample the `combine-gvcf` can be skipped and the `gene-calling` and `genome-calling` configs can point to single sample `.g.vcf`. 

## Example
To get test data that van be used from the alignment to calling and VQSR phase download the GiaB dataset here
```
wget -O NA12878_R1.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget -O NA12878_R2.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

## To do
1. Pull all software as Docker containers (tabix, samtools, GATK). 
1. Currently setup for `b37`. Setup a swith in nextflow to allow for processing against build `b38`. This would required switching between chromosome names and reference databases
1. All output alignments in the `align` step should be in CRAM
