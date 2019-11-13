# Intro

The repos contains individual pipelines to process large cohorts of human genome samples through alignment, calling and joint calling.

## To run

Each folder contains part of the pipeline or scripts to prepare the data for a specific part of the pipeline.


### Main pipeline
* **align** - Need to provide a sample sheet containing paths to forward and reverse reads. BWA, Picard MarkDuplicates and GATK BaseRecalibrator are run on the samples.
* **generate-gvcf** - Need to provide a sample sheet with paths to the BAMs. GATKs HaplotypeCaller is run in gVCF mode on the samples.
* **combine-gvcfs** - Need to provide a sample sheet with paths to the gVCFs. GATKs CombineGVCFs is run on all samples.
* **genomics-db-import** - Need to provide a sample sheet with paths to the gVCFs. GATKs GenomicsDBImport is run on all samples. A new database can be created or an existing database can be updated.
* **genome-calling** - Need to provide the path to the combined gVCF. GenotypeGVCFs are run jointly on all samples on full genome. VQSR is also applied to chromosome 1 to 22, X, Y and MT. GenotypeGVCFs are run on chromosome level where SNP and INDEL VQSR are run on genome level (**according to GATKs best practices**).

**Note**: Can use **combine-gvcfs** or **genomics-db-import** for combining the gVCFs. **genomics-db-import** is however the latest method and allows for updating of an existing database / combined set.

### Separate pipelines
* **bam-to-cram** - Need to provide a sample sheet with paths to the BAMs. BAMs are converted to CRAM (v3), indexed, stats are calculated and md5sums are generated.
* **cram-to-fastq** - Need to provide a sample sheet with paths to the CRAMs. CRAMs are converted to Fastq (forward and reverse pair).
* **index-bams** - Need to provide a sample sheet with paths to the BAMs. BAMs are indexed.
* **filter-vcf** - Filter final VCFs based on the PASS flag in the FILTER column.
* **bam-flagstat** - Need to provide a sample sheet with paths to the BAMs/CRAMS. Flagstat reporsts are created.

## Note
* The main pipeline can be followed for a single sample or joint calling sample. For a single sample the `combine-gvcfs` or `genomics-db-import` can be skipped and the `genome-calling` configs can point to single sample `.g.vcf`.

## Example
To get test data that van be used from the alignment to calling and VQSR phase download the GiaB dataset here
```
wget -O NA12878_R1.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget -O NA12878_R2.fastq.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

## Tested on Ilifu SLURM cluster
- [x] **align**
- [x] **generate-gvcf**
- [x] **combine-gvcf**
- [x] **genomics-db-import**
- [x] **genome-calling**
- [x] **cram-to-fastq**
- [x] **bam-to-cram**
- [x] **index-bams**
- [x] **filter-vcf**
- [x] **bam-flagstat**
