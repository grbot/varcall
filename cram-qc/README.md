# Intro

This nextflow pipeline runs `samtools flagstat`, `samtools stats`, `mosdepth` and `verifybamid2` on BAMs/CRAMs and combines the reports in a MultiQC report.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  Only the SampleID and BAM columns are required

- BAM/CRAM column should contain the file path to the BAM/CRAM.
- All collumns not used in this step (Gender, FastqR1, FastqR2, gVCF) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM/CRAM | gVCF |
| -------- | ------ | ------- | ------- | -------- | ---- |
| A01      | .      | .       | .       | /pathto/A01.cram | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv` for the GiaB dataset.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /scratch2/users/gerrit/scratch/NA12878/nextlow-work -c nextflow.NA12878.config main.nf
```

