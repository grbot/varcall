# Intro

This nextflow pipeline converts the BAMs to CRAMs, runs samtools stats on the CRAM and creates an md5sum of the CRAM.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the BAM to CRAM conversion only the SampleID and BAM columns are required

- BAM column should contain the flll path to the AM.
- All collumns not used in this step (Gender, FastqR1, FastqR2, CRAM, gVCF) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM | CRAM | gVCF |
| -------- | ------ | ------- | ------- | --- | --- | --- |
| A01      | .      | .       | .       | /pathto/A01.bam | . | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv` for the GiaB dataset.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-work -c /home/gerrit/code/recalling/bam-to-cram/nextflow.NA12878.config /home/gerrit/code/recalling/bam-to-cram/main.nf -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```

