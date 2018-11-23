# Intro

This nextflow pipeline converts the CRAMs to Fastq. A forward and reverse read set will be created.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the CRAM to Fastq conversion only the SampleID and BAM columns are required

- BAM column should contain the flll path to the CRAM.
- All collumns not used in this step (Gender, FastqR1, FastqR2, gVCF) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | /pathto/A01.cram | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `awigen.samplesheet.tsv` for the 100 sample AwiGEN dataset.
2) Modify your `nextflow.config` to read the `awigen.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/awigen/nextflow-work -c /home/gerrit/projects/recalling/cram-to-fastq/nextflow.config.awigen /home/gerrit/projects/recalling/cram-to-fastq/main.nf -with-report awigen.report.html -with-trace -with-timeline awigen.timeline.html -profile wits_slurm -resume
```

