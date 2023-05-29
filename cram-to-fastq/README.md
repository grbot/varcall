# Intro

This nextflow pipeline converts the CRAMs or BAMs to Fastq. A forward and reverse read set will be created.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the BAM/CRAM to Fastq conversion only the SampleID and BAM columns are required

- BAM/CRAM column should contain the flll path to the BAM/CRAM.
- All collumns not used in this step (Gender, FastqR1, FastqR2, gVCF) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM/CRAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | /pathto/A01.cram | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `ggvp.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `ggvp.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /spaces/gerrit/ggvp/nextflow-work -c /home/gerrit/ggvp/varcall/cram-to-fastq/ggvp.nextflow.config /home/gerrit/ggvp/varcall/cram-to-fastq/main.nf -with-report ggvp.report.html -with-trace -with-timeline ggvp.timeline.html -profile ilifu -resume
```

