# Intro

This nextflow pipeline merges a multi lane sample BAM/CRAM into a single CRAM per sample.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the BAM to merging only the SampleID and BAM columns are required. Note the BAM column contains a directory to where the multi BAMs/CRAMs are located.


| SampleID | BAMdir |
| -------- | ------ |
| A01      | /pathto/A01 |


## To run

For each dataset
1) Create your sample sheet. E.g. `ggvp.samplesheet.tsv` for the Gambian Genome Variation Project dataset.
2) Modify your `nextflow.config` to read the `ggvp.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/ggvp/nextflow-work -c /home/gerrit/code/recalling/combine-lanes/nextflow.ggvp.config /home/gerrit/code/recalling/combine-lanes/main.nf -with-report ggvp.report.html -with-trace -with-timeline ggvp.timeline.html -resume
```

