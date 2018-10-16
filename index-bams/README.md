# Intro

Some BAMs we received were not indexed. This nextflow pipeline indexes the BAMs (indexes are created in the location where the BAMs are stored).


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the BAM indexing only the SampleID and BAM columns are required

- BAM column should contain the flll path to the BAM.
- All collumns not used in this step (Gender, FastqR1, FastqR2, gVCF) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | /pathto/A01.bam | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `1kg_african.samplesheet.tsv` for the 1KG African dataset.
2) Modify your `nextflow.config` to read the `1kg_african.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/1kg_african/nextflow-work -c nextflow.config.1kg_african main.nf -profile wits -with-report 1kg_african.report.html -with-trace -with-timeline 1kg_african.timeline.html -resume
```

## Note
After the run the index files needs to be moved back to where the original BAMs are located.
