# Intro

This nextflow pipeline indexes the VCFs/gVCFs.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the VCF/gVCF indexing only the SampleID and gVCF columns are required

- gVCFF column should contain the flll path to the BAM.
- All collumns not used in this step (Gender, FastqR1, FastqR2,i BAM) should be filled in with a "." 


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | . | /pathto/A01.gvcf |


## To run

For each dataset
1) Create your sample sheet. E.g. `samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w work -c config -profile ilifu -with-report report.html -with-trace -with-timeline timeline.html -resume
```

## Note
After the run the index files needs to be moved back to where the original VCFs/gVCFs are located.
