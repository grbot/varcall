# Intro

The Nextflow script runs CombineGVCFs on a list of gVCFS.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the CombineGVCFs run, SampleID and gVCFs columns are required.

- gVCFs column should contain the full path to the gVCF (whole genome). 
- All collumns not used in this step (Fastq1R1, FastqR2, BAM) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | .   | /pathto/A01.g.vcf.gz |


## To run

For each dataset
1) Create your sample sheet. E.g. `samplesheet.tsv`. 
2) Modify your `nextflow.config` to read the `samplesheet.tsv`, specify the output directory and the reference build.
3) Run the workflow
```
nextflow -log nextflow.log run -w /work -c nextflow.config main.nf -profile ilifu -with-report report.html -with-trace -with-timeline timeline.html -resume
```

## Output

The output directory will contain  a directory with the following files.

1. `cohort_id.g.vcf.gz` and `cohort_id.g.vcf.gz.tbi` - the combined gVCF files. 
