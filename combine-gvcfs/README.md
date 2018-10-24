# Intro

The Nextflow script runs CombineGVCFs on a list of gVCFS.

Please see `nextflow.conf` for GATK version and references databases used.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the CombineGVCFs run, SampleID and gVCFs columns are required.

- gVCFs column should contain the full path to the gVCF (whole genome). 
- All collumns not used in this step (Fastq1R1, FastqR2, BAM) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | F      | .       | .       | .   | /pathto/A01.g.vcf.gz |


## To run

For each dataset
1) Create your sample sheet. E.g. `sahgp.samplesheet.tsv` for the SAHGP dataset. 
2) Modify your `nextflow.config` to read the `sahgp.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/adme/datasets/sahgp/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/sahgp/nextflow-work -c nextflow.config main.nf -profile wits -with-report sahgp.report.html -with-trace -with-timeline sahgp.timeline.html -resume
```

## Output

The output directory will contain  a cohort directory e.g. `sahgp`. With the following files.

1. `cohort_id.g.vcf.gz` and `cohort_id.g.vcf.gz.tbi` - the combined gVCF files. 
