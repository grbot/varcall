# Intro

The Nextflow script runs GATKs ValidateVaraints on gVCFs.

Please see `nextflow.conf` for GATK version and references databases used.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the run, SampleID, Gender and gVCF columns are required.

- Gender column should contain M for male and F for Female.
- gVCF column should contain the flll path to the gVCF.
- All collumns not used in this step (Fastq1R1, FastqR2, BAM) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | F      | .       | .       | . | /pathto/gvcf |


## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.txt` for the NA12878 dataset.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/adme/datasets/sahgp/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/NA12878/nextflow-work -c nextflow.config main.nf -profile wits -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```

## Output

The output directory will contain per sample directories. Each sample directory will contain

1. A `sample_id.validatevariants` file
