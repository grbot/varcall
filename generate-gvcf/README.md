# Intro

The Nextflow script run HaplotypeCaller on BAMs. It has been tested on the Wits cluster and contains a profile setting for Wits. For CBIO we just need to add another profile. For the processing of 1KG on CHPC the code would need to be rewritten.

Please see `nextflow.conf` for GATK version and references databases used.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the HaplotypeCaller run, SampleID, Gender and BAM columns are required.

- Gender column should contain M for male and F for Female.
- BAM column should contain the flll path to the BAM.
- All collumns not used in this step (Fastq1R1, FastqR2, gVCF) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | F      | .       | .       | /pathto/A01.bam | . |


## To run

For each dataset
1) Create your sample sheet. E.g. `sahgp.samplesheet.txt` for the SAHGP dataset.
2) Modify your `nextflow.config` to read the `sahgp.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/adme/datasets/sahgp/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/sahgp/nextflow-work -c nextflow.config main.nf -profile wits -with-report sahgp.report.html -with-trace -with-timeline sahgp.timeline.html -resume
```

## Output

The output directory will contain per sample directories. Each sample directory will contain

1) `sample_id.[1-22].g.vcf.gz` files
2) `sample_id.MT.g.vcf.gz`
3)
 a) Female samples, a `sample_id.X.g.vcf.gz`
 b) Male samples, a  `sample_id.X_PAR1.g.vcf.gz,sample_id.X_PAR2.g.vcf.gz,sample_id.X_nonPAR.g.vcf.gz,sample_id.Y_PAR1.g.vcf.gz,sample_id.Y_nonPAR.g.vcf.gz`
4) A `sample_id.g.vcf.gz` that is a combined, sorted and indexed file of the files above.  Would be easiest to use this file for downstream processing.
