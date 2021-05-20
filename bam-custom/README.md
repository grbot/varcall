# Intro

This nextflow pipeline runs `samtools view ${bam} ${region}` on BAMs or CRAMs.


## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  For the BAM/CRAM `samtools view ${bam} ${region}` only the SampleID and BAM/CRAM columns are required

- BAM/CRAM column should contain the file path to the BAM/CRAM.

| SampleID | BAM/CRAM |
| -------- |--------- |
| A01      | /pathto/A01.bam |

## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv` for the GiaB dataset.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv`
3) Run the workflow

```
nextflow -log nextflow.log run -w /scratch3/users/gerrit/scratch/NA12878/nextflow-work -c /users/gerrit/projects/refimpute/varcall/bam-custom/nextflow.NA12878.config /users/gerrit/projects/refimpute/varcall/bam-custom/main.nf -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -profile ilifu -resume
```

