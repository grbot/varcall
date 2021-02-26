# Intro

The Nextflow script runs mitochondrial calling on BAMs using Mutect2. (Currently setup for b37, need to add b38 compatibility)

Please see `nextflow.conf` for GATK version and references databases used.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimmted text file and should be specified in `nextflow.config`.  SampleID and BAM columns are required.

- BAM column should contain the flll path to the BAM.
- All collumns not used in this step (Fastq1R1, FastqR2, gVCF) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | /pathto/A01.bam | . |


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

1. `sample_id.MT.bam` - BAM file only containing MT reads
1. `sample_id.MT.vcf.gz` - Raw calls 
1. `sample_id.MT.vcf.gz.tbi` - Index for above
1. `sample_id.MT.vcf.gz.stats` - Call stats for above
1. `sample_id.MT.filtered.vcf.gz` - VCF flagged for filtering
1. `sample_id.MT.filtered.vcf.gz.tbi` - Index for above
1. `sample_id.MT.filtered-pass.vcf.gz` - VCF with only sites marked as PASS (final VCF)
1. `sample_id.MT.filtered-pass.vcf.gz.tbi - Index of above

