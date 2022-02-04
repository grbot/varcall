# Intro

Per sample data from SAMRC are normally received as groups of data distinguised by barcode though on the same lane. For example sample `101` whould have the following datasets.

`V350028760_L03_13_1.fq.gz`, `V350028760_L03_13_2.fq.gz`,`V350028760_L03_14_1.fq.gz`,`V350028760_L03_14_2.fq.gz`,`V350028760_L03_15_1.fq.gz`,`V350028760_L03_15_2.fq.gz`,`V350028760_L03_16_1.fq.gz`,`V350028760_L03_16_2.fq.gz`

Four barcode sets per read direction. Forward reads end with `_1.fq.gz` and reverse reads with `_2.fq.gz`

These sets needs to be combined into one set for each read direction. The workflow does exactly that.

Please see `nextflow.conf` for settings required.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimited text file and should be specified in `nextflow.config`.  For the alignment run, SampleIDi and FastqDir are required.

- SampleID contains the sample id (form naming purposes).
- FastqDir contains the full path to the sample Fastq reads.

| SampleID | FastqDir |
| -------- | ------ | 
| A01      | /pathto/A01 | 

## Examples

Attached is `example.samplesheet.tsv` an example sheet.

## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-workdir -c /home/gerrit/projects/varcall/combine-samrc-fastq/nextflow.config.NA12878.b37 /home/gerrit/projects/varcall/combine-samrc-fastq/main.nf -with-report NA12878.report.html -with-timeline NA12878.timeline.html -profile wits_slurm -resume
```

## Output

The output directory will contain per sample directories. Each sample directory will contain

1. `NA12878_f1.fq.gz` - Combined forward reads
2. `NA12878_f2.fq.gz` - Combined reverse reads
