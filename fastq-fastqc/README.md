# Intro

The Nextflow script runs FastQC on a list of sample reads and combine the results in a MultiQC report.

Please see `nextflow.conf` for FastQC and MultiQC versions used.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimited text file and should be specified in `nextflow.config`.  For the FastqQC run, SampleID, FastqR1 and FastqR2 columns are required.

- FastqR1 and FastqR2 should contain the full path to the sample Fastq reads.
- All collumns not used in this step (e.g. Gender, Flowcell, Lane, BAM, gVCF) should be filled in with a "." or should not be included


| SampleID | Gender | FastqR1 | FastqR2 | Flowcell | Lane | BAM | gVCF |
| -------- | ------ | ------- | ------- | -------- | ---- | --- | --- |
| A01      | .      | /pathto/A01_R1.fastq.gz       | /pathto/A01_R2.fastq.gz  | . |  . | .  | . |

## Examples

Attached is `NA12878.samplesheet.tsv` an example sheet and `nextflow.config.NA12878.b37` and `nextflow.config.NA12878.b38` to show configuration settings for different reference genomes.

## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-workdir -c /home/gerrit/projects/varcall/align/nextflow.config.NA12878.b37 /home/gerrit/projects/varcall/align/main.nf -with-report NA12878.report.html -with-timeline NA12878.timeline.html -profile wits_slurm -resume
```

## Output

The output directory will contain per sample directories with FastQC files a MultiQC report containing FastQC info from all the samples.
