# Intro

The Nextflow script runs bwa-mem, Picard markduplicates, GATK BQSR, BAM to CRAM and samtools flagstat on a list of sample reads.

Please see `nextflow.conf` for bwa, samtools and GATK and Docker/Singularity versions and references databases used.

Sentieon is now included as an option for alignment. In config set `sentieon = true` to enable or `sentieon = false` to run normal BWA/GATK workflow. See `nextflow.config.NA12878.b38` for additional settings that are needed for a Sentieon run.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimited text file and should be specified in `nextflow.config`.  For the alignment run, SampleID, FastqR1, FastqR2, Flowcell and Lane columns are required.

- FastqR1 and FastqR2 should contain the full path to the sample Fastq reads.
- Flowcell should contain the flowcell id
- Lane should contain the lane id. **Note** if the lane id = 0, we assume a sample that was sequenced across lanes. The BAM/CRAM naming in this case would be e.g. `sample_id.bam` where in cases where you have multiple lane sequenced data your BAM/CRAM naming would be `sample_id-flowcell.id.bam` 
- All collumns not used in this step (Gender, BAM, gVCF) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | Flowcell | Lane | BAM | gVCF |
| -------- | ------ | ------- | ------- | -------- | ---- | --- | --- |
| A01      | .      | /pathto/A01_R1.fastq.gz       | /pathto/A01_R2.fastq.gz  | HHTN2BBXX |  6 | .  | . |

## Examples

Attached is `NA12878.samplesheet.tsv` an example sheet and `nextflow.config.NA12878.b37` and `nextflow.config.NA12878.b38` to show configuration settings for different reference genomes.

## To run

For each dataset
1) Create your sample sheet. E.g. `NA12878.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `NA12878.samplesheet.tsv` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-out"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/1kg/datasets/NA12878/nextflow-workdir -c /home/gerrit/projects/varcall/align/nextflow.config.NA12878.b37 /home/gerrit/projects/varcall/align/main.nf -profile wits_slurm -resume
```

## Output

The output directory will contain per sample directories. Each sample directory will contain (Naming will differ depending on the type of lane sequencing that were done, see discussion above on lane)

1. `NA12878.bam` - BAM file after BWA alignment (soft linked)
2. `NA12878.md.bam` - BAM file after running MarkDuplicates (soft linked)
3. `NA12878.md.recal.bam` - BAM file after running GATK BQSR (final BAM) (soft linked)
4. `NA12878.md.recal.cram` - CRAM file converted from `NA12878.md.recal.bam` (copy)
5. `NA12878.md.recal.cram.flagstat` - Samtools stats on `NA12878.md.recal.cram` (copy)
6. `NA12878.md.recal.cram.md5sum` - Md5sum of `NA12878.md.recal.cram` (copy)
