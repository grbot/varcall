# Intro

The Nextflow script calls on genome level (runs GenotypeGVCFs) and does calling of VSQR SNPs and INDELs per genome not chromosome. Currently only doing it for chr1 -> 22, X, Y and MT. X, Y and MT will probably not be recalibrated.

Please see `nextflow.conf` for GATK version and references databases used. Path to combined gVCF file is also specified in nextflow config.

## To run

For each dataset
1) Modify your `nextflow.config` to read the `gvcf_file` and specify the output directory e.g. `out_dir = "/global/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-out"`
2) Run the workflow
```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/NA12878-work -c /home/gerrit/code/recalling/genome-calling/nextflow.config.NA12878 /home/gerrit/code/recalling/genome-calling/main.nf -profile cbio -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```

## Output

The output directory will contain
1. A per chromosome VCF file (X, Y and MT will only have these files)
1. A recalibrated SNP per chromsome VCF file
1. A recalibrated SNP and INDEL per chromsome VCF file
1. A combined genome recalibrated SNP and INDEL VCF file

