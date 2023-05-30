# Intro

Three Nextflow scripts are provided to run GenotypeGVCFs and VQSR. They behave differently and depending on what is required the one would be more suitable than the other.

1) `main.nf` - Can run on both output from CombineGVCFs or GenomicsDBImport. Runs GenotypeGVCFs in parallel per chromosome and VQSR on the whole gnome. Suitable for <= 1000 samples.

2) `main.intervals.nf` - More appropriate for output from GenomicsDBImport. Runs GenotypeGVCFs in parallel per intervals and VQSR on the whole gnome. Suitable for <= 1000
samples.

3) `main.intervals.vqsr-parallel.nf` - More appropriate for output from GenomicsDBImport. Runs GenotypeGVCFs in parallel per intervals and VQSR in perallel per chromosome. Suitable for > 1000 samples.

Please see `nextflow.conf` for GATK version and references databases used. Path to combined gVCF file is also specified in nextflow config.

## To run

For each dataset
1) Modify your `nextflow.config` to read the `gvcf_file` or `db_path` and specify the output directory e.g. `out_dir = "/global/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-out"`

See two examples for runs from output from CombineGVCFs or GenomicsDBImport

- CombineGVCFs - `nextflow.config.NA12878`, `gvcf_file` is set, `db_import = "no" and `db_path` is empty
- GenomicsDBImport - `nextflow.intervals.v6.config`, `gvcf_file` is empty,  `db_import = "yes"` and `db_path is set.

2) Run the workflow
```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/NA12878-work -c /home/gerrit/code/recalling/genome-calling/nextflow.config.NA12878 /home/gerrit/code/recalling/genome-calling/main.nf -profile cbio -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```

## Output

The output directory will contain
1. A per autosome and X, Y and MT VCF file
1. A recalibrated SNP genome VCF file
1. A recalibrated SNP and INDEL genome VCF file

