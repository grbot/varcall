# Intro

The Nextflow script does Genotype Refinement (runs CalculateGenotypePosteriors, VariantFiltration and VariantAnnotator ) on the whole genome (chr1 -> 22, X, Y and MT).

Please see `nextflow.conf` for GATK version and references databases used. Path to VQSR filtered VCF file is also specified in nextflow config.

## To run

For each dataset
1) Modify your `nextflow.config` to read the `vcf` and specify the output directory e.g. `out_dir = "/global/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-out"`
2) Run the workflow
```
nextflow -log nextflow.log run -w /global5/scratch/gerrit/projects/adme/datasets/NA12878-work -c /home/gerrit/code/varcall/genotype-refinement/nextflow.config.NA12878 /home/gerrit/code/varcall/genotype-refinement/main.nf -profile cbio -with-report NA12878.report.html -with-trace -with-timeline NA12878.timeline.html -resume
```

## Output

The output directory will contain
1. CalculateGenotypePosteriors VCFs - *.cgp.vf.*.vcf.gz*
1. VariantFiltration VCFs - *.cgp.vf.*.vcf.gz*
1. VariantAnnotator VCFs - *.cgp.vf.va.*.vcf.gz*
1. Genome combined processed VCF - *.cgp.vf.va.vcf.gz*
