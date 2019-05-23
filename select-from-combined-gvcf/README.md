# Intro

The Nextflow script removes samples from a gvcf previously created by CombineGVCFs.

Please see `nextflow.conf` for GATK version and references databases used.

## To run

For each dataset
1) Modify your nextflow.config to read the `combined_gvcf`, specify the path to the list of samples to be removed (`exclude_samples`) and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/adme/combine-gvcfs/sahgp-baylor-sgdp-1kg_african-caroline-awigen/nextflow-out/"`
3) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/datasets/sbs1ca/nextflow-work -c nextflow.config main.nf -profile wits -with-report sbs1ca.report.html -with-trace -with-timeline sbs1ca.timeline.html -resume
```

## Output

The output directory will contain a gvcf file with the samples removed. E.g. `sbs1ca.g.vcf.gz`.

