# Intro

The Nextflow script calls on gene and genome level provided a combined gVCF file.

Please see `nextflow.conf` for GATK version and references databases used. Path to combined gVCF file is also specified in nextflow config.


## To run

For each dataset
1) Modify your `nextflow.config` to read the `gVCF_file` and specify the output directory e.g. `out_dir = "/spaces/gerrit/projects/adme/gene-calling/sahgp-baylor-agvp-sgdp-1kg_african/nextflow-out"`
2) Run the workflow
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/adme/gene-calling/sahgp-baylor-agvp-sgdp-1kg_african/nextflow-workdir -c /home/gerrit/projects/recalling/gene-calling/nextflow.config.sbas1.high /home/gerrit/projects/recalling/gene-calling/main.nf -with-report sbas1.report.html -with-trace -with-timeline sbas1.timeline.html -profile wits_slurm -resume
```

## Output

The output directory will contain a called file on gene level for each chromosome (`*_genes.vcf.gz`) and a called file on genome level for each chromosome (`*_genes.vcf.gz`).

