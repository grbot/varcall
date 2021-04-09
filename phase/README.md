# Intro

The Nextflow script uses eagle2 to phase chr1 -> 22. The code works but is still in an experimental phase in terms of structure. Currently phasing is done without a reference and if phasing with a reference wants to be done it just needs to be uncommented.

Please see `nextflow.conf` for tool version and references databases used. 

## To run

For each dataset
1) Modify your `nextflow.config` to read the `vcf` and specify the output directory e.g. `out_dir = "/scratch3/users/gerrit/projects/adme/phase.without-reference-biallelic-split"`
2) Run the workflow
```
nextflow -log nextflow.log run -c /users/gerrit/projects/refimpute/varcall/phase/nextflow.config.adme /users/gerrit/projects/refimpute/varcall/phase/main.nf -w /scratch3/users/gerrit/projects/adme/nextflow-work-phase  --with-report adme.phase.report.html --with-timeline adme.phase.timeline.html -profile ilifu -resume
```

## Output

The output directory will contain
1. Per chromosme phased data - *.phased.chr*.vcf.gz
1. Combined phased data - *.phased.vcf.gz
