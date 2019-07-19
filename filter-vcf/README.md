# Intro

The Nextflow script does two types of filtering on a VCF file. The filtering is for example relelvant on a final VCF file that was genotyped and VQSRed on short SNPs and INDELs.

1) It creates a filtered VCF file that contains only the SNP, INDELs, and MIXED sites that passed filtering. This is where the FILTER column is set to PASS. MIXED sites are sites that contain a combination of short SNPs and INDELS.
2) A second file is created that contains all sites that are flagged with `.` in their FILTER column.
 
Pleas  see `nextflow.conf` for GATK version and references databases used. Path to combined VCF file is also specified in nextflow config.

## To run

For each dataset
1) Modify your `nextflow.config` to read the `gvcf_file` and specify the output directory e.g. `out_dir = "/global/scratch/gerrit/projects/adme/datasets/NA12878/nextflow-out"`
2) Run the workflow
```
nextflow run /home/gerrit/code/recalling/filter-vcf/main.nf  -c /home/gerrit/code/recalling/filter-vcf/nextflow.config.NA12878
```

## Output

The output directory will contain
1. A VCF file containing high quality short SNPs and INDELS - `*.filter-pass.vcf.gz`
1. A VCF file containing MNPs - `*filter-other.vcf.gz`
