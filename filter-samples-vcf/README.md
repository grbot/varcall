# Intro

The nextflow script selects or removes a list of samples from a VCF. Samples are selected/removed on chromosome level and then later combined into a complete VCF.
 

## To run

For each dataset
1) Modify your `nextflow.config`. Specify `in_files` to where your VCF is located. Specify `filter_samples` to where the list of samples to be selected/removed are located (list most contain a sample id per line). Also specify the output directory in `out_dir`.
2) Run the workflow
```
nextflow run /home/gerrit/code/recalling/filter-samples-vcf/main.nf  -c /home/gerrit/code/recalling/filter-samples-vcf/nextflow.config
```

## Output

The output directory will contain
1. A VCF file for each chromosome with samples selected/removed as well as a combined VCF with samples selected/removed
