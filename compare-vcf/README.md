# Intro

The workflow calulates the overlap and intersections between VCFs based on a common list of samples. It also calculates the novelty by comparing the overlap and intersections with dbSNP. 

## To run

Specify the two VCFs to compare: `vcf_1` and `vcf_2`. Specify the common set of samples: `sample_list`. Specify the dbSNP VCF: `dbsnp`

See example `nextflow.config` 


Run `nextflow -log nextflow.log run /users/gerrit/projects/refimpute/varcall/compare-vcf/main.nf -w /scratch3/users/gerrit/projects/refimpute/v6.high-cov.sentieon/compare-vcf-work -c /users/gerrit/projects/refimpute/varcall/compare-vcf/nextflow.v6.high-cov.sentieon.config -profile ilifu -resume
`

## Output

Final output is in `final.stats.txt`

`NRD SNPs (%)` - Non-reference discordance % between SNPs for samples in `vcf_1` and `vcf_2`
`NRD INDELs (%)` - Non-reference discordance % between INDELs for samples in `vcf_1` and `vcf_2`
`Records unique to 1st VCF` - Number of records unique to `vcf_1`
`Records unique to 2nd VCF` - Number of records unique to `vcf_2`
`Records that overlap` - Number of records that overlap between `vcf_1` and `vcf_2`
`Records in 0000.sorted.vcf.gz_vs_dbsSNP.bcftools.stats that are not novel` - Number of records that are unique to `vcf_1` that are found in dbSNP
`Records in 0001.sorted.vcf.gz_vs_dbsSNP.bcftools.stats that are not novel` - Number of records that are unique to `vcf_2` that are found in dbSNP
`Records in 0002.sorted.vcf.gz_vs_dbsSNP.bcftools.stats that are not novel` - Number of records that are in overlap with `vcf_1` and `vcf_2` that are found in dbSNP
`Records in 0003.sorted.vcf.gz_vs_dbsSNP.bcftools.stats that are not novel` - Number of records that are in overlap with `vcf_1` and `vcf_2` that are found in dbSNP (should be same as one above)
