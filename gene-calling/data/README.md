* Get all genes from [http://grch37.ensembl.org/biomart/martview/](http://grch37.ensembl.org/biomart/martview/)
* Get chromosome, gene start bp, and gene end bp: `awk 'NR>1 {print $1"\t"$2"\t"$3"\n"}' ensembl_genes_b37_20Nov2018_annotation.bed > ensembl_genes_b37_20Nov2018_annotation_final.bed`
* Sort bed by chromosome and gene start bp: `sort -k1,1 -k2,2n ensembl_genes_b37_20Nov2018_annotation_final.bed  | uniq > ensembl_genes_b37_20Nov2018_annotation_final.sorted.bed`
* Merge overlapping poisiton using bedtools: `bedtools merge -i ensembl_genes_b37_20Nov2018_annotation_final.sorted.bed > ensembl_genes_b37_20Nov2018_annotation_final.sorted.merged.bed`
* Add 25kb flanking regions up and down stream: `awk -F"\t" '{print $1"\t"$2-25"\t"$3+25"}' ensembl_genes_b37_20Nov2018_annotation_final.sorted.merged.bed > ensembl_genes_b37_20Nov2018_annotation_final.sorted.merged.slop25kb.bed`
