# Intro

The Nextflow script runs GATK's GenimicsImportDB on a set of gVCF files. Genome calling is then done on the generated database. There are two ways to create a database
* A new database can be created. In the nextflow config set `db_update = "no"` and specify a `db_path`.
* An existing database be updated. In the nextflow config set `db_update = "yes"` and specify a `db_path`.

## Sample sheet format

Below is the sample sheet format. The sample sheet should be a tab delimited text file and should be specified in `nextflow.config`.  For the genomics-db-import run, SampleID and gVCF columns are required.

- gVCF should contain the full path to the sample gVCF file.
- All collumns not used in this step (Gender, FastqR1, FastqR2, BAM) should be filled in with a "."


| SampleID | Gender | FastqR1 | FastqR2 | BAM | gVCF |
| -------- | ------ | ------- | ------- | --- | ---- |
| A01      | .      | .       | .       | /pathto/A01.g.vcf.gz | . |

## To run

For each dataset
1) Create your sample sheet. E.g. `subset.samplesheet.tsv`.
2) Modify your `nextflow.config` to read the `subset.samplesheet.tsv`i, specify the output directory and set the `db_update` and `db_path` settings.
3) Run the workflow e.g.
```
nextflow -log nextflow.log run -w /spaces/gerrit/projects/subset/nextflow-workdir -c /home/gerrit/projects/varcall/genomics-db-import/nextflow.config /home/gerrit/projects/varcall/genomics-db-import/main.nf -with-report subset.report.html -with-timeline subset.timeline.html -profile wits_slurm -resume
```

## Output

Each output directory will contain

1. `tool.gatk.version` - Tool version file (soft linked)
2. `gvcf.list` - List of gVCF files (soft linked)
3. `dbimport` - New or updated db (soft linked). Note the db path in the Nexflow config are also updated.
4. `dbimport.backup` - Backup of the old db if the db was updated (copy)

## Reasons for different `main.*`
- `main.nf` - GenomicsDBImport running in parallel on chromomes (Update DB available)
- `main.intervals.nf` - GenomicsDBImport running in parallel on intervals (No update on DB available. Need to add update DB from example in `main.nf` if there is a need.)

## get-fragmanets.py
For creating `intervals.list` for `main.intervals.nf`. Need to specify the sample size and a file with chromosome sizes. 

### b37
`cat /cbio/dbs/gatk/2.8/b37/human_g1k_v37_decoy.fasta.fai | cut -f 1,2 | grep -P "^1\t|^2\t|^3\t|^4\t|^5\t|^6\t|^7\t|^8\t|^9\t|^10\t|^11\t|^12\t|^13\t|^14\t|^15\t|^16\t|^17\t|^18\t|^19\t|^20\t|^21\t|^22\t|^X\t|^Y\t|^MT\t" > genome.sizes`

### b38
`cat /cbio/dbs/gatk/hg38/Homo_sapiens_assembly38.fasta.fai | cut -f 1,2 | grep -P "^chr1\t|^chr2\t|^chr3\t|^chr4\t|^chr5\t|^chr6\t|^chr7\t|^chr8\t|^chr9\t|^chr10\t|^chr11\t|^chr12\t|^chr13\t|^chr14\t|^chr15\t|^chr16\t|^chr17\t|^chr18\t|^chr19\t|^chr20\t|^chr21\t|^chr22\t|^chrX\t|^chrY\t|^chrM\t" > genome.sizes`

Then run `get-gragments.py`
```
get-fragments.py --sample-size sample-size --genome-size-file genome.sizes --intervals-file output-intervals-file
```
