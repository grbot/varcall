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
