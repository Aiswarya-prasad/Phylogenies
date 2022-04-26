# Phylogenies

Create phylogenies from MAG and isolate genomes. The aim is to see how MAGs compare to isolate genomes in a core genome phylogeny and to understand how much resolution they offer to study these communities. The MAGs need to have taxonomic annotations.

To get started there needs to be a config file and some metadata files with their paths specified in the config file. If the column order needs to be changed in the metadata files, the scripts, parsing functions and/or rules need to be modified accordingly.

The pipeline is currently being tested on one of the Engel lab workstations. However, since the conda environment is available, it can also be run in the cluster with minimal setting up (lines to provide slurm parameters must be uncommented and paths in config file should be changed).

## Metadata and config files

The config file, `config.yaml` provides the following information:

* Project_identifier: **Name of the project**
* Project_path: **Absolute path to the project directory**
* Phylo_groups: **These can be any names such as phylotype names. It can also be an arbitrary name grouping the MAGs and genomes into phyologenetically similar groups for which orthologs must be inferred together and the tree plotted. These names should be the ones found in the column `Group` on the metadata file.**
* Genome_info_path: **Absolute path to metadata file.**

The metadata file, `GenomeInfo.csv` contains the following columns:

* ID: Unique name for each genome. Could be same as locus tag or strain name for example. **(Required)**
* Accession: Accession number for the database that it came from.
* Locus_tag:
* Strain_name: **(Required)**
* Phylotype: Phylotype if known from earlier studies.
* SDP: If known from previous calculations.
* Species: As annotated for MAGs.
* Host: Source organism.
* Study:
* Origin: Location it was sampled from.
* Source_database
* Cluster: Relevant for MAGs if they were dereplicated and multiple MAGs per cluster are used.
* Group: **(Required)**

There are no spaces between columns. Entries that contain commas are written within quotes

## Rules

### Download genomes

All the MAGs and isolate genomes listed in the metadata file should be placed in the directory within the project directory named `00_genomes`. If any isolate genomes are missing they will be downloaded. However this is not too robust because of the way the link to the genome on the NCBI ftp server is obtained especially if the strain name is not unique although this works for most geneomes.

`genbank_path` is parsed from the downloaded `assembly_summary_genbank.txt` file. However, to avoid having multiple urls if there are
multiple genomes for the strain it is sorted by date column (format "yyyy/mm/dd") in reversed order (so most recent first) and the top entry is used.

### Annotate genomes

We re-annotate them all in this pipeline (even though the annotations should be available online for the isolate genomes) to maintain consistency.

use force because snakemake creates output dirs and this confuses it

### Orthofinder



## Scripts

Putting together scripts borrowed from Garance into individual scripts or `Snakemake` rules.
