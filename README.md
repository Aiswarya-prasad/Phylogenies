# Phylogenies

Create phylogenies from MAG and isolate genomes. The aim is to see how MAGs compare to isolate genomes in a core genome phylogeny and to understand how much resolution they offer to study these communities. The MAGs need to have taxonomic annotations.

To get started there needs to be a config file and some metadata files with their paths specified in the config file. If the column order needs to be changed in the metadata files, the scripts, parsing functions and/or rules need to be modified accordingly.

The pipeline is currently being tested on one of the Engel lab workstations. However, since the conda environment is available, it can also be run in the cluster with minimal setting up (lines to provide slurm parameters must be uncommented and paths in config file should be changed).

## Choosing MAGs

The quality of MAGs mentioned below are according to this [paper](https://doi.org/10.1038/s41586-019-0965-1), mentioned in this [review](https://doi.org/10.1016/j.csbj.2021.11.028)

There are a total of 1029 redundant MAGs. drep clustered these into 86 unique clusters.

* 343 MAGs have are > 95% complete and have < 5% contamination. These would be considered **high quality** MAGs.
* 20 MAGs have are > 50% but < 95% complete and have < 10% but > 5% contamination. These would be considered **medium quality** MAGs.
* 877 MAGs have $completeness - 5 x Contamination > 50$ could be acceptable (mentioned in this [review](https://doi.org/10.1016/j.csbj.2021.11.028)).

152 MAGs will have to be dropped according to the most relaxed measures. Only high quality MAGs are used for making phylogenies.

Only **343** MAGs are at least medium quality out of 1029. So we dropped, 686 MAGs. These MAGs are listed in the metadata table along with the isolate genomes considered.

## Metadata and config files

The config file, `config.yaml` provides the following information:

* Project_identifier: **Name of the project**
* Project_path: **Absolute path to the project directory**
* Phylo_groups: **These can be any names such as phylotype names. It can also be an arbitrary name grouping the MAGs and genomes into phyologenetically similar groups for which orthologs must be inferred together and the tree plotted. These names should be the ones found in the column `Group` on the metadata file.**
* Genome_info_path: **Absolute path to metadata file.**


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
* Other columns can be added **after** this column without affecting the parsing of this file.

It is a `.tsv` file. In order to avoid parsing issues, export your metadata as a csv file and prepare the tsv file using the script at `scripts/csv_to_tsv.py`. Most other ways of making tsv files results in tabs being written as a variable number of spaces.

**Groups** are based on the family or order (annotated for the MAGs) for the most part. The isolate genomes are grouped based on their phylotype and/or the family they fall into. The names of the groups are:

* **bifido** all the bifidos

* **snod**: All the f__Neisseriaceae are snodgrasella so all in snod. Also f__Burkholderiaceae included with snod because same order o__Burkholderiales.

* **entero**: Everything with f__Enterobacteriaceae is entero. includes gilli and fper, hafnia (few) and londalia (few).

* **com**: All the f__Acetobacteraceae are com they are all commensalibacter.

* **api**: All the f__Weeksellaceae are apibacter so these are api.

* **lacto**: Includes all the order o__Lactobacillales which includes f__Lactobacillaceae and f__Enterococcaceae (few. foul-brood!).

* **bapis**: Includes all the bartonella. They fall under __Rhizobiaceae_A.

The other MAGs belonged to familes not mentioned above. They are not included for now.
**pseudo**: includes the following families : o__Pseudomonadales includes f__Pseudomonadaceae and f__Moraxellaceae (acinetobacter) just a few MAGs. For now, leaving them out.
**sphingo** includes about 6 MAGs. Can be combined with bapis or com as they all belong to the same class c__Alphaproteobacteria. For now, leaving them out.




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
