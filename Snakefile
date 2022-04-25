import os

configfile: "config.yaml"

ProjectIdentifier = config["ProjectIdentifier"]
ProjectPath = config["ProjectPath"]
GenomeInfo = config["Genome_info_path"]

"""
define functions

- read genome dictonary as needed
"""

def convertToMb(string):
    """
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    """
    if string.endswith("G"):
        number = int(string.split("G")[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    """
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    """
    days = string.split("-")[0]
    hrs = string.split("-")[1].split(":")[0]
    min = string.split("-")[1].split(":")[1]
    sec = string.split("-")[1].split(":")[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

"""
rule all:
    needs just the final trees!

Start with genome nucleotide sequences
groups for Shinilab MAGs is based on family
(and lower lowest non-ambiguous choice picked) assigned by GTDBtk
the group columns are arbitrarily based on the gtdbtk results.
    - lacto - firm5 firm4 lkun and all the f__Lactobacillaceae
    - gillif - gilliamella and fper
    - bapis, bifido, bom, com, api, snod
    for now ignore:
    sphingo, syntro, moraxella,
    cutibac, burkhold, pseudo, entero, rhodocy, enteroc

rule annotate (prokka)
    re-annotating to maintain consisitency
    check Garance's doAnnotation_${Genus}.sh

rule prep for Orthofinder
    organise faa files by phylotype

rule Orthofinder
    run for each group (phylotype)
    check Garance's OrthoFinder.sh
    ├─ OrthoFinder :
	│	└─ Results_${Genus}
	│		├─ Single_Copy_Orthologue_Sequences : directory of interest

rule multiple sequence alignement on Orthogroups
    check Garances's AlignOG.sh

rule pruned alignments
    check Garance's ../scripts/SAGE_FetchAndConcat2.py or execPythonScript.sh
    it makes
    │	├─ CoreGeneAlignment.fasta : concatenated OG per species -> to use for phylogeny
    │	└─ OGxx_aligned_prined.fasta files : single OG aligned and pruned

rule execIQTREE

    check Garance's execIQTREE.sh
    └─ phylogeny :
		│  cmd ```sbatch execIQTREE.sh ```
		│
		├─ WGP_${Familiy}.contree : consensus tree with renamed samples.
		│  not only locus but complete species name + strain
		│
		└─ .contree ; .treefile .log ...
    #!/bin/bash
    ## $1 path to OG alignments

    for OG in $(ls $1) ; do
    	iqtree -s $1$OG -st AA -alrt 1000 -b 500 -seed 1234 -m GTR20 -redo -ntmax 20 -o Swig_F0424.faa,Psui_DSM24744.faa,Aaes_DSM22365.faa
    done

scripts by Garance to refer to:
+ scripts
│
├── AlignOG.sh : to align sequences
│	│
│	└─ parameters :
│		--amino : use amino-acids
│		--inputorder
│		--localpair
│		--maxiterate 1000
│
├── doAnnotation_${Genus}.sh
│	│
│	├─ arrays :
│	│	- array SAMPLE : .fna file names
│	│	- array LocusTag : locus tag associated to .fna files (thy have to be in the same order as .fna file names)
│	│
│	└─ parameters :
│
│		--compliant
│		--proteins proteins.faa
│		--evalue 0.01
│
├── execIQTREE.sh :
│	│
│	├─ array :
│	│	- Genus : Genus names corresponding to the ${Genus} directory name
│	│
│	└─ parameters :
│		-st AA
│		-nt 16
│		-bb 16
│		-seed 1234
│		-m TEST (to find the best model for tge data)
│		-pre ${Genus}_Phylogeny
│
├── execPythonScript.sh :
│	│
│	│  For this script to work, we assume that the headers have a structure :
│	│  ${LocusTag}_000xx or ${LocusTag}_000xx|strain_name
│	│
│	│
│	├─ array :
│	│	- Genus : Genus names corresponding to the ${Genus} directory name
│	└── parameters :
│		- arg1 : script
│		- arg2 : input directory -> path to 04_Aligned_OG
│		- arg3 : output directory -> path to 05_Whole_Alignments (scripts creates the outdir)
│		- arg4 : pipeNames -> set to TRUE if headers contains a pipe.
│
├── OrthoFinder.sh :
│	│
│	├─ array :
│	│	- Genus : Genus names corresponding to the ${Genus} directory name
│	│
│	└─ parameters :
│		-f ${FilesLocation} : input files location, path to 03_GenesFaa
│		-t 16 : threads
│		-n ${Genus} : name for the output diresctory -> ${Genus}_Results
│
├── SAGE_FetchAndConcat2.py :
│	│  python script to prune the OG alignments.
│	│  Positions in the alignments with >50% gaps "-" are removed.
│	│  Writes all the pruned OGs
│	│  Writes the concatenated alignment
│	└── parameters : (in the execPythonScript)
│		- arg1 : script
│		- arg2 : input directory -> path to 04_Aligned_OG
│		- arg3 : output directory -> path to 05_Whole_Alignments (scripts creates the outdir)
│		- arg4 : pipeNames -> set to TRUE if headers contains a pipe.
│
└── pySED.py :
	│
	│  Sed python function do change a SINGLE string in a file OR
	│  to replace multiple patterns by others with a correspondance file.
	│
	│  Desinged for phylogenetic trees, column name of the old patter : 'Locus_Tag'
	│  New pattern is the concatenation of 'Specie' + 'Strain_name' columns
	│
	└─ parameters :
		- arg1 : True if multiple patterns and correspondance file
				 False if a single pattern to replace
		- arg2 : File path : where to change pattern
		- arg3 : output file name : path will be the same as the input file
		- arg4 : Table with patterns correspondance if TRUE
				 Replacement string if FALSE

"""

rule annotate:
    input


rule run_orthofinder:
    input:
        faa_files = expand("database/faa_files/{{phylotype}}/{genome}.faa", genome=get_g_list_by_phylo(DBs["microbiome"], "{phylotype}")),
        faa_dir = "database/faa_files/{phylotype}"
    output:
        ortho_file = "database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt"
    conda: "envs/core-cov-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    threads: 5
    log: "logs/{phylotype}_run_orthofinder.log"
    shell:
        """
        # so that nothing created in advance messes up the name
        rm -rf "database/faa_files/{wildcards.phylotype}/OrthoFinder"
        orthofinder -og -t {threads} -n dir -f {input.faa_dir}
        """

rule get_single_ortho:
    input:
        ortho_file = "database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt"
    output:
        ortho_single = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho.txt"
    params:
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{phylotype}_get_single_ortho.log"
    threads: 5
    script:
        "scripts/get_single_ortho.py"

rule extract_orthologs:
    input:
        ortho_single = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho.txt",
    output:
        ortho_seq_dir = directory("database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_ortho_sequences/")
    params:
        ffn_dir = "database/ffn_files",
        faa_dir = "database/faa_files/{phylotype}",
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{phylotype}_extract_orthologs.log"
    threads: 5
    script:
        "scripts/extract_orthologs.py"

rule calc_perc_id:
    input:
        ortho_seq_dir = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_ortho_sequences/"
    output:
        perc_id = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_perc_id.txt"
    params:
        scripts_dir = os.path.join(os.getcwd(), "scripts"),
        meta = os.path.join(os.getcwd(), "database/"+DBs["microbiome"]+"_metafile.txt"),
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-3:10:00"),
    resources:
        mem_mb = 8000
    conda: "envs/core-cov-env.yaml"
    log: "logs/{phylotype}_calc_perc_id.log"
    threads: 5
    shell:
        """
        cd {input.ortho_seq_dir}
        bash {params.scripts_dir}/aln_calc.sh {params.scripts_dir} {params.meta} ../{wildcards.phylotype}_perc_id.txt
        """
