import os

configfile: "config.yaml"

ProjectIdentifier = config["ProjectIdentifier"]
ProjectPath = config["ProjectPath"]
GenomeInfo = config["Genome_info_path"]

"""
define functions

- read genome dictonary as needed
"""

def get_genomes(path):
    """
    reads genome info file (must be csv!)
    and return list of genomes
    """
    genomeList = []
    with open(path, "r") as info_fh:
        for line in info_fh:
            if line.startswith("ID"):
                continue
            genome = line.split(",")[0]
            genome = genome.strip()
            # exclude the unpublished genome!
            if genome == "Ga0418777":
                continue
            genomeList.append(genome)
    return(genomeList)

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

rule all:
    input:
        # genomes = expand("00_genomes/{genome}.fa", genome=get_genomes(GenomeInfo)),
        prokka_faa = expand("01_prokka/{genome}/{genome}.faa", genome=get_genomes(GenomeInfo)),

rule download_genome:
    # input:
    #     info = "GenomeInfo.csv",
    output:
        genome = "00_genomes/{genome}.fa",
        genome_gz = temp("00_genomes/{genome}.fa.gz"),
    threads: 1
    params:
        ftp_summary = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt",
        assembly_summary_genbank = "assembly_summary_genbank.txt",
        info = "GenomeInfo.csv",
    log: "logs/{genome}_download.log"
    shell:
        """
        # temporarily download the assembly summary to get paths
        if [ -f {output.genome} ]; then
            # if snakemake is working, this should never happen!
            echo "{output.genome} exists." 2>&1 | tee -a {log}
        else
            echo "cannot find {output.genome}. Going to try to download!" 2>&1 | tee -a {log}
            wget -N "{params.ftp_summary}" 2>&1 | tee -a {log}
            strain=$(cat {params.info} | grep {wildcards.genome} | cut -f4 -d",")
            echo "The strain is ${{strain}}" 2>&1 | tee -a {log}
            # all these steps are to avoid having multiple urls if there are
            # multiple genomes for the strain
            # we sort by date column (format "yyyy/mm/dd") in reversed order
            # (so most recent first) and take the top one
            genbank_path=$(cat {params.assembly_summary_genbank} | grep "strain="${{strain}} | cut -f15,20 | sort -k1 -r | head -1 | cut -f2)
            echo "${{genbank_path}}_here" 2>&1 | tee -a {log}
            if [ -z ${{genbank_path}} ]; then
                echo "No genbank path found for genome {wildcards.genome}. Exiting." 2>&1 | tee -a {log}
            else
                echo "Downloading from ${{genbank_path}}/$(basename "${{genbank_path}}")_genomic.fna.gz" 2>&1 | tee -a {log}
                wget -O {output.genome}.gz "${{genbank_path}}/$(basename "${{genbank_path}}")_genomic.fna.gz" 2>&1 | tee -a {log}
                ls {output.genome}.gz
                echo "Downloaded! unzipping {output.genome}.gz" 2>&1 | tee -a {log}
                # gzip -d {output.genome}.gz 2>&1 | tee -a {log}
                gunzip < {output.genome}.gz > {output.genome} 2>&1 | tee -a {log}
                echo "Unzipped! Editing header to containg only locus tag" 2>&1 | tee -a {log}
                echo ">{wildcards.genome}" > {output.genome}.temp
                cat {output.genome} | grep -v ">" >> {output.genome}.temp
                mv {output.genome}.temp {output.genome}
            fi
        fi
        """

rule annotate:
    input:
        genome = "00_genomes/{genome}.fa"
    output:
        faa = "01_prokka/{genome}/{genome}.faa",
    params:
        outdir = "01_prokka/{genome}/"
    #     mailto="aiswarya.prasad@unil.ch",
    #     account="pengel_spirit",
    #     runtime_s=convertToSec("0-2:10:00"),
    # resources:
    #     mem_mb = 8000
    threads: 4
    log: "logs/{genome}_annotate.log"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        start=${{SECONDS}}
        # use force because snakemake creates output dirs and this confuses it
        prokka --compliant --force \
            --outdir {params.outdir} \
            --locustag {wildcards.genome} \
            --prefix {wildcards.genome} \
            --evalue 0.001 \
            {input.genome} 2>&1 | tee -a {log}
            #? --proteins proteins.faa \
        duration=$(( SECONDS - start ))
        echo -e "The script ran for "${{duration}} "seconds" 2>&1 | tee -a {log}
        """


# rule run_orthofinder:
#     input:
#         faa_files = expand("database/faa_files/{{phylotype}}/{genome}.faa", genome=get_g_list_by_phylo(DBs["microbiome"], "{phylotype}")),
#         faa_dir = "database/faa_files/{phylotype}"
#     output:
#         ortho_file = "database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt"
#     conda: "envs/core-cov-env.yaml"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     threads: 5
#     log: "logs/{phylotype}_run_orthofinder.log"
#     shell:
#         """
#         # so that nothing created in advance messes up the name
#         rm -rf "database/faa_files/{wildcards.phylotype}/OrthoFinder"
#         orthofinder -og -t {threads} -n dir -f {input.faa_dir}
#         """
#
# rule get_single_ortho:
#     input:
#         ortho_file = "database/faa_files/{phylotype}/OrthoFinder/Results_dir/Orthogroups/Orthogroups.txt"
#     output:
#         ortho_single = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho.txt"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     log: "logs/{phylotype}_get_single_ortho.log"
#     threads: 5
#     script:
#         "scripts/get_single_ortho.py"
#
# rule extract_orthologs:
#     input:
#         ortho_single = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_single_ortho.txt",
#     output:
#         ortho_seq_dir = directory("database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_ortho_sequences/")
#     params:
#         ffn_dir = "database/ffn_files",
#         faa_dir = "database/faa_files/{phylotype}",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     log: "logs/{phylotype}_extract_orthologs.log"
#     threads: 5
#     script:
#         "scripts/extract_orthologs.py"
#
# rule calc_perc_id:
#     input:
#         ortho_seq_dir = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_ortho_sequences/"
#     output:
#         perc_id = "database/OrthoFiles_"+DBs["microbiome"]+"/{phylotype}_perc_id.txt"
#     params:
#         scripts_dir = os.path.join(os.getcwd(), "scripts"),
#         meta = os.path.join(os.getcwd(), "database/"+DBs["microbiome"]+"_metafile.txt"),
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-3:10:00"),
#     resources:
#         mem_mb = 8000
#     conda: "envs/core-cov-env.yaml"
#     log: "logs/{phylotype}_calc_perc_id.log"
#     threads: 5
#     shell:
#         """
#         cd {input.ortho_seq_dir}
#         bash {params.scripts_dir}/aln_calc.sh {params.scripts_dir} {params.meta} ../{wildcards.phylotype}_perc_id.txt
#         """
