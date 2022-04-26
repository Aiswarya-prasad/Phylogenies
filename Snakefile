import os

configfile: "config.yaml"

Project_identifier = config["Project_identifier"]
Project_path = config["Project_path"]
GenomeInfo = config["Genome_info_path"]
Groups = config["Phylo_groups"]

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
            genome = line.split("\t")[0]
            genome = genome.strip()
            # exclude the unpublished genome!
            if genome == "Ga0418777":
                # unpublished genome
                continue
            genomeList.append(genome)
    return(genomeList)

def get_g_list_by_group(path, group_name):
    """
    Returns list of genomes in given group based on the Genome info file (path)
    """
    g_list = []
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r") as info_fh:
        for line in info_fh:
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            if genome == "Ga0418777":
                # unpublished genome
                continue
            group = line.split("\t")[12]
            if group == group_name:
                g_list.append(genome)


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
        prokka_faa = expand("01_prokka/{genome}/{genome}.faa", genome=get_genomes(GenomeInfo)),
        out_dir = expand("03_aligned_orthogroups/{group}/", group=Groups),

rule download_genome:
    input:
        info = GenomeInfo,
    output:
        genome = "00_genomes/{genome}.fa",
        genome_gz = temp("00_genomes/{genome}.fa.gz"),
    threads: 1
    params:
        ftp_summary = "https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt",
        assembly_summary_genbank = "assembly_summary_genbank.txt",
    log: "logs/{genome}_download.log"
    shell:
        """
        if [ -f {output.genome} ]; then
            # if snakemake is working, this part should never run!
            echo "{output.genome} exists." 2>&1 | tee -a {log}
        else
            echo "cannot find {output.genome}. Going to try to download!" 2>&1 | tee -a {log}
            wget -N "{params.ftp_summary}" 2>&1 | tee -a {log}
            strain=$(cat {input.info} | grep {wildcards.genome} | cut -f4)
            echo "The strain is ${{strain}}" 2>&1 | tee -a {log}
            genbank_path=$(cat {params.assembly_summary_genbank} | grep "strain="${{strain}} | cut -f15,20 | sort -k1 -r | head -1 | cut -f2)
            echo "${{genbank_path}}_here" 2>&1 | tee -a {log}
            if [ -z ${{genbank_path}} ]; then
                echo "No genbank path found for genome {wildcards.genome}. Exiting." 2>&1 | tee -a {log}
            else
                echo "Downloading from ${{genbank_path}}/$(basename "${{genbank_path}}")_genomic.fna.gz" 2>&1 | tee -a {log}
                wget -O {output.genome}.gz "${{genbank_path}}/$(basename "${{genbank_path}}")_genomic.fna.gz" 2>&1 | tee -a {log}
                ls {output.genome}.gz
                echo "Downloaded! unzipping {output.genome}.gz" 2>&1 | tee -a {log}
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
        faa = temp("01_prokka/{genome}/{genome}.faa"),
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
        prokka --compliant --force \
            --outdir {params.outdir} \
            --locustag {wildcards.genome} \
            --prefix {wildcards.genome} \
            --evalue 0.001 \
            {input.genome} 2>&1 | tee -a {log}
            #? --proteins proteins.faa \
        """

rule prepare_faa:
    input:
        faa_files = expand("01_prokka/{genome}/{genome}.faa", genome=get_genomes(GenomeInfo)),
        info = GenomeInfo
    output:
        faa_files = expand("02_orthofinder/{{group}}/{genome}.faa", genome=get_g_list_by_group(GenomeInfo, "{group}")),
        ortho_dir = "02_orthofinder/{group}/"
    log: "logs/{group}_prepare_faa"
    shell:
        """
        for file in {input.faa_files}
        do
            cp ${{file}} {output.ortho_dir}  2>&1 | tee -a {log}
        done
        """

rule run_orthofinder:
    input:
        faa_dir = "02_orthofinder/{group}/",
        faa_files = expand("02_orthofinder/{{group}}/{genome}.faa", genome=get_g_list_by_group(GenomeInfo, "{group}")),
    output:
        ortho_file = "02_orthofinder/{group}/OrthoFinder/Results_{group}/Orthogroups/Orthogroups.txt",
    conda: "envs/phylogenies-env.yaml"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
    threads: 4
    log: "logs/{group}_run_orthofinder.log"
    shell:
        """
        orthofinder -og -t {threads} -n {wildcards.group} -f {input.faa_dir}  2>&1 | tee -a {log}
        """

rule align_orthologs:
    input:
        ortho_dir = "02_orthofinder/{group}/OrthoFinder/Results_{group}/Single_Copy_Orthologue_Sequences/"
    output:
        out_dir = directory("03_aligned_orthogroups/{group}/"),
    conda: "envs/phylogenies-env.yaml"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
    threads: 4
    log: "logs/{group}_run_orthofinder.log"
    shell:
        """
        for OG in $(ls {input.ortho_dir})
        do
            echo "starting alignment for ${{OG}}" 2>&1 | tee -a {log}
            mafft --amino --inputorder --localpair --maxiterate 1000 {input.ortho_dir}/${{OG}} > {output.out_dir}/${{OG/.fa/_aligned.fa}} 2>&1 | tee -a {log}
        done
        """
