"""
Run using

snakemake -p -r --use-conda --conda-prefix /home/aiswarya/anaconda3/phylogenies-environment --cores 8

or

snakemake -p -r --use-conda --conda-prefix /home/aiswarya/anaconda3/phylogenies-environment --cores 8  --rerun-incomplete --keep-going
"""

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
            # unpublished genome
            if genome == "Ga0418777":
                continue
            genomeList.append(genome)
    return(genomeList)

def get_g_list_by_group(path, group_name):
    """
    Returns list of genomes in given group based on the Genome info file (path)
    """
    g_list = []
    print(f"getting genome list for {group_name}")
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r") as info_fh:
        for line in info_fh:
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            # unpublished genome
            if genome == "Ga0418777":
                continue
            group = line.split("\t")[12]
            # only include groups of interest!
            if group not in Groups:
                continue
            if group == group_name:
                g_list.append(genome)
        return(g_list)

def get_g_dict_for_groups(path):
    """
    Returns dictionary of genomes and groups with each value being a list of
    genomes corresponding to a given group
    """
    g_list_dict = {}
    for group in Groups:
        g_list_dict[group] = []
    if os.path.isfile(path):
        pass
    else:
        print(f"Could not find file at {path}")
    with open(path, "r") as info_fh:
        for line in info_fh:
            if line.startswith("ID"):
                continue
            genome = line.split("\t")[0]
            # unpublished genome
            if genome == "Ga0418777":
                continue
            group = line.split("\t")[12]
            # only include groups of interest!
            if group not in Groups:
                continue
            g_list_dict[group].append(genome)
    return(g_list_dict)


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

rule all:
    input:
        alignment = expand("04_pruned_and_concat_alignments/{group}/CoreGeneAlignment.fasta", group=Groups),

rule download_genome:
    output:
        genome = "00_genomes/{genome}.fa",
        genome_gz = temp("00_genomes/{genome}.fa.gz"),
    threads: 1
    params:
        info = GenomeInfo,
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
            strain=$(cat {params.info} | grep {wildcards.genome} | cut -f4)
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
        faa = "01_prokka/{genome}/{genome}.faa",
    params:
        outdir = "01_prokka/{genome}/"
    #     mailto="aiswarya.prasad@unil.ch",
    #     mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
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
        """

rule prepare_faa:
    input:
        faa_files = expand("01_prokka/{genome}/{genome}.faa", genome=get_genomes(GenomeInfo)),
        info = GenomeInfo,
    output:
        # faa_files = expand("02_orthofinder/{{group}}/{genome}.faa", genome=get_g_dict_for_groups(GenomeInfo)["{group}"]),
        faa_files = [expand("02_orthofinder/"+group+"/{genome}.faa", genome=get_g_dict_for_groups(GenomeInfo)[group]) for group in Groups]
    log: "logs/prepare_faa.log"
    shell:
        """
        for file in {output.faa_files}
        do
            genome_name=$(basename ${{file}})
            if [ ! -f ${{file}} ]; then
                cp 01_prokka/${{genome_name/.faa/}}/${{genome_name}} ${{file}}
            fi
        done
        """

rule run_orthofinder:
    input:
        faa_files = lambda wildcards: expand("02_orthofinder/{{group}}/{genome}.faa", genome=get_g_dict_for_groups(GenomeInfo)[wildcards.group]),
    output:
        orthogroups_dir = directory("02_orthofinder/{group}/OrthoFinder/Results_{group}/Single_Copy_Orthologue_Sequences/"),
        orthofinder_dir = directory("02_orthofinder/{group}/OrthoFinder/Results_{group}/"),
    conda: "envs/phylogenies-env.yaml"
    # params:
    #     mailto="aiswarya.prasad@unil.ch",
    #     mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
    #     account="pengel_spirit",
    #     runtime_s=convertToSec("0-2:10:00"),
    #  resources:
    #      mem_mb = 8000
    threads: 4
    log: "logs/{group}_run_orthofinder.log"
    shell:
        """
        rm -rf {output.orthofinder_dir}*
        orthofinder -og -t {threads} -n {wildcards.group} -f $(dirname {input.faa_files[0]})  2>&1 | tee -a {log}
        """

rule align_orthologs:
    input:
        orthogroups_dir = "02_orthofinder/{group}/OrthoFinder/Results_{group}/Single_Copy_Orthologue_Sequences/"
    output:
        out_dir = directory("03_aligned_orthogroups/{group}/"),
        done = touch("03_aligned_orthogroups/{group}/mafft.done"),
    conda: "envs/phylogenies-env.yaml"
    # params:
    #     mailto="aiswarya.prasad@unil.ch",
    #     mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
    #     account="pengel_spirit",
    #     runtime_s=convertToSec("0-2:10:00"),
    # resources:
    #     mem_mb = 8000
    threads: 4
    log: "logs/{group}_align_orthologs.log"
    shell:
        """
        mkdir -p {output.out_dir}
        for OG in $(ls {input.orthogroups_dir})
        do
            echo "starting alignment for ${{OG}}" 2>&1 | tee -a {log}
            mafft --amino --inputorder --localpair --maxiterate 1000 {input.orthogroups_dir}/${{OG}} > {output.out_dir}/${{OG/.fa/_aligned.fa}}
        done
        """

# rule orthologs_summary:
# summarise how may orthogroups across groups to make sure we have enough
# we want ~ 500? or so

rule prune_and_concat:
    input:
        aligned_dir = "03_aligned_orthogroups/{group}/",
        done = "03_aligned_orthogroups/{group}/mafft.done",
    output:
        pruned_dir = directory("04_pruned_and_concat_alignments/{group}/"),
        pruned_cat = "04_pruned_and_concat_alignments/{group}/CoreGeneAlignment.fasta",
    params:
        pipe_names = False
    #     mailto="aiswarya.prasad@unil.ch",
    #     mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
    #     account="pengel_spirit",
    #     runtime_s=convertToSec("0-2:10:00"),
    # resources:
    #     mem_mb = 8000
    threads: 2
    log: "logs/{group}_prune_and_concat.log"
    script:
        "scripts/prune_and_concat_alns.py"

rule make_tree:
    input:
        pruned_cat = "04_pruned_and_concat_alignments/{group}/CoreGeneAlignment.fasta"
    output:
        pdf = "05_IQTree/{group}/{group}_Phylogeny.iqtree",
        iqlog = "05_IQTree/{group}/{group}_Phylogeny.log"
    params:
        outdir = "05_IQTree/{group}/"
    #     mailto="aiswarya.prasad@unil.ch",
    #     mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
    #     account="pengel_spirit",
    #     runtime_s=convertToSec("0-2:10:00"),
    # resources:
    #     mem_mb = 8000
    threads: 16
    log: "logs/{group}_make_tree.log"
    conda: "envs/phylogenies-env.yaml"
    shell:
        """
        mkdir -p 05_IQTree
        mkdir -p {params.outdir}
        iqtree -s {input.pruned_cat} \
                -st AA -nt {threads} \
                -bb 1000 -seed 12345 -m TEST \
                -pre {params.outdir}{wildcards.group}_Phylogeny
        """
