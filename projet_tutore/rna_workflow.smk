########################################################################################################################
# Workflow purpose :                                                                                                   #
#     Treats Chlamydomonas reinhardtii timeseries RNAseq datasets from their raw state through quality check and       #
#     alignment on an indexed reference genome to a single gene counts table with all the librairies.                  #
# Author : Domitille Jarrige                                                                                           #
# Date   : 2021-06-12                                                                                                  #
#----------------------------------------------------------------------------------------------------------------------#
#                     CHANGE PARAMETERS AND DIRECTORIES IN THE CONFIG FILE BEFORE RUNNING                              #
#----------------------------------------------------------------------------------------------------------------------#
# A test dataset is available at https://github.com/JarrigeD/dubii2021/tree/master/projet_tutore                       #
# to run workflow:       snakemake --drmaa  --snakefile {SNAKEFILE} --configfile {CONFIG_FILE}                         #        
########################################################################################################################

SAMPLES, = glob_wildcards(config["data_directory"] + "{sample}.fastq.gz")
SAMPLES.sort()
FILES = ["Aligned.sortedByCoord.out.bam", "Log.final.out", "Log.out", "Log.progress.out", "SJ.out.tab"]
INDEX = ["chrLength.txt", "geneInfo.tab", "sjdbInfo.txt", "chrNameLength.txt", "Genome", "sjdbList.fromGTF.out.tab", "chrName.txt", "genomeParameters.txt", "sjdbList.out.tab", "chrStart.txt", "Log.out", "transcriptInfo.tab", "exonGeTrInfo.tab", "SA", "exonInfo.tab", "SAindex"]


rule all:
    input: 
        expand(config["results_directory"] + "{sample}_on_Crev5.6/{file}", sample = SAMPLES, file = FILES),
        expand(config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv", sample = SAMPLES),
        config["results_directory"] + "total_counts.csv",
        expand(config["quality_control_directory"] + "{sample}_fastqc.html", sample = SAMPLES),
        expand(config["quality_control_directory"] + "{sample}_fastqc.zip", sample = SAMPLES),
        config["quality_control_directory"] + "multiqc_report.html"

rule fastqc:
    input: config["data_directory"] + "{sample}.fastq.gz"
    output:
        html = config["quality_control_directory"] + "{sample}_fastqc.html",
        zip = config["quality_control_directory"] + "{sample}_fastqc.zip"
    threads: config["parameters"]["threads"]
    params:
        fastqc_version = config["parameters"]["fastqc_version"],
        outdir = config["quality_control_directory"]
    shell: """
        mkdir -p {params.outdir}
        module load {params.fastqc_version}
        fastqc {input} --threads={threads} -o {params.outdir}
    """

rule multiqc:
    input: config["quality_control_directory"][:-1]
    output: config["quality_control_directory"] + "multiqc_report.html"
    params:
        multiqc_version = config["parameters"]["multiqc_version"],
        outdir = config["quality_control_directory"][:-1]
    shell: """
        module load {params.multiqc_version}
        multiqc {input} -o {params.outdir} -f
    """


rule STAR_genome_indexation:
    input:
        nuclear = config["nuclear_genome"],
        organelles = config["organelles_genome"],
        annotation = config["annotation"]
    output: expand(config["genome_directory"] + "star/{index}", index = INDEX)
    threads: config["parameters"]["threads"]
    params:
        star_version = config["parameters"]["star_version"]
    shell : """
        module load {params.star_version}
        STAR --runMode genomeGenerate --runThreadN 8 --genomeDir {output} --genomeFastaFiles {input.nuclear} {input.organelles} --sjdbGTFfile {input.annotation} --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 --sjdbGTFtagExonParentGene ID --sjdbGTFtagExonParentGeneName Name 
    """



rule STAR_alignment:
    input: 
        genome = config["genome_directory"] + "star", 
        fastq = config["data_directory"] + "{sample}.fastq.gz"
    output: [config["results_directory"] + "{sample}_on_Crev5.6/Aligned.sortedByCoord.out.bam", config["results_directory"] + "{sample}_on_Crev5.6/Log.final.out", config["results_directory"] + "{sample}_on_Crev5.6/Log.out", config["results_directory"] + "{sample}_on_Crev5.6/Log.progress.out", config["results_directory"] + "{sample}_on_Crev5.6/SJ.out.tab"]
    params:
        parentdir = config["results_directory"], 
        outdir = config["results_directory"] + "{sample}_on_Crev5.6/",
        star_version = config["parameters"]["star_version"]
    threads: config["parameters"]["threads"]
    shell: """
        module load {params.star_version}
        mkdir -p {params.parentdir}
        STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.fastq} --outFileNamePrefix {params.outdir} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 30 --alignIntronMin 20 --alignIntronMax 10000
    """
    
rule feature_counts:
    input: 
        annotation = config["annotation"],
        bam_file = config["results_directory"] + "{sample}_on_Crev5.6/Aligned.sortedByCoord.out.bam"
    output: config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv"
    threads: config["parameters"]["threads"]
    params: 
        subread_version = config["parameters"]["subread_version"]
    shell: """
    module load {params.subread_version}
    featureCounts -a {input.annotation} -o {output} {input.bam_file} -g "ID" -T {threads} -M --fraction -O --fraction -t "gene" --extraAttributes "Name"
    """
    
 
rule generate_temp_tables:
    input:         
        count_table = config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv",
        sample_ref = config["sample_ref"]
    output: config["results_directory"] + "tmp_{sample}.table"
    params:
        parentdir = config["results_directory"]
    shell: """
        mkdir -p {params.parentdir}
        sample_id=$(basename -s _counts.csv {input.count_table})
        sample_condition=$(grep ${{sample_id}} {input.sample_ref} | cut -f 5)
        echo ${{sample_condition}} > {output}
        cut -f 8 {input.count_table} | tail -n +3 >> {output}
    """


rule generate_final_table:
    input: 
        tables = expand(config["results_directory"] + "tmp_{sample}.table", sample=SAMPLES),
        gene_names = config["gene_names"]
    output: config["results_directory"] + "total_counts.csv"
    params:
        parentdir = config["results_directory"]
    shell: """
        #cut -f 1  > gene_names
        paste {input.tables} > tmp_out
        paste {input.gene_names} tmp_out > {output}
        rm {params.parentdir}/tmp_SRR*
        rm tmp_out

    """




    