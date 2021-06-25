#####################################################################################################################
# Workflow purpose :                                                                                                #
#     Treats Chlamydomonas reinhardtii timeseries RNAseq datasets from their raw state through quality check and    #
#     alignment on an indexed reference genome to a single gene counts table with all the librairies.               #
# Author : Domitille Jarrige                                                                                        #
# Date   : 2021-06-25                                                                                               #
# ------------------------------------------------------------------------------------------------------------------#
#                     CHANGE PARAMETERS AND DIRECTORIES IN THE CONFIG FILE BEFORE RUNNING                           #
# ------------------------------------------------------------------------------------------------------------------#
# A test dataset is available at https://github.com/JarrigeD/dubii2021/tree/master/projet_tutore                    #
# to run workflow:   snakemake --use-conda --cores={YOURCHOICE} --snakefile {SNAKEFILE} --configfile {CONFIG_FILE}  #
#                      Do not forget to build your project architecture before running!                             #
#####################################################################################################################

(SAMPLES,) = glob_wildcards(config["data_directory"] + "{sample}.fastq.gz")
SAMPLES.sort()
FILES = [
    "Aligned.sortedByCoord.out.bam",
    "Log.final.out",
    "Log.out",
    "Log.progress.out",
    "SJ.out.tab",
]
INDEX = [
    "chrLength.txt",
    "geneInfo.tab",
    "sjdbInfo.txt",
    "chrNameLength.txt",
    "Genome",
    "sjdbList.fromGTF.out.tab",
    "chrName.txt",
    "genomeParameters.txt",
    "sjdbList.out.tab",
    "chrStart.txt",
    "Log.out",
    "transcriptInfo.tab",
    "exonGeTrInfo.tab",
    "SA",
    "exonInfo.tab",
    "SAindex",
]


rule all:
    input:
        expand(
            config["results_directory"] + "{sample}_on_Crev5.6/{file}",
            sample=SAMPLES,
            file=FILES,
        ),
        expand(
            config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv",
            sample=SAMPLES,
        ),
        config["results_directory"] + "total_counts.tsv",
        config["results_directory"] + "gene_name.tab",
        expand(
            config["quality_control_directory"] + "{sample}_fastqc.html", sample=SAMPLES
        ),
        expand(
            config["quality_control_directory"] + "{sample}_fastqc.zip", sample=SAMPLES
        ),
        config["quality_control_directory"] + "multiqc_report.html",
        config["rhythm_directory"] + "rhythmic_genes_cosinor.csv",
        config["rhythm_directory"] + "qvalues_cs.png",


rule fastqc:
    input:
        config["data_directory"] + "{sample}.fastq.gz",
    output:
        html=config["quality_control_directory"] + "{sample}_fastqc.html",
        zip=config["quality_control_directory"] + "{sample}_fastqc.zip",
    threads: config["parameters"]["threads"]
    params:
        outdir=config["quality_control_directory"]
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input} --threads={threads} -o {params.outdir}
        """


rule multiqc:
    input:
        config["quality_control_directory"],
    output:
        config["quality_control_directory"] + "multiqc_report.html",
    params:
        outdir=config["quality_control_directory"][:-1],
    conda:
        "envs/envBash.yml"
    shell:
        """
        multiqc {input} -o {params.outdir} -f
        """


rule STAR_genome_indexation:
    input:
        nuclear=config["nuclear_genome"],
        organelles=config["organelles_genome"],
        annotation=config["annotation"],
    output:
        expand(config["genome_directory"] + "star/{index}", index=INDEX),
    threads: config["parameters"]["threads"]
    conda:
        "envs/envBash.yml"
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} --genomeFastaFiles {input.nuclear} {input.organelles} --sjdbGTFfile {input.annotation} --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 --sjdbGTFtagExonParentGene ID --sjdbGTFtagExonParentGeneName Name 
        """


rule STAR_alignment:
    input:
        genome=config["genome_directory"] + "star",
        fastq=config["data_directory"] + "{sample}.fastq.gz",
    output:
        [
            config["results_directory"]
            + "{sample}_on_Crev5.6/Aligned.sortedByCoord.out.bam",
            config["results_directory"] + "{sample}_on_Crev5.6/Log.final.out",
            config["results_directory"] + "{sample}_on_Crev5.6/Log.out",
            config["results_directory"] + "{sample}_on_Crev5.6/Log.progress.out",
            config["results_directory"] + "{sample}_on_Crev5.6/SJ.out.tab",
        ],
    params:
        parentdir=config["results_directory"],
        outdir=config["results_directory"] + "{sample}_on_Crev5.6/",
    threads: config["parameters"]["threads"]
    conda:
        "envs/envBash.yml"
    shell:
        """
        mkdir -p {params.parentdir}
        STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.fastq} --outFileNamePrefix {params.outdir} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 30 --alignIntronMin 20 --alignIntronMax 10000
        """


rule feature_counts:
    input:
        annotation=config["annotation"],
        bam_file=(
            config["results_directory"]
            + "{sample}_on_Crev5.6/Aligned.sortedByCoord.out.bam"
        ),
    output:
        config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv",
    threads: config["parameters"]["threads"]
    conda:
        "envs/envBash.yml"
    shell:
        """
        featureCounts -a {input.annotation} -o {output} {input.bam_file} -g "ID" -T {threads} -M --fraction -O --fraction -t "gene" --extraAttributes "Name"
        """


rule generate_temp_tables:
    input:
        count_table=(
            config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv"
        ),
        sample_ref=config["sample_ref"],
    output:
        config["results_directory"] + "tmp_{sample}.table",
    params:
        parentdir=config["results_directory"],
    shell:
        """
        mkdir -p {params.parentdir}
        sample_id=$(basename -s _counts.csv {input.count_table})
        sample_condition=$(grep ${{sample_id}} {input.sample_ref} | cut -f 5)
        echo ${{sample_condition}} > {output}
        cut -f 8 {input.count_table} | tail -n +3 >> {output}
        """


rule generate_ref_gene_file:
    input:
        expand(
            config["results_directory"] + "{sample}_on_Crev5.6/{sample}_counts.csv",
            sample=SAMPLES,
        ),
    output:
        config["results_directory"] + "gene_name.tab",
    shell:
        """
        tail -n +2 {input[1]} | cut -f 1 > {output}
        """


rule generate_final_table:
    input:
        tables=expand(
            config["results_directory"] + "tmp_{sample}.table", sample=SAMPLES
        ),
        gene_names=config["results_directory"] + "gene_name.tab",
    output:
        config["results_directory"] + "total_counts.tsv",
    params:
        parentdir=config["results_directory"],
        rhythm_dir=config["rhythm_directory"],
    shell:
        """
        paste {input.tables} > tmp_out
        paste {input.gene_names} tmp_out > {output}
        rm {params.parentdir}/tmp_SRR*
        rm tmp_out
        mkdir -p {params.rhythm_dir}
        """
        

rule find_rhythmic_genes:
    input: 
        meta = config["sample_ref"],
        counts = config["results_directory"] + "total_counts.tsv"
    output:
        config["rhythm_directory"] + "rhythmic_genes_cosinor.csv",
        config["rhythm_directory"] + "qvalues_cs.png"
    params:
        dir = config["rhythm_directory"]
    threads: config["parameters"]["threads"]
    conda:
        "envs/envR.yml"
    script: 
        "scripts/find_rhythmic_gene_script.R"
    