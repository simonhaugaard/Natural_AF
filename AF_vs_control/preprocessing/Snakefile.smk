# SNAKEMAKE FILE FOR PAIRED END READS AND STRANDED LIBRARY SETTINGS

# The implemented steps are:
# - Quality assessment with FastQC
# - Trimming and adapter removal with fastp
# - Alignment with STAR against Equus caballus genome 
# - Feature count of gene expression using reversely stranded library setting

# All rules are adapted from the wrappers documentation
# https://snakemake-wrappers.readthedocs.io/en/stable/

# Dataset: Define samples and their identifiers
dataset = ["Dorado_LA_post", "Dorado_LVFW", "Dorado_RA_post", "Dorado_RVFW", 
           "Im_A_Mets_Fan_LA_post", "Im_A_Mets_Fan_LVFW", "Im_a_Mets_Fan_RA_post", "Im_a_Mets_Fan_RVFW", 
           "Jytte_LA_post", "Jytte_LVFW", "Jytte_RA_post", "Jytte_RVFW", 
           "Kevin_Cook_LA_post", "Kevin_Cook_LVFW", "Kevin_Cook_RA_post", "Kevin_Cook_RVFW", 
           "MAP11_LA_post", "MAP11_LV_midwall", "MAP11_RApost", "MAP11_RV_endo", 
           "MAP3_LA_post", "MAP3_LV_midwall", "MAP3_RApost", "MAP3_RV_endo", 
           "MAP4_LA_post", "MAP4_LV_midwall", "MAP4_RApost", "MAP4_RV_endo", 
           "MAP6_LA_post", "MAP6_LV_midwall", "MAP6_RApost", "MAP6_RV_endo", 
           "MAP8_LA_post", "MAP8_LV_midwall", "MAP8_RApost", "MAP8_RV_endo", 
           "MAP9_LA_post", "MAP9_LV_midwall", "MAP9_RApost", "MAP9_RV_endo", 
           "San_Diego_LA_post", "San_Diego_LVFW", "San_Diego_RA_post", "San_Diego_RVFW", 
           "Styled_LA_post", "Styled_LVFW", "Styled_RA_post", "Styled_RVFW"]

# Rule to define final outputs without wildcards
rule all:
    input:
        "qc-raw-report/multiqc.html",
        "qc-clean-report/multiqc.html",
        "mapping-report/multiqc.html",
        "count-report/multiqc.html",
        "gene-expression-all-reverse-stranded-countReadPairs.tsv",
        "bam_collection.txt",
        "qualimap_collection.txt"

################################################################################
############################### QUALITY CONTROL ################################
################################################################################

# Assess quality of reads with FastQC
rule fastqc:
    input:
        "{folder}-data/{file}.fastq"
    output:
        html = "fastqc/{folder}-data/{file}.html",
        zip = "fastqc/{folder}-data/{file}_fastqc.zip"
    log:
        "logs/fastqc/{folder}/{file}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v3.3.3/bio/fastqc"

# MultiQC report for the raw data
rule raw_data:
    input:
        expand("fastqc/raw-data/{sample}_R1_fastqc.zip", sample=dataset),
        expand("fastqc/raw-data/{sample}_R2_fastqc.zip", sample=dataset)
    output:
        "qc-raw-report/multiqc.html",
        directory("qc-raw-report/multiqc_data")
    log:
        "logs/multiqc-raw.log"
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
############################ READ PRE-PROCESSING ###############################
################################################################################

# Pre-process reads with fastp
rule fastp_pe:
    input:
        sample = ["raw-data/{sample}_R1.fq.gz", "raw-data/{sample}_R2.fq.gz"]
    output:
        trimmed = ["clean-data/{sample}_1.fq.gz", "clean-data/{sample}_2.fq.gz"],
        html = "clean-data/{sample}.html",
        json  = "clean-data/{sample}-fastp.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra = " --average_qual=20 --qualified_quality_phred=20 --unqualified_percent_limit=10 --length_required=40 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    threads: 2
    wrapper:
        "v3.3.3/bio/fastp"

# MultiQC report for the cleaned data
rule clean_data:
    input:
        expand("fastqc/clean-data/{sample}_1_fastqc.zip", sample=dataset),
        expand("fastqc/clean-data/{sample}_2_fastqc.zip", sample=dataset),
        expand("clean-data/{sample}-fastp.json", sample=dataset)
    output:
        "qc-clean-report/multiqc.html",
        directory("qc-clean-report/multiqc_data")
    log:
        "logs/multiqc-clean.log"
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
############################## ALIGNMENT & MAPPING #############################
################################################################################

# Index the genome for STAR
rule star_index:
    input:
        fasta = "data/Equus_caballus.EquCab3.0.dna.toplevel.fa"
    output:
        directory("Equus_caballus.EquCab3.0.dna.toplevel_index")
    threads: 8
    log:
        "logs/star_index.log"
    wrapper:
        "v3.3.3/bio/star/index"

# Map cleaned reads against the indexed genome using STAR
rule star_pe_multi:
    input:
        fq1 = "clean-data/{sample}_1.fq.gz",
        fq2 = "clean-data/{sample}_2.fq.gz",
        idx = "Equus_caballus.EquCab3.0.dna.toplevel_index"
    output:
        aln = "star/{sample}.bam",
        log = "star/{sample}.Log.out",
        log_final = "star/{sample}.Log.final.out",
        sj  = "star/{sample}.sj.out.tab"
    log:
        "logs/star/{sample}.log"
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"

# Report on the mapping
rule map_data:
    input:
        expand("star/{sample}.Log.final.out", sample=dataset)
    output:
        "mapping-report/multiqc.html",
        directory("mapping-report/multiqc_data")
    log:
        "logs/multiqc-mapping.log"
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
############################## QUANTIFICATION ##################################
################################################################################

# Quantify gene expression using featureCounts with reverse-stranded setting
rule feature_counts_reverse_stranded:
    input:
        samples = "star/{sample}.bam",
        annotation = "data/Equus_caballus.EquCab3.0.111.gtf"
    output:
        multiext(
            "featureCount/countReadPairs/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        )
    threads: 2
    params:
        strand = 2,  # Set to 2 for reversely stranded
        extra = " --largestOverlap -p -B -P -t exon -F GFF --countReadPairs"
    log:
        "logs/reverse_stranded/{sample}.log"
    wrapper:
        "v3.3.3/bio/subread/featurecounts"

# Combine feature counts into a single expression matrix for all samples
rule feature_counts_combine_reverse_stranded:
    input:
        expand("featureCount/countReadPairs/{sample}.featureCounts", sample=dataset)
    output:
        "gene-expression-all-reverse-stranded-countReadPairs.tsv"
    shell:
        """
        tmp=$(mktemp -d)
        tail -n +2 {input[0]} | cut -f 1-6 > $tmp/00_annot
        for i in {input}; do
            bsn=$(basename $i .featureCounts)
            echo $bsn > $tmp/$bsn
            tail -n +3 $i | cut -f 7 >> $tmp/$bsn
        done
        paste $tmp/* > {output}
        rm -rf $tmp
        """

# MultiQC report on quantification
rule count_data:
    input:
        expand("featureCount/countReadPairs/{sample}.featureCounts.summary", sample=dataset)
    output:
        "count-report/multiqc.html",
        directory("count-report/multiqc_data")
    log:
        "logs/multiqc-count.log"
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
############################### FINAL REPORTS ##################################
################################################################################

# Rule for collecting BAM files
rule collect_bam:
    input:
        expand("star/{sample}.sorted.bam", sample=dataset)
    output:
        "bam_collection.txt"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))

# Rule for collecting Qualimap results
rule collect_qualimap:
    input:
        expand("qc-qualimap/{sample}", sample=dataset)
    output:
        "qualimap_collection.txt"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))
