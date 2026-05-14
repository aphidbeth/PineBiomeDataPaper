# Snakemake file for RNA file preparation prior to nextflow pipelines

# SNAKEMAKE SETUP ===========================

configfile: "RNA_preprocessing_config.yaml"

print('Config details:')
print(config)

scratch_dir = config["scratch_dir"]
proj_dir = config["proj_dir"]
SAMPLES = config["samples"]

#============================================

# DEFINE OUTPUT

rule all:
    """
    Defines the final output files we want from this pipeline. In this case we want a cleaned file for both host (psylvestris) and non-host (metagenomic) RNA. 
    """
    input: 
        expand("{scratch_dir}/data/processed/{sample}_host_cleaned_R1.fastq.gz", scratch_dir=config["scratch_dir"], sample=SAMPLES),
        expand("{scratch_dir}/data/processed/{sample}_host_cleaned_R2.fastq.gz", scratch_dir=config["scratch_dir"], sample=SAMPLES),
        expand("{scratch_dir}/data/processed/{sample}_nonhost_cleaned_R1.fastq.gz", scratch_dir=config["scratch_dir"], sample=SAMPLES),
        expand("{scratch_dir}/data/processed/{sample}_nonhost_cleaned_R2.fastq.gz", scratch_dir=config["scratch_dir"], sample=SAMPLES)

# STEPS ======================================


# 1) Quality control: removal of low quality, short and poly tail reads

rule fastp_clean:
    """
    Some basic quality control using fastp to trim adapters and low quality reads.
    """
    input:
        fwd = "{scratch_dir}/data/raw/{sample}_R1.fq.gz",
        rev = "{scratch_dir}/data/raw/{sample}_R2.fq.gz"
    output:
        fwd = "{scratch_dir}/data/processed/{sample}_cleaned_R1.fastq.gz",
        rev = "{scratch_dir}/data/processed/{sample}_cleaned_R2.fastq.gz",
        html = "{scratch_dir}/logs/fastp/{sample}_fastp_report.html",
        json = "{scratch_dir}/logs/fastp/{sample}_fastp_report.json"
    params:
        qual_thresh = config["fastp_qual_thresh"],
        length_thresh = config["fastp_length_thresh"]
    threads:
        config["fastp_threads"]
    log:
        "{scratch_dir}/logs/fastp/{sample}_fastp.log"
    resources:
        mem_mb = 4000
    conda: 
        config["conda_envs"]["fastp"]
    shell:
        """
        set -euo pipefail

        fastp -i {input.fwd} -I {input.rev} \
            -o {output.fwd} -O {output.rev} \
              --thread {threads} \
              --qualified_quality_phred {params.qual_thresh} \
              --length_required {params.length_thresh} \
              --trim_poly_x \
              --html {output.html} --json {output.json} \
              > {log} 2>&1
        """

# 2) Removal of rRNA reads

rule sortmerna_filter:
    """
    Removal of rRNA reads using SortMeRNA. Note the database must be prepared beforehand - see the Prepare_databases.sh script.
    """
    input:
        fwd = rules.fastp_clean.output.fwd,
        rev = rules.fastp_clean.output.rev
    output:
        fwd = "{scratch_dir}/data/processed/{sample}_no_rRNA_R1.fastq.gz",
        rev = "{scratch_dir}/data/processed/{sample}_no_rRNA_R2.fastq.gz"
    params:
        dbs = config["sortmerna_dbs"]
    threads:
        config["sortmerna_threads"]
    log:
        "{scratch_dir}/logs/sortmerna/{sample}_sortmerna.log"
    resources:
        mem_mb = 20000
    conda: 
        config["conda_envs"]["sortmerna"]
    shell:
        """
        set -euo pipefail

        sortmerna \
        -ref {params.dbs} \
        -reads {input.fwd} \
        -reads {input.rev} \
        -fastx \
        -other \
        -paired_out \
        -out2 \
        --zip-out \
        --blast '1 cigar qcov' \
        -threads {threads} \
        -workdir $TMPDIR/{wildcards.sample}_sortmeRNA

        # rename & move output files
        mv $TMPDIR/{wildcards.sample}_sortmeRNA/out/other_fwd.fq.gz {output.fwd}  
        mv $TMPDIR/{wildcards.sample}_sortmeRNA/out/other_rev.fq.gz {output.rev}
    
        """

# 3) Host and non-host read separation
# NOTE we ended up running this same code as an array to speed up processing but params remained the same

rule hisat2_get_host_reads:
    """
    Mapping reads to host genome using HISAT2 and separating host and non-host reads. Note the database must be prepared beforehand - see the prepare_databases.sh script.
    """
    input:
        fwd = rules.sortmerna_filter.output.fwd,
        rev = rules.sortmerna_filter.output.rev
    output:
        host_fwd = "{scratch_dir}/data/processed/{sample}_host_cleaned_R1.fastq.gz",
        host_rev = "{scratch_dir}/data/processed/{sample}_host_cleaned_R2.fastq.gz",
    params:
        hisat2_index = config["hisat2_host_index"]
    threads:
        config["hisat2_threads"]
    log:
        "{scratch_dir}/logs/hisat2/{sample}_hisat2_host.log"
    resources:
        mem_mb = 60000
    conda: 
        config["conda_envs"]["hisat2"]
    shell:
        """
        set -euo pipefail

        hisat2 \
        --seed 7 \
        --threads {threads} \
        -x {params.hisat2_index} \
        -1 {input.fwd} \
        -2 {input.rev} \
        --score-min "L,0,-0.2"  \
        -S $TMPDIR/{wildcards.sample}_host.sam \
        --no-unal \
        --al-conc-gz $TMPDIR/{wildcards.sample}.targetreads.%.fq.gz \
        --summary-file {log}

        mv $TMPDIR/{wildcards.sample}.targetreads.1.fq.gz {output.host_fwd}
        mv $TMPDIR/{wildcards.sample}.targetreads.2.fq.gz {output.host_rev}
        """

rule hisat2_get_nonhost_reads:
    """
    Mapping reads to a index made up of host, phiX, host mitochondrial and host chloroplast genomes to leave non-host reads
    that have potential for assembly into microbial transcriptomes.
    """
    input:
        fwd = rules.sortmerna_filter.output.fwd,
        rev = rules.sortmerna_filter.output.rev
    output:
        nonhost_fwd = "{scratch_dir}/data/processed/{sample}_nonhost_cleaned_R1.fastq.gz",
        nonhost_rev = "{scratch_dir}/data/processed/{sample}_nonhost_cleaned_R2.fastq.gz",
    params:
        hisat2_index = config["hisat2_nonhost_index"]
    threads:
        config["hisat2_threads"]
    log:
        "{scratch_dir}/logs/hisat2/{sample}_hisat2_nonhost.log"
    resources:
        mem_mb = 60000
    conda: 
        config["conda_envs"]["hisat2"]
    shell:
        """
        set -euo pipefail

        hisat2 \
        --seed 7 \
        --threads {threads} \
        -x {params.hisat2_index} \
        -1 {input.fwd} \
        -2 {input.rev} \
        --score-min "L,0,-0.2"  \
        -S $TMPDIR/{wildcards.sample}_nonhost.sam \
        --no-unal \
        --al-conc-gz $TMPDIR/{wildcards.sample}.alignedreads.%.fq.gz \
        --un-conc-gz $TMPDIR/{wildcards.sample}.targetreads.%.fq.gz \
        --summary-file {log}

        mv $TMPDIR/{wildcards.sample}.targetreads.1.fq.gz {output.nonhost_fwd}
        mv $TMPDIR/{wildcards.sample}.targetreads.2.fq.gz {output.nonhost_rev}
        """
