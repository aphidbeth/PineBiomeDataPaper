# Snakemake file for PineBiome Metabarcoding

# SNAKEMAKE SETUP =========================
configfile: "config.yaml" 

print('Config details:')
print(config)

#============================================

# DEFINE OUTPUT

rule all:
    """
    Defines the final output files we want from this pipeline. In this case we want an ASV table per sample. 
    """
    input: 
        expand("{proj_dir}/manifests/{run_ID}_{amplicon}_manifest.csv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon="16S"),
        expand("{proj_dir}/manifests/{run_ID}_{amplicon}_manifest.csv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon="ITS"),
        expand("{proj_dir}/artifacts/raw_seqs_{run_ID}_{amplicon}.qza", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/artifacts/raw_seqs_{run_ID}_{amplicon}_trimmed.qza", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/sequence_quality_plots_{run_ID}_{amplicon}.qzv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/denoising_stats_{run_ID}_{amplicon}.qzv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/raw_seqs_{run_ID}_{amplicon}.qzv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/raw_seqs_{run_ID}_{amplicon}-trimmed.qzv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/denoising_stats_{run_ID}_{amplicon}.tsv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/artifacts/merged_runs/{featuretab}_{timepoint}_{amplicon}.qza", proj_dir=config["proj_dir"], amplicon=["ITS","16S"], featuretab = ["table", "rep-seqs"], timepoint=config["timepoint"]),
        expand("{proj_dir}/artifacts/merged_runs/rep-seqs_{timepoint}_{amplicon}-nochlor_nomito.qza", proj_dir=config["proj_dir"], amplicon=["ITS","16S"], timepoint=config["timepoint"]),
        expand("{proj_dir}/artifacts/merged_runs/merged_taxonomy_{timepoint}_{amplicon}.qza", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/raw_read_counts_{run_ID}_{amplicon}.tsv", proj_dir = config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/post_fastp_filtering_read_counts_{run_ID}_{amplicon}.tsv", proj_dir = config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/post_ITSxpress_read_counts_{run_ID}_{amplicon}.tsv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/freq_table_seq_{run_ID}_{amplicon}.tsv", proj_dir=config["proj_dir"], run_ID=config["run_ID"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/merged_runs/freq_table_{timepoint}_{amplicon}-nochlor_nomito.tsv", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/stats/merged_runs/freq_table_{timepoint}_{amplicon}-nochlor_nomito_nounassigned.tsv", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/merged_runs/barplots_{timepoint}_{amplicon}.qzv", proj_dir=config["proj_dir"],  timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/visualisations/merged_runs/seqs-and-taxonomy_{timepoint}_{amplicon}.qzv", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/outputs/seqs-and-taxonomy_{timepoint}_{amplicon}.tsv", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"]),
        expand("{proj_dir}/outputs/phylogeny/tree_{timepoint}_{amplicon}.qza", proj_dir=config["proj_dir"], timepoint=config["timepoint"], amplicon=["ITS","16S"])

#=============================================
# DATA PREPROCESSING STEPS:

# The following files must be run before this pipeline is executed: 
# 1) Repeat-demultiplexing.sh 
# 2) Rename_T2.sh
# 3) Rename_T3_run5.sh
# 4) Rename_T3_run6.sh      # These first four will create a set of folders T1_run1_new, T1_run2_new, T2_run3_new, T2_run4_new, T3_run5_new, T3_run6_new.
# 5) Pine-Metabar-QC-generic.sh  # Run once per folder, with the folder name given as a postional arguement 
# 6) Run-itsxpress-generic.sh   # Run once per folder, with the folder name given as a postional arguement 

# =============================================
# START QIIME PIPELINE


rule make_16S_manifests:
    """
    Creates a manifest file for importing reads into qiime2 by searching through a specified cleaned reads directory and creating list consisting of filename, filepath and read orientation.
    """
    params: 
        cleaned_reads_dir=config["cleaned_reads_16S"]
    output: 
        expand("{proj_dir}/manifests/{{run_ID}}_{amplicon}_manifest.csv", proj_dir=config["proj_dir"], amplicon="16S")
    shell: 
        """
        reads_in={params.cleaned_reads_dir}/{wildcards.run_ID}
        cd $reads_in

        # Define the bit you want to cut from the file names, leaving only the sample name
        cutMEr1="_cleaned.1.fastq.gz" 
        cutMEr2="_cleaned.2.fastq.gz" 

        for item in `ls *16S*.1.fastq.gz` ; do echo `basename $item $cutMEr1` ;done > $TMPDIR/{wildcards.run_ID}_16S_names
        for item in `ls *16S*.2.fastq.gz` ; do echo `basename $item $cutMEr2` ;done >>  $TMPDIR/{wildcards.run_ID}_16S_names

        for item in `ls *16S*.1.fastq.gz` ; do printf $reads_in/$item'\n'; done >  $TMPDIR/{wildcards.run_ID}_16S_reads
        for item in `ls *16S*.2.fastq.gz` ; do printf $reads_in/$item'\n'; done >>  $TMPDIR/{wildcards.run_ID}_16S_reads

        for item in `ls *16S*.1.fastq.gz`  ; do printf "forward\n" ; done >  $TMPDIR/{wildcards.run_ID}_16S_direction
        for item in `ls *16S*.2.fastq.gz`  ; do printf "reverse\n" ; done >>  $TMPDIR/{wildcards.run_ID}_16S_direction

        printf "sample-id,absolute-filepath,direction\n" > $TMPDIR/{wildcards.run_ID}_16S_tempfile.csv
        paste -d,  $TMPDIR/{wildcards.run_ID}_16S_names  $TMPDIR/{wildcards.run_ID}_16S_reads  $TMPDIR/{wildcards.run_ID}_16S_direction >>  $TMPDIR/{wildcards.run_ID}_16S_tempfile.csv

        mv $TMPDIR/{wildcards.run_ID}_16S_tempfile.csv {output}

        rm  $TMPDIR/{wildcards.run_ID}_16S* 

        """

rule make_ITS_manifests:
    """
    Creates a manifest file for importing reads into qiime2 by searching through a specified cleaned reads directory and creating list consisting of filename, filepath and read orientation.
    """
    params: 
        cleaned_reads_dir=config["cleaned_reads_ITS"],
        amplicon="ITS"
    output: 
        expand("{proj_dir}/manifests/{{run_ID}}_{amplicon}_manifest.csv", proj_dir=config["proj_dir"], amplicon="ITS")
    shell: 
        """
        reads_in={params.cleaned_reads_dir}/{wildcards.run_ID}
        cd $reads_in

        # Define the bit you want to cut from the file names, leaving only the sample name
        cutMEr1="_trimmed_reads.1.fastq.gz" 
        cutMEr2="_trimmed_reads.2.fastq.gz" 

        for item in `ls *ITS*.1.fastq.gz` ; do echo `basename $item $cutMEr1` ;done > $TMPDIR/{wildcards.run_ID}_ITS_names
        for item in `ls *ITS*.2.fastq.gz` ; do echo `basename $item $cutMEr2` ;done >>  $TMPDIR/{wildcards.run_ID}_ITS_names

        for item in `ls *ITS*.1.fastq.gz` ; do printf $reads_in/$item'\n'; done >  $TMPDIR/{wildcards.run_ID}_ITS_reads
        for item in `ls *ITS*.2.fastq.gz` ; do printf $reads_in/$item'\n'; done >>  $TMPDIR/{wildcards.run_ID}_ITS_reads

        for item in `ls *ITS*.1.fastq.gz`  ; do printf "forward\n" ; done >  $TMPDIR/{wildcards.run_ID}_ITS_direction
        for item in `ls *ITS*.2.fastq.gz`  ; do printf "reverse\n" ; done >>  $TMPDIR/{wildcards.run_ID}_ITS_direction

        printf "sample-id,absolute-filepath,direction\n" > $TMPDIR/{wildcards.run_ID}_ITS_tempfile.csv
        paste -d,  $TMPDIR/{wildcards.run_ID}_ITS_names  $TMPDIR/{wildcards.run_ID}_ITS_reads  $TMPDIR/{wildcards.run_ID}_ITS_direction >>  $TMPDIR/{wildcards.run_ID}_ITS_tempfile.csv

        mv $TMPDIR/{wildcards.run_ID}_ITS_tempfile.csv {output}

        rm  $TMPDIR/{wildcards.run_ID}_ITS* 

        """


rule import_reads: 
    """
    Uses the manifest to import reads into the qiime2 artifact
    """
    input: 
        expand("{proj_dir}/manifests/{{run_ID}}_{{amplicon}}_manifest.csv", proj_dir=config["proj_dir"])
    output:
        expand("{proj_dir}/artifacts/raw_seqs_{{run_ID}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    threads: 2
    resources:
        mem_mb= 1000
    shell: 
        """
        qiime tools import \
        --type SampleData[PairedEndSequencesWithQuality] \
        --input-path {input} \
        --input-format PairedEndFastqManifestPhred33 \
        --output-path {output}
        """

#=============================================

# QC AND CLEANING

rule cutadapt:
    """
    Remove any adapters from reads
    """
    input: 
        rules.import_reads.output
    output:
        expand("{proj_dir}/artifacts/raw_seqs_{{run_ID}}_{{amplicon}}_trimmed.qza", proj_dir=config["proj_dir"])
    log:
        expand("{proj_dir}/logs/cutadapt_logs/{{run_ID}}_{{amplicon}}.out", proj_dir=config["proj_dir"])
    params: 
        forwardseq=lambda wildcards: config["primers"][wildcards.amplicon + "_forward"],
        reverseseq=lambda wildcards: config["primers"][wildcards.amplicon + "_reverse"]
    conda:
        config["conda_env"]
    threads: 2
    resources:
        mem_mb= 1000
    shell:
        """
            echo "forward sequence being used is {params.forwardseq}"
            echo "reverse sequence being used {params.reverseseq}"

            qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input} \
            --p-front-f {params.forwardseq} \
            --p-front-r {params.reverseseq} \
            --o-trimmed-sequences {output} \
            --verbose > {log}
        """

rule check_seq_quality: 
    """
    Output visualation of read quality per run used decide parameters for the subsequent denoising step. 
    """
    input: 
        rules.cutadapt.output
    output: 
         expand("{proj_dir}/visualisations/sequence_quality_plots_{{run_ID}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime demux summarize \
        --i-data {input} \
        --o-visualization {output}
        """
    

rule dada2_denoise: 
    """
    This runs dada2denoise plugin. Dada2 trims reads, removes low quality reads, merges the forward and reverse reads then removes chimeric reads. 
    """
    input: 
        rules.cutadapt.output
    output: 
        repseqs = expand("{proj_dir}/artifacts/rep-seqs_{{run_ID}}_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        table = expand("{proj_dir}/artifacts/table_{{run_ID}}_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        stats = expand( "{proj_dir}/artifacts/denoising_stats_{{run_ID}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    params: 
        run_specific_params = lambda wildcards: config["denoise_params"][wildcards.amplicon][wildcards.run_ID]
    threads: 6
    resources:
        mem_mb= 15000
    conda:
        config["conda_env"]
    shell:
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input} \
        {params.run_specific_params} \
        --p-trunc-q 2 \
        --p-n-threads {threads} \
        --o-table {output.table} \
        --o-representative-sequences {output.repseqs} \
        --o-denoising-stats {output.stats}
        """

rule denoising_stats:
    input: 
        rules.dada2_denoise.output.stats
    output: 
        expand("{proj_dir}/visualisations/denoising_stats_{{run_ID}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"])
    conda: 
        config["conda_env"]
    shell: 
        """
        qiime metadata tabulate \
        --m-input-file {input} \
        --o-visualization {output}
        """

rule convert_denoising_stats_to_tsv:
    input:
        rules.denoising_stats.output
    output: 
         expand("{proj_dir}/stats/denoising_stats_{{run_ID}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    params: 
        folderdir = expand("{proj_dir}/qiime_tools_exports/denoising_stats_{{run_ID}}_{{amplicon}}", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime tools export \
        --input-path {input} \
        --output-path {params.folderdir}

        # rename stats file and remove parent folder
        mv {params.folderdir}/*.tsv {output}
        rm -R {params.folderdir} 
        """

#=============================================
# MOVE AND MERGE RUNS 

rule split_T1_repeats: 
    """
    Spilts the T3 tables based on metadata column into T3_run6 and T1 repeats.
    """
    input: 
        T3_run6_tab = expand("{proj_dir}/artifacts/table_T3_run6_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run6_metadat = expand("{proj_dir}/manifests/T3_run6_new_{{amplicon}}_metadata.tsv", proj_dir=config["proj_dir"])
    output: 
        T3_run6_noT1_tab = expand("{proj_dir}/artifacts/table_T3_run6_new_{{amplicon}}_noT1.qza", proj_dir=config["proj_dir"]),
        T1_repeats_tab = expand("{proj_dir}/artifacts/tabrepeats_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
    conda: 
        config["conda_env"]
    shell:
        """
        qiime feature-table filter-samples\
        --i-table {input.T3_run6_tab} \
        --m-metadata-file {input.T3_run6_metadat} \
        --p-where '[timepoint]="T3"' \
        --o-filtered-table {output.T3_run6_noT1_tab}

        qiime feature-table filter-samples \
        --i-table {input.T3_run6_tab} \
        --m-metadata-file {input.T3_run6_metadat} \
        --p-where '[timepoint]="T1"'\
        --o-filtered-table {output.T1_repeats_tab}

        """

rule split_T1_repeats_repseqs: 
    """
    Create new repseqs based on the tables
    """
    input: 
        T1_repeats_tab = expand("{proj_dir}/artifacts/tabrepeats_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run6_noT1_tab = expand("{proj_dir}/artifacts/table_T3_run6_new_{{amplicon}}_noT1.qza", proj_dir=config["proj_dir"]),
        T3_run6_rs = expand("{proj_dir}/artifacts/rep-seqs_T3_run6_new_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    output: 
        T3_run6_noT1_rs = expand("{proj_dir}/artifacts/rep-seqs_T3_run6_new_{{amplicon}}_noT1.qza", proj_dir=config["proj_dir"]),
        T1_repeats_rs = expand("{proj_dir}/artifacts/rep-seqsrepeats_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    conda: 
        config["conda_env"]
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data  {input.T3_run6_rs} \
        --i-table  {input.T3_run6_noT1_tab}  \
        --o-filtered-data {output.T3_run6_noT1_rs}

        qiime feature-table filter-seqs \
        --i-data  {input.T3_run6_rs} \
        --i-table  {input.T1_repeats_tab}  \
        --o-filtered-data {output.T1_repeats_rs}

        """

rule merge_tables:
    """
    Merges tables based on timepoint
    """
    input:
        T1_run1 = expand("{proj_dir}/artifacts/table_T1_run1_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T1_run2 = expand("{proj_dir}/artifacts/table_T1_run2_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2_run3 = expand("{proj_dir}/artifacts/table_T2_run3_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2_run4 = expand("{proj_dir}/artifacts/table_T2_run4_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run5 = expand("{proj_dir}/artifacts/table_T3_run5_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run6_noT1 = expand("{proj_dir}/artifacts/table_T3_run6_new_{{amplicon}}_noT1.qza", proj_dir=config["proj_dir"]),
        T1_repeats = expand("{proj_dir}/artifacts/tabrepeats_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    output: 
        T1 = expand("{proj_dir}/artifacts/merged_runs/table_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2 = expand("{proj_dir}/artifacts/merged_runs/table_T2_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3 = expand("{proj_dir}/artifacts/merged_runs/table_T3_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table merge \
        --i-tables {input.T1_run1} {input.T1_run2} {input.T1_repeats} \
        --o-merged-table {output.T1}

        qiime feature-table merge \
        --i-tables {input.T2_run3} {input.T2_run4} \
        --o-merged-table {output.T2}

        qiime feature-table merge \
        --i-tables {input.T3_run5} {input.T3_run6_noT1} \
        --o-merged-table {output.T3}
        """

rule merge_seqs:
    """
    Merges seqs based on timepoint
    """
    input:
        T1_run1 = expand("{proj_dir}/artifacts/rep-seqs_T1_run1_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T1_run2 = expand("{proj_dir}/artifacts/rep-seqs_T1_run2_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2_run3 = expand("{proj_dir}/artifacts/rep-seqs_T2_run3_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2_run4 = expand("{proj_dir}/artifacts/rep-seqs_T2_run4_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run5 = expand("{proj_dir}/artifacts/rep-seqs_T3_run5_new_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3_run6_noT1 = expand("{proj_dir}/artifacts/rep-seqs_T3_run6_new_{{amplicon}}_noT1.qza", proj_dir=config["proj_dir"]),
        T1_repeats = expand("{proj_dir}/artifacts/rep-seqsrepeats_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    output: 
        T1 = expand("{proj_dir}/artifacts/merged_runs/rep-seqs_T1_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T2 = expand("{proj_dir}/artifacts/merged_runs/rep-seqs_T2_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        T3 = expand("{proj_dir}/artifacts/merged_runs/rep-seqs_T3_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table merge-seqs \
        --i-data {input.T1_run1} {input.T1_run2} {input.T1_repeats} \
        --o-merged-data {output.T1}

        qiime feature-table merge-seqs \
        --i-data {input.T2_run3} {input.T2_run4} \
        --o-merged-data {output.T2}

        qiime feature-table merge-seqs \
        --i-data {input.T3_run5} {input.T3_run6_noT1} \
        --o-merged-data {output.T3}
        """


#=============================================

# TAXONOMIC ASSIGNMENT

rule classify_reads:
    input: 
        repseqs = expand("{proj_dir}/artifacts/merged_runs/rep-seqs_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    output: 
        expand("{proj_dir}/artifacts/merged_runs/merged_taxonomy_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    params:
        classifier = lambda wildcards: config["classifier_dir"][wildcards.amplicon]
    threads: 6
    resources:
        mem_mb= 60000
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-classifier classify-sklearn \
        --p-n-jobs {threads} \
        --i-classifier {params.classifier}/trained_classifier.qza \
        --i-reads {input.repseqs} \
        --o-classification {output}
        """

# Filter for mitochondria, chloroplast and unassigned sequences

rule filter_table_mt_chlor:
    input: 
        table = expand("{proj_dir}/artifacts/merged_runs/table_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        taxonomy = rules.classify_reads.output
    output: 
        expand("{proj_dir}/artifacts/merged_runs/table_{{timepoint}}_{{amplicon}}-nochlor_nomito.qza", proj_dir=config["proj_dir"])
    params:
        exclude = "mitochondria,chloroplast"
    conda:
        config["conda_env"]
    shell:
        """
        qiime taxa filter-table \
        --i-table {input.table} \
        --i-taxonomy {input.taxonomy} \
        --p-exclude {params.exclude} \
        --o-filtered-table {output}
        """


rule filter_repseqs_mt_chlor:
    input: 
        table = rules.filter_table_mt_chlor.output,
        repseqs = expand("{proj_dir}/artifacts/merged_runs/rep-seqs_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    output: 
        expand("{proj_dir}/artifacts/merged_runs/rep-seqs_{{timepoint}}_{{amplicon}}-nochlor_nomito.qza", proj_dir=config["proj_dir"])
    params:
        exclude = "mitochondria,chloroplast"
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data  {input.repseqs} \
        --i-table  {input.table}  \
        --o-filtered-data {output}
        """

rule filter_table_unassigned:
    input: 
        table = rules.filter_table_mt_chlor.output,
        taxonomy = rules.classify_reads.output
    output: 
        expand( "{proj_dir}/artifacts/merged_runs/table_{{timepoint}}_{{amplicon}}-nochlor_nomito_nounassigned.qza", proj_dir=config["proj_dir"])
    params:
        exclude = "Unassigned"
    conda:
        config["conda_env"]
    shell:
        """
        qiime taxa filter-table \
        --i-table {input.table} \
        --i-taxonomy {input.taxonomy} \
        --p-exclude {params.exclude} \
        --o-filtered-table {output}
        """


rule filter_repseqs_unassigned:
    input: 
        table = rules.filter_table_unassigned.output,
        repseqs = rules.filter_repseqs_mt_chlor.output
    output: 
        expand( "{proj_dir}/artifacts/merged_runs/rep-seqs_{{timepoint}}_{{amplicon}}-nochlor_nomito_nounassigned.qza", proj_dir=config["proj_dir"])
    params:
        exclude = "Unassigned"
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data  {input.repseqs} \
        --i-table  {input.table}  \
        --o-filtered-data {output}
        """
    

#========================================================
# GATHER READ FILTERING STATISTICS

# READ COUNTS

rule gather_raw_and_fastp_cleaned_read_counts:
    """
    Output some inital read statistics per sample
    """
    params: 
        raw_reads = config["raw_reads_dir"],
        fastp_out = config["cleaned_reads_16S"]
    output: 
        raw_read_stats = expand("{proj_dir}/stats/raw_read_counts_{{run_ID}}_{{amplicon}}.tsv", proj_dir = config["proj_dir"]),
        cleaned_read_stats = expand("{proj_dir}/stats/post_fastp_filtering_read_counts_{{run_ID}}_{{amplicon}}.tsv", proj_dir = config["proj_dir"])
    shell: 
        """
         cd {params.raw_reads}/{wildcards.run_ID}
        
        for item in `ls *{wildcards.amplicon}*_R1_*.fastq.gz` ; do 
            echo `basename $item _R1_001.fastq.gz` >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.names
            echo {wildcards.run_ID} >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.runID
            echo {wildcards.amplicon} >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.amplicon
            echo $(zcat $item |wc -l)/4|bc >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.read_counts
        done

        paste $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.names $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.runID $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.amplicon $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.read_counts > $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.raw_read_counts
        
        mv $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.raw_read_counts {output.raw_read_stats}

        rm $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.read_counts $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.runID $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.amplicon $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.names

        cd {params.fastp_out}/{wildcards.run_ID}
        
        for item in `ls *{wildcards.amplicon}*_cleaned.1.fastq.gz` ; do 
            echo `basename $item _cleaned.1.fastq.gz` >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.names
            echo $(zcat $item |wc -l)/4|bc >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.read_counts
            echo {wildcards.run_ID} >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.runID
            echo {wildcards.amplicon} >> $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.amplicon
        done

        paste $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.names $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.runID $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.amplicon $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.read_counts > $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.cleaned_read_counts
        
        mv $TMPDIR/{wildcards.run_ID}.{wildcards.amplicon}.cleaned_read_counts {output.cleaned_read_stats}

        """

rule post_ITSxpress_read_counts_vis: 
    """
    Output read counts of our intial qiime 2 object - only the ITS read counts should differ to the fast-p read counts as this is post ITSxpress filtering
    """
    input: 
        rules.import_reads.output
    output: 
         expand("{proj_dir}/visualisations/raw_seqs_{{run_ID}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime demux summarize \
            --i-data {input} \
            --o-visualization {output}
        """

rule gather_post_ITSxpress_read_counts: 
    """
    Output read counts of our intial qza object - only the ITS read counts should differ to the fast-p read counts as this is post ITSxpress filtering
    """
    input: 
        rules.post_ITSxpress_read_counts_vis.output
    params: 
        expand("{proj_dir}/outputs/qiime_tools_exports/post_ITSxpress_seq_summary_{{run_ID}}_{{amplicon}}", proj_dir=config["proj_dir"])
    output: 
        expand("{proj_dir}/stats/post_ITSxpress_read_counts_{{run_ID}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {params}

        mv {params}/per-sample-fastq-counts.tsv {output}

        rm -R {params}
        """

rule post_cutadapt_read_counts_vis: 
    """
    Output read counts of our intial qiime 2 object - only the ITS read counts should differ to the fast-p read counts as this is post ITSxpress filtering
    """
    input: 
        rules.cutadapt.output
    output: 
         expand("{proj_dir}/visualisations/raw_seqs_{{run_ID}}_{{amplicon}}-trimmed.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime demux summarize \
            --i-data {input} \
            --o-visualization {output}
        """

rule gather_post_cutadapt_read_counts: 
    """
    Output read counts of our intial qza object - only the ITS read counts should differ to the fast-p read counts as this is post ITSxpress filtering
    """
    input: 
        rules.post_cutadapt_read_counts_vis.output
    params: 
        expand("{proj_dir}/outputs/qiime_tools_exports/post_cutadapt_seq_summary_{{run_ID}}_{{amplicon}}", proj_dir=config["proj_dir"])
    output: 
        expand("{proj_dir}/stats/post_cutadapt_read_counts_{{run_ID}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {params}

        mv {params}/per-sample-fastq-counts.tsv {output}

        rm -R {params}
        """

# ASV AND READ COUNTS

rule tabulate_frequencies_raw:
    input: 
        raw_table = rules.dada2_denoise.output.table
    output: 
        raw_table_freq = expand( "{proj_dir}/artifacts/freq_table_seq_{{run_ID}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table tabulate-sample-frequencies \
        --i-table {input.raw_table} \
        --o-sample-frequencies {output.raw_table_freq}
        """

rule tabulate_frequencies:
    input: 
        mt_ch_table = rules.filter_table_mt_chlor.output,
        unassigned_table = rules.filter_table_unassigned.output
    output: 
        mt_ch_table_freq = expand( "{proj_dir}/artifacts/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito.qza", proj_dir=config["proj_dir"]),
        unassigned_table_freq = expand( "{proj_dir}/artifacts/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito_nounassigned.qza", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime feature-table tabulate-sample-frequencies \
        --i-table {input.mt_ch_table} \
        --o-sample-frequencies {output.mt_ch_table_freq}

        qiime feature-table tabulate-sample-frequencies \
        --i-table {input.unassigned_table} \
        --o-sample-frequencies {output.unassigned_table_freq}

        """

rule convert_freq_to_vis_raw:
    input: 
        rules.tabulate_frequencies_raw.output
    output: 
        expand( "{proj_dir}/visualisations/merged_runs/freq_table_seq_{{run_ID}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"]),
    conda:
        config["conda_env"]
    shell:
        """
        qiime metadata tabulate \
        --m-input-file {input} \
        --o-visualization {output}

        """

rule convert_freq_to_vis:
    input: 
        mt_ch_table_freq  = rules.tabulate_frequencies.output.mt_ch_table_freq,
        unassigned_table_freq = rules.tabulate_frequencies.output.unassigned_table_freq
    output: 
        mt_ch_table_freq_vis = expand( "{proj_dir}/visualisations/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito.qzv", proj_dir=config["proj_dir"]),
        unassigned_table_freq_vis = expand( "{proj_dir}/visualisations/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito_nounassigned.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """

        qiime metadata tabulate \
        --m-input-file {input.mt_ch_table_freq} \
        --o-visualization {output.mt_ch_table_freq_vis}

        qiime metadata tabulate \
        --m-input-file {input.unassigned_table_freq} \
        --o-visualization {output.unassigned_table_freq_vis}

        """
   
rule export_freq_vis_raw:
    input: 
        rules.convert_freq_to_vis_raw.output
    output: 
        expand("{proj_dir}/stats/freq_table_seq_{{run_ID}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    params:
        expand("{proj_dir}/outputs/qiime_tools_exports/raw_counts_seq_{{run_ID}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
            qiime tools export \
                --input-path {input} \
                --output-path {params}

            mv {params}/metadata.tsv {output}
            
            #rm -R {params}
        """


        
rule export_freq_vis:
    input: 
        mt_ch_table_freq_vis = rules.convert_freq_to_vis.output.mt_ch_table_freq_vis,
        unassigned_table_freq_vis = rules.convert_freq_to_vis.output.unassigned_table_freq_vis
    output: 
        mt_ch_table_tsv = expand( "{proj_dir}/stats/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito.tsv", proj_dir=config["proj_dir"]),
        unassigned_table_tsv = expand( "{proj_dir}/stats/merged_runs/freq_table_{{timepoint}}_{{amplicon}}-nochlor_nomito_nounassigned.tsv", proj_dir=config["proj_dir"])
    params:
        folder_dir2 = expand("{proj_dir}/outputs/qiime_tools_exports/merged_runs/post_filter_mt_ch_counts_{{timepoint}}_{{amplicon}}", proj_dir=config["proj_dir"]),
        folder_dir3 = expand("{proj_dir}/outputs/qiime_tools_exports/merged_runs/post_filter_unassigned_counts_{{timepoint}}_{{amplicon}}", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell:
        """
        qiime tools export \
            --input-path {input.mt_ch_table_freq_vis} \
            --output-path {params.folder_dir2}

        mv {params.folder_dir2}/metadata.tsv {output.mt_ch_table_tsv}

        qiime tools export \
            --input-path {input.unassigned_table_freq_vis} \
            --output-path {params.folder_dir3}

        mv {params.folder_dir3}/metadata.tsv {output.unassigned_table_tsv}

        #rm -R {params.folder_dir2}
        #rm -R {params.folder_dir3}
        """


#=================================================================================================================================
# OUTPUT BARCHARTS & TREES:


rule make_barplots:
    input: 
        table = rules.filter_table_mt_chlor.output,
        taxonomy = rules.classify_reads.output
    output:
        vis = expand("{proj_dir}/visualisations/merged_runs/barplots_{{timepoint}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell: 
        """
        qiime taxa barplot \
        --i-table {input.table} \
        --i-taxonomy {input.taxonomy} \
        --o-visualization {output.vis}

        """

        
# Build trees: 

rule run_iqtree: 
    """
    Builds a phylogeny for each of the timepoints
    """
    input:
        rules.filter_repseqs_unassigned.output
    output: 
        tree = expand("{proj_dir}/outputs/phylogeny/tree_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"]),
        rooted_tree = expand("{proj_dir}/outputs/phylogeny/rooted_tree_{{timepoint}}_{{amplicon}}.qza", proj_dir=config["proj_dir"])
    params: 
        output_dir = expand("{proj_dir}/outputs/phylogeny/tree_{{timepoint}}_{{amplicon}}", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell: 
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences {input} \
        --output-dir {params.output_dir}

        cp {params.output_dir}/tree.qza {output.tree}
        cp {params.output_dir}/rooted_tree.qza {output.rooted_tree}
        """

rule output_fullseq_list: 
    """
    Outputs repseqs against taxcnomy 
    """
    input:
        repseqs = rules.filter_repseqs_unassigned.output,
        taxonomy = rules.classify_reads.output
    output: 
        expand("{proj_dir}/visualisations/merged_runs/seqs-and-taxonomy_{{timepoint}}_{{amplicon}}.qzv", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell: 
        """
        qiime metadata tabulate \
        --m-input-file {input.repseqs} \
        --m-input-file {input.taxonomy} \
        --o-visualization {output}

        """

rule output_fullseq_list_tsv: 
    """
    Outputs repseqs against taxonomy 
    """
    input:
        rules.output_fullseq_list.output
    output: 
        expand("{proj_dir}/outputs/seqs-and-taxonomy_{{timepoint}}_{{amplicon}}.tsv", proj_dir=config["proj_dir"])
    params:
        folder_dir = expand("{proj_dir}/outputs/qiime_tools_exports/seqs_and_taxonomy_{{timepoint}}_{{amplicon}}", proj_dir=config["proj_dir"])
    conda:
        config["conda_env"]
    shell: 
        """
        qiime tools export \
            --input-path {input} \
            --output-path {params.folder_dir}

        mv {params.folder_dir}/metadata.tsv {output}

        """  