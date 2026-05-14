#!/bin/bash

#===========================================================
# Description
#===========================================================
# Script to prepare databases for RNAseq host and non-host
# sequence separation using SortMeRNA and HISAT2.

# Closest related genome to Pinus sylvestis available is 
# Pinus tabuliformis so I use this as host genome

# Written by Beth Moore 2026 
# with assistance from github copilot

#===========================================================
# Prepare folder structure
#===========================================================
#scratch_dir= #redacted


# For SortMeRNA
sortmeRNA_dir="$scratch_dir/rRNA_databases_v4"
mkdir -p $sortmeRNA_dir

# For Genome download
ptabu_dir=$scratch_dir/genomes/ptabuliformis_genome
phix_dir=$scratch_dir/genomes/phix_genome
pchloro_dir=$scratch_dir/genomes/psylvestris_chloroplast
pslyv_mito_dir=$scratch_dir/genomes/psylvestris_mitochondria
mkdir -p $ptabu_dir
mkdir -p $phix_dir
mkdir -p $pchloro_dir
mkdir -p $pslyv_mito_dir

# For HISAT2 indexes
hisat2_index_dir=$scratch_dir/hisat2_indexes
mkdir -p $hisat2_index_dir

#================================================================================
# Download and index databases for SortMeRNA 
#================================================================================
# echo "Downloading SortMeRNA databases"

# # SortMeRNA databases
# wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz -O $sortmeRNA_dir/database.tar.gz
# tar -xvf $sortmeRNA_dir/database.tar.gz -C $sortmeRNA_dir
# rm $sortmeRNA_dir/database.tar.gz

# echo "SortMeRNA databases downloaded"

#================================================================================
# Download and index databases for HISAT2
#================================================================================

# echo "Downloading Pinus tabuliformis genome"
# cd $ptabu_dir

# conda activate base

# # Download Pinus tabuliformis genome from NCBI datasets
# cd $ptabu_dir
# datasets download genome accession GCA_031772625.1 --filename ptauliformis.zip
# unzip ptauliformis.zip -d $ptabu_dir/ptauliformis

# mv $ptabu_dir/ncbi_dataset/data/GCA_031772625.1/GCA_031772625.1_CPIR1_genomic.fna $ptabu_dir/ptabuliformis.fna
# rm ncbi_dataset -R
# rm md5sum.txt
# rm README.md

# conda deactivate

# echo "Pinus tabuliformis genome downloaded"

# Build HISAT2 index for host genome (Pinus tabuliformis)
#-----------------------------------------------------------------------

source activate TestPineRNA


# hisat2-build -p $SLURM_CPUS_PER_TASK \
#     $ptabu_dir/ptabuliformis.fna \
#     $hisat2_index_dir/ptabuliformis_index \
#     --large-index \
#     --seed 7


# echo "HISAT2 index for Pinus tabuliformis built"

# gzip $ptabu_dir/ptabuliformis.fna

# HISAT2 non-host genome index
# Note we concatenate here the host genome and additional mt and chl DNA sequences to avoid misclassification of organellar reads as non-host

# # Download PhiX genome
# #-----------------------
# echo "dowloading PhiX genome"

# cd $phix_dir
# curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000819615.1/download?include_annotation_type=GENOME_FASTA&filename=GCF_000819615.1.zip" -H "Accept: application/zip"
# unzip GCF_000819615.1.zip
# mv $phix_dir/ncbi_dataset/data/GCF_000819615.1/GCF_000819615.1_ViralProj14015_genomic.fna $phix_dir
# ls | grep -xv "GCF_000819615.1_ViralProj14015_genomic.fna" | xargs rm -R
# gzip GCF_000819615.1_ViralProj14015_genomic.fna
# mv GCF_000819615.1_ViralProj14015_genomic.fna.gz $phix_dir/phix_genome.gz

# echo "PhiX genome downloaded" 

# # Download Pinus sylvestris chloroplast genome
# #----------------------------------------------
# echo "Downloading Pinus chloroplast genome"

# wget https://www.ebi.ac.uk/ena/browser/api/fasta/KR476379.1?download=true  -O $pchloro_dir/KR476379.1.fna 
# gzip $pchloro_dir/KR476379.1.fna
# mv $pchloro_dir/KR476379.1.fna.gz $pchloro_dir/Pinus_sylvestris_chloroplast_genome.fasta.gz

# echo "Pinus chloroplast genome downloaded"


# # Download Pinus sylvestris mitochondrial genome
# #----------------------------------------------
# echo "Downloading Pinus mitochondrial genome"

# cd $pslyv_mito_dir
# wget https://www.ebi.ac.uk/ena/browser/api/fasta/KY302806.1?download=true
# mv 'KY302806.1?download=true' pinus_sylvestris_mitochondrial_genome.fasta
# gzip pinus_sylvestris_mitochondrial_genome.fasta

# echo "Pinus mitochondrial genome downloaded"

# Contatenate host genome and non target (organellar & phiX genomes)
#-----------------------------------------------------------------------

zcat $ptabu_dir/ptabuliformis.fna.gz \
    $pchloro_dir/Pinus_sylvestris_chloroplast_genome.fasta.gz \
    $phix_dir/phix_genome.gz \
    $pslyv_mito_dir/pinus_sylvestris_mitochondrial_genome.fasta.gz > $scratch_dir/genomes/pinus_host_non_target_combined.fa

# Build HISAT2 index for metatranscriptome analysis (non-host genome)
#-----------------------------------------------------------------------

hisat2-build -p $SLURM_CPUS_PER_TASK\
             $scratch_dir/genomes/pinus_host_non_target_combined.fa \
             $hisat2_index_dir/pinus_nonhost_index \
            --large-index \
            --seed 7


echo "HISAT2 index for non-host genome built"

gzip $scratch_dir/genomes/pinus_host_non_target_combined.fa
echo "Host and non target genomes concatenated"
