#!/bin/bash

#SBATCH # header redacted

# Filepaths redacted
proj_dir= # raw data files
read_dir= # scratch output directory


# # Copying across files from the project directory to scratch dir with flattened structure
# # Also merges reads over multiple lanes into single R1 and R2 files per sample

# mkdir -p $read_dir

# for sample_dir in "$proj_dir"/*; do
#     sample=$(basename "$sample_dir")  # e.g., T3_7156

#     # Find all R1 and R2 files for this sample
#     r1_files=("$sample_dir"/"$sample"*"_1.fq.gz")
#     r2_files=("$sample_dir"/"$sample"*"_2.fq.gz")

#     # Determine output names (first 7 characters)
#     out_r1="$read_dir/${sample:0:7}_R1.fq.gz"
#     out_r2="$read_dir/${sample:0:7}_R2.fq.gz"

#     # Concatenate if more than one, otherwise just copy
#     if [ "${#r1_files[@]}" -gt 1 ]; then
#         echo "Concatenating $r1_files for $sample -> $out_r1"
#         cat "${r1_files[@]}" > "$out_r1"
#     else
#         echo "Copying single R1 file for $sample -> $out_r1"
#         cp "${r1_files[0]}" "$out_r1"
#     fi

#     if [ "${#r2_files[@]}" -gt 1 ]; then
#         echo "Concatenating $r2_files for $sample -> $out_r2"
#         cat "${r2_files[@]}" > "$out_r2"
#     else
#         echo "Copying single R2 file for $sample -> $out_r2"
#         cp "${r2_files[0]}" "$out_r2"
#     fi
# done

# Run snakemake pipeline
cd #snakemake pipeline directory redacted

source activate snakemake

snakemake --cores 20 --resources mem_mb=80000 --use-conda --snakefile "RNA_preprocessing.smk" --latency-wait 60 --rerun-incomplete
