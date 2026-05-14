#!/bin/bash

#SBATCH # Header redacted

source activate Nextseq_soil

# Define timepoint => this needs changing for each pool
if [ ! $@ ]; then
    printf 'Positional argument needed to define the run folder containing raw reads i.e. "T1_run2".'
fi

timepoint=$1

# Define folders
# Filepaths redacted
homedir= # raw reads directory
tempdir= # scratch directory


mkdir -p $tempdir/PineBiome-Metabarcoding/HPC/outputs/fastp/cleaned_reads/$timepoint
reads_out=$tempdir/outputs/fastp/cleaned_reads/$timepoint
echo "read directory: $reads_out"
	
mkdir -p $tempdir/outputs/fastp/reports/$timepoint
reports_out=$tempdir/outputs/fastp/reports/$timepoint
echo "reports directory: $reports_out"

# Run fastp
cd $homedir

# PINE METABARCODING QC
# Written by Beth Moore 08/24

# Run fastp

for file1 in $homedir/*_R1_001.fastq.gz
do
        file2=${file1%%_R1_001.fastq.gz}"_R2_001.fastq.gz" 
		nosuffix=${file1%_R*} 
		sample_name=${nosuffix#*"$timepoint"\/}
	echo "Running $sample_name"
	fastp -i $file1 \
		  -I $file2 \
		  --trim_poly_g \
		  --length_required 90 \
		  -o $reads_out/"$sample_name"_cleaned.1.fastq.gz \
		  -O $reads_out/"$sample_name"_cleaned.2.fastq.gz\
		  -h $reports_out/"$sample_name".fastp.html \
		  -j $reports_out/"$sample_name".fastp.json \
          --thread 4
	>> $tempdir/logs/fastp_$timepoint.out
echo "Completed $sample_name"
done

# Integrate htmls using MultiQC

multiqc $tempdir/outputs/fastp/reports/$timepoint/*16S* --outdir $tempdir/outputs/fastp/reports/$timepoint/ --title cleaned_"$timepoint"_16S_reads_v1
multiqc $tempdir/outputs/fastp/reports/$timepoint/*ITS* --outdir $tempdir/outputs/fastp/reports/$timepoint/ --title cleaned_"$timepoint"_ITS_reads_v1
