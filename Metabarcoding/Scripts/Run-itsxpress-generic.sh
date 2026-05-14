#!/bin/bash

#SBATCH #header redacted

# Define timepoint => this needs changing for each pool
if [ ! $@ ]; then
    printf 'Positional argument needed to define the run folder containing raw reads i.e. "T1_run2".'
fi

timepoint=$1

source activate itsxpressenv

# Define version number of run - to be changed with any major revisions
runno="24_09_24"
amplicon="ITS" # we can look at inputting these as an array for the final run of the pipeline but for now just submit a script per batch

# File paths redacted 
tempdir= # scratch dir
reads_in= # fastp cleaned reads

#==========================================================================================================

mkdir -p $tempdir/outputs/itsxpress_out/$timepoint
mkdir -p $tempdir/logs/itsxpress_out # output log files here

# CODE FOR RUNNING ITS EXPRESS SEPERATELY=================
cd $reads_in

cutMEr1="_cleaned.1.fastq.gz" 
cutMEr2="_cleaned.2.fastq.gz" 

for item in `ls *"$amplicon"*.1.fastq.gz` ; do
echo `basename $item $cutMEr1` ;
sample=$(basename $item $cutMEr1) ;
itsxpress --fastq "$sample"_cleaned.1.fastq.gz --fastq2 "$sample"_cleaned.2.fastq.gz --region ITS2 \
--taxa Fungi --log $tempdir/logs/itsxpress_out/"$timepoint"_itsexpress_"$sample"_logfile.txt \
--threads 10 \
--outfile $tempdir/outputs/itsxpress_out/"$timepoint"/"$sample"_trimmed_reads.1.fastq.gz \
--outfile2 $tempdir/outputs/itsxpress_out/"$timepoint"/"$sample"_trimmed_reads.2.fastq.gz ;
echo "$sample done" ;
done 
