#!/bin/bash

#SBATCH # header redacted

# array params:
while getopts ":t:" opt; do
case ${opt} in
t )
filename=${OPTARG}
;;
* ) echo " Option -t is required [-t, <samplelist>]"
;;
esac
done

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $filename)

IFS='' read -r -a array <<< "$LINE"

# Notes: 
# Whilst 80G is not needed for most of the pipeline, the tabuliformis index is massive and needs this space to load. 

# Set our reference database path and file path as variables 
# filepaths redacted
hisat2_host_index= # path to host genome index
hisat2_nonhost_index= # path to non-host genome index
READ_DIR= # path to processed data directory
LOG_DIR= # path to logs directory

#Load conda enviroment
source activate TestPineRNA

echo "Aligning reads from sample ${array[0]} to ptabu index" 

echo "Aligning reads to target genomes and retaining both mathcing and non-matching reads"

        hisat2 \
        --seed 7 \
        --threads $SLURM_CPUS_PER_TASK \
        -x $hisat2_host_index \
        -1 $READ_DIR/${array[0]}_no_rRNA_R1.fastq.gz \
        -2 $READ_DIR/${array[0]}_no_rRNA_R2.fastq.gz \
        -S $TMPDIR/${array[0]}.sam \
        --score-min "L,0,-0.2"  \
        --no-unal \
        --al-conc-gz $TMPDIR/${array[0]}.hostreads.%.fq.gz \
        --un-conc-gz $TMPDIR/${array[0]}.nonhostreads.%.fq.gz \
        --summary-file $LOG_DIR/${array[0]}_hisat2_host.log
        
        echo "Finished aligning paired reads for ${array[0]}"
        echo "Reads in TMPDIR:"
        ls $TMPDIR

        echo "Copying across paired host reads and renaming"
        mv $TMPDIR/${array[0]}.hostreads.1.fq.gz $READ_DIR/${array[0]}_host_cleaned_R1.fastq.gz
        mv $TMPDIR/${array[0]}.hostreads.2.fq.gz $READ_DIR/${array[0]}_host_cleaned_R2.fastq.gz

        rm $TMPDIR/${array[0]}.sam

# Repeat for the Non host reads:

echo "Aligning reads from sample ${array[0]} to ptabu index + other non target genomes and retaining the reads that don't match to the non host index (i.e. putative microbial reads)" 

    hisat2 \
        --seed 7 \
        --threads $SLURM_CPUS_PER_TASK \
        -x $hisat2_nonhost_index \
        -1 $TMPDIR/${array[0]}.nonhostreads.1.fq.gz \
        -2 $TMPDIR/${array[0]}.nonhostreads.2.fq.gz \
        -S $TMPDIR/${array[0]}.sam \
        --score-min "L,0,-0.2"  \
        --no-unal \
        --al-conc-gz $TMPDIR/${array[0]}.nontargetreads.%.fq.gz \
        --un-conc-gz $TMPDIR/${array[0]}.targetreads.%.fq.gz \
        --summary-file $LOG_DIR/${array[0]}_hisat2_nonhost.log
        
        
        mv $TMPDIR/${array[0]}.targetreads.1.fq.gz $READ_DIR/${array[0]}_nonhost_cleaned_R1.fastq.gz
        mv $TMPDIR/${array[0]}.targetreads.2.fq.gz $READ_DIR/${array[0]}_nonhost_cleaned_R2.fastq.gz
      
        rm $TMPDIR/${array[0]}.sam
        rm $TMPDIR/${array[0]}.nontargetreads.1.fq.gz
        rm $TMPDIR/${array[0]}.nontargetreads.2.fq.gz

        echo "Finished extracting putative microbial reads for sample ${array[0]}"
        