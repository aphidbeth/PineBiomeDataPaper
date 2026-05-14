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
SCRATCHDIR= # redacted

#Load conda enviroment
source activate TestPineRNA

echo "Aligning reads from sample ${array[0]}" 

echo "Aligning reads to non target genomes and retaining non-matching reads"

        hisat2 \
        --seed 7 \
        --threads $SLURM_CPUS_PER_TASK \
        -x $SCRATCHDIR/hisat2_out/indexes/hisat_non_target_index\
        -1 $SCRATCHDIR/sortmeRNA_out/${array[0]}_sortmeRNA_nonaligned_fwd.fq.gz \
        -2 $SCRATCHDIR/sortmeRNA_out/${array[0]}_sortmeRNA_nonaligned_rev.fq.gz \
        -S $TMPDIR/${array[0]}.sam \
        --score-min "L,0,-0.2"  \
        --no-unal \
        --al-conc-gz $TMPDIR/${array[0]}.nontargetreads.%.fq.gz \
        --un-conc-gz $TMPDIR/${array[0]}.targetreads.%.fq.gz \
        --summary-file $SCRATCHDIR/slurm_logs/Hisat2/${array[0]}.hisat_align.log
        
        echo "Finished aligning paired reads for ${array[0]}"
        echo "Reads in TMPDIR:"
        ls $TMPDIR

        echo "Copying across paired reads and renaming"
        cp $TMPDIR/${array[0]}.targetreads.1.fq.gz $SCRATCHDIR/hisat2_out/target_reads/${array[0]}.targetreads.1.fq.gz
        cp $TMPDIR/${array[0]}.targetreads.2.fq.gz $SCRATCHDIR/hisat2_out/target_reads/${array[0]}.targetreads.2.fq.gz
        
        cp $TMPDIR/${array[0]}.nontargetreads.1.fq.gz $SCRATCHDIR/hisat2_out/non_target_reads/${array[0]}.nontargetreads.1.fq.gz
        cp $TMPDIR/${array[0]}.nontargetreads.2.fq.gz $SCRATCHDIR/hisat2_out/non_target_reads/${array[0]}.nontargetreads.2.fq.gz
      
