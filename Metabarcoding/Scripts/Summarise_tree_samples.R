# SUMMARISE THE TREE SAMPLES IN THE PINEBIOME DATA PAPER

# Author: Beth Moore
# Updated: 15/01/2026

# Input -> Information on read counts of the successfully sequenced trees
#       -> List of the planned tree samples

# Output -> Summary table of number of trees

#==============================================================================================================

# SET UP 

library(tidyverse)
library(ggplot2)

rm(list=ls())

#===================================================================================================================
# READ IN DATA

# Planned list of control samples: 
planned_samples <- read.csv("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Sample_information/planned_tree_samples.csv")

# Extract the population and family information:

# Read in the trees which were sequenced: 
seq_samples <- read_tsv("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/Metabarcoding_Reads_per_sample_summary_table_all_steps.tsv")

#====================================================================================================================
# CURATE THE DATASET

# Merge the seq information
all_samples_m <- left_join(planned_samples, seq_samples, by = "Sample_ID") %>% select(!Sampling_timepoint) # drop the duplicate column

# Add notes on the drop out stage
all_samples_m$Notes <- rep("Reads retained in sample throughout processing", times = nrow(all_samples_m))
all_samples_m$Notes[is.na(all_samples_m$Final_reads_per_sample)|all_samples_m$Final_reads_per_sample== 0 ] <- "Dropped out during removal of outlying samples"
all_samples_m$Notes[is.na(all_samples_m$no_unassigned_reads)|all_samples_m$no_unassigned_reads== 0 ] <- "Dropped out during removal of non-fungal and unassigned reads"
all_samples_m$Notes[is.na(all_samples_m$Non_syn_comm)|all_samples_m$Non_syn_comm == 0 ] <- "Dropped out during removal of synthetic community (only applicable to non-synthetic community samples)"
all_samples_m$Notes[is.na(all_samples_m$Non_mito_chlor)|all_samples_m$Non_mito_chlor == 0 ] <- "Dropped out during removal of mitochondrial and chloroplast reads"
all_samples_m$Notes[is.na(all_samples_m$Non_chimeric)|all_samples_m$Non_chimeric == 0 ] <- "Dropped out during removal of chimeras"
all_samples_m$Notes[is.na(all_samples_m$Merged)|all_samples_m$Merged == 0 ] <- "Dropped out during merging"
all_samples_m$Notes[is.na(all_samples_m$Denoised)|all_samples_m$Denoised == 0 ] <- "Dropped out during denoising"
all_samples_m$Notes[is.na(all_samples_m$Filtered)|all_samples_m$Filtered == 0 ] <- "Dropped out during quality filtering"
all_samples_m$Notes[is.na(all_samples_m$Raw_reads)|all_samples_m$Raw_reads == 0 ] <- "Dropped out due to failed library preparation or sequencing"

# Count the catagories
table(all_samples_m$Notes)

# List all dropouts:
all_samples_m$Sample_ID[!(all_samples_m$Notes == "Reads retained in sample throughout processing")]

sample_dropouts <- all_samples_m[!(all_samples_m$Notes == "Reads retained in sample throughout processing"),]
sample_successful <- all_samples_m[all_samples_m$Notes == "Reads retained in sample throughout processing",]


#=====================================================================================================================
# OUTPUTS

# Create summary outputs from the paper

# Table of the samples which we have sequence data for 

sequenced_samples <- all_samples_m %>%
  filter(Notes == "Reads retained in sample throughout processing") %>%
  count(Timepoint) %>% 
  pivot_wider(names_from = Timepoint, values_from = n)


write_tsv(all_samples_m, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Sample_information/All_samples_information.tsv")

write_tsv(sequenced_samples, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Sample_information/Sequenced_samples_summary_table.tsv")

write_tsv(sample_dropouts, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Sample_information/Sample_dropouts.tsv")
write_tsv(sample_successful, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Sample_information/Sample_successful.tsv")
