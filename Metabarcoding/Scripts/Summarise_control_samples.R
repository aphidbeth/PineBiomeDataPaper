# SUMMARISE THE CONTROL SAMPLES IN THE PINEBIOME DATA PAPER

# Author: Beth Moore
# Updated: 15/01/2026

# Input -> Information on read counts of the successfully sequenced controls
#       -> List of the planned controls

# Output -> Summary table of number of controls
#==============================================================================================================

# SET UP 

library(tidyverse)
library(ggplot2)

rm(list=ls())

#===================================================================================================================
# READ IN DATA

# Planned list of control samples: 
planned_controls <- read.csv("~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Control_information/planned_controls.csv")

# per control type
table(planned_controls$Control_type)

# per plate
table(planned_controls$Extraction_plate)

# Read in the controls which were sequenced: 
seq_controls <- read_tsv("~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/Metabarcoding_Reads_per_control_summary_table_all_steps.tsv")

#====================================================================================================================
# CURATE THE DATASET

# Merge the seq information
all_controls_m <- left_join(planned_controls, seq_controls, by = "Sample_ID") %>% select(!Sampling_timepoint) # drop the duplicate column

# Add notes on the drop out stage
all_controls_m$Notes <- rep("Reads retained in control throughout processing", times = nrow(all_controls_m))
all_controls_m$Notes[is.na(all_controls_m$Final_reads_per_sample)|all_controls_m$Final_reads_per_sample== 0 ] <- "Dropped out during removal of non-fungal and unassigned reads"
all_controls_m$Notes[is.na(all_controls_m$Non_syn_comm)|all_controls_m$Non_syn_comm == 0 ] <- "Dropped out during removal of synthetic community (only applicable to non-synthetic community controls)"
all_controls_m$Notes[is.na(all_controls_m$Non_mito_chlor)|all_controls_m$Non_mito_chlor == 0 ] <- "Dropped out during removal of mitochondrial and chloroplast reads"
all_controls_m$Notes[is.na(all_controls_m$Non_chimeric)|all_controls_m$Non_chimeric == 0 ] <- "Dropped out during removal of chimeras"
all_controls_m$Notes[is.na(all_controls_m$Merged)|all_controls_m$Merged == 0 ] <- "Dropped out during merging"
all_controls_m$Notes[is.na(all_controls_m$Denoised)|all_controls_m$Denoised == 0 ] <- "Dropped out during denoising"
all_controls_m$Notes[is.na(all_controls_m$Filtered)|all_controls_m$Filtered == 0 ] <- "Dropped out during quality filtering"
all_controls_m$Notes[is.na(all_controls_m$Raw_reads)|all_controls_m$Raw_reads == 0 ] <- "Dropped out due to failed library preparation or sequencing"

# Count the catagories
table(all_controls_m$Notes)

#=====================================================================================================================
# OUTPUTS

# Create summary outputs from the paper

# Table of the controls which we have sequence data for 

sequenced_controls <- all_controls_m %>%
  filter(Notes == "Reads retained in control throughout processing") %>%
  count(Timepoint, Control_type) %>% 
  pivot_wider(names_from = Timepoint, values_from = n)
sequenced_controls$Total <- sequenced_controls$T1 + sequenced_controls$T2 + sequenced_controls$T3


write_tsv(all_controls_m, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Control_information/All_controls_information.tsv")

write_tsv(sequenced_controls, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Control_information/Sequenced_controls_summary_table.tsv")
