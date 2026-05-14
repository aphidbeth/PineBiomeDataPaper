# COLLATION OF READ AND OTU FILTERING ACROSS SAMPLES IN THE PINEBIOME DATA PAPER

# Author: Beth Moore
# Updated: 15/01/2026

# Input -> files from HPC counting number of reads in each sample file across the pipeline steps + 
# phyloseq objects for reads in samples after filtering

# Output -> A single table per timepoint with the sample ID, the tree ID and read counts across pipeline stages
#==============================================================================================================

# SET UP 

library(tidyverse)
library(ggplot2)
library(plotly)

rm(list=ls())

runs <- c("T1_run1_new", "T1_run2_new", "T2_run3_new", "T2_run4_new", "T3_run5_new", "T3_run6_new")
amplicons <- c("ITS") # just ITS for now

#====================================================================================================

# 1) RAW READS OFF SEQ MACHINE

my_list <- list()

for (run in runs){
  
  for (amplicon in amplicons){
    
    runname <- paste("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/raw_read_counts_", paste(run, amplicon, sep = "_"), ".tsv", sep = "")
    
    my_df <- read.table(runname)
    
    my_list[[runname]] <- my_df
  }
  
}

raw_reads_df <- do.call(rbind, my_list)

names(raw_reads_df) <- c("filename", "seq_run", "amplicon", "raw_reads")

# From these we want to keep the filename and the raw_reads information
selected_cols <- raw_reads_df %>% select(filename, raw_reads)

#====================================================================================================

# 2) PER SAMPLE READS AFTER FASTP POLYG TAIL AND SHORT READ LENGTH FILTERING

my_list <- list()

for (run in runs){
  
  for (amplicon in amplicons){
    
    runname <- paste("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/post_fastp_filtering_read_counts_", paste(run, amplicon, sep = "_"), ".tsv", sep = "")
    
    my_df <- read.table(runname)
    
    my_list[[runname]] <- my_df
  }
  
}

fastp_reads_df <- do.call(rbind, my_list)

names(fastp_reads_df) <- c("filename", "seq_run", "amplicon", "reads_post_fastp_filters")


#====================================================================================================

# 3) PER SAMPLE READS AFTER ITEXPRESS FILTERING (ITS reads impacted only)

my_list <- list()

for (run in runs){
  
  for (amplicon in amplicons){
    
    runname <- paste("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/post_ITSxpress_read_counts_", paste(run, amplicon, sep = "_"), ".tsv", sep = "")
    
    my_df <- read.table(runname, skip = 1)
    
    my_list[[runname]] <- my_df
  }
  
}

itsxpress_reads_df <- do.call(rbind, my_list)

names(itsxpress_reads_df) <- c("sample_id", "post_itsexpress_read_counts", "rev_counts") # post_itsexpress_read_counts = forward counts

itsxpress_reads_df <- itsxpress_reads_df[,1:2] 

#====================================================================================================

# 4) PER SAMPLE READS ACROSS DADA2 FILTERING STEPS

my_list <- list()

for (run in runs){
  
  for (amplicon in amplicons){
    
    runname <- paste("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/denoising_stats_", paste(run, amplicon, sep = "_"), ".tsv", sep = "")
    
    my_df <- read.table(runname, skip = 2) # removes metadata type that qiime2 automatically outputs
    
    my_df$seq_run <- rep(run, times = nrow(my_df))
    
    my_list[[runname]] <- my_df
  }
  
}

dada2_denoise_df <- do.call(rbind, my_list)

names(dada2_denoise_df) <- c("sample_id", "reads_post_cutadapt_filter", "reads_post_qual_filter", "p_reads_post_qual_filter", "reads_post_denoise", "reads_post_merge", "p_reads_merged", "reads_post_chimera_filter", "p_reads_retained", "seq_run")

#====================================================================================================

# 7) PER SAMPLE READS AFTER CLASSIFICATION

my_list <- list()

for (run in runs){
  
  for (amplicon in amplicons){
    
    runname <- paste("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/freq_table_seq_", paste(run, amplicon, sep = "_"), ".tsv", sep = "")
    
    my_df <- read.table(runname, skip = 2) # removes metadata type that qiime2 automatically outputs
    
    my_df$seq_run <- rep(run, times = nrow(my_df))
    
    my_list[[runname]] <- my_df
  }
  
}

post_class_df <- do.call(rbind, my_list)

post_class_df[, 2:3] <- lapply(post_class_df[, 2:3], function(x) as.numeric(gsub(",", "", x)))

names(post_class_df ) <- c("sample_id", "reads_post_classifcation", "ASVs_post_classifcation", "seq_run")

#====================================================================================================

# 6) PER SAMPLE READS & ASVS AFTER FILTERING MITOCHONDRIA & CHLOROPLAST

mtchlor_reads <- read.table("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/mtchlor_depleted_reads_df.tsv", header = T)
str(mtchlor_reads)

#====================================================================================================

# 7) Integrate the number dropped after removing contaminants

sc_reads <- read.table("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/syncomm_depleted_reads_df.tsv", header = T)
str(sc_reads)

#====================================================================================================

# 8) Integrate the number dropped after removal of non-fungal reads and outliers
no_unassigned_reads <- read.table("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/no_unassigned_reads_df.tsv", header = T)
str(no_unassigned_reads)


#====================================================================================================
# 9) Integrate the number dropped after removal of outliers
final_reads <- read.table("~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/final_reads_df.tsv", header = T)
str(final_reads)


#====================================================================================================
# CURATION

# Remove the sample number suffix from some datasets ()

raw_reads_df <-  raw_reads_df %>% 
  mutate(sample_id = str_remove_all(filename, regex("_R1_001.fastq.gz"))) %>% 
  mutate(sample_id = str_remove_all(sample_id, regex("_S\\d{1,3}")))

itsxpress_reads_df <- itsxpress_reads_df %>% mutate(sample_id = str_remove_all(sample_id, regex("_S\\d{1,3}"))) 

dada2_denoise_df <- dada2_denoise_df %>% mutate(sample_id = str_remove_all(sample_id, regex("_S\\d{1,3}"))) 

post_class_df <- post_class_df %>% mutate(sample_id = str_remove_all(sample_id, regex("_S\\d{1,3}"))) 



#==============================================================
# MERGING

t1df <- left_join(raw_reads_df, itsxpress_reads_df, by = c("sample_id"))
t2df <- left_join(t1df, dada2_denoise_df, by = c("sample_id", "seq_run"))
t3df <- left_join(t2df, post_class_df, by = c("sample_id", "seq_run"))
t4df <- left_join(t3df, mtchlor_reads, by = c("sample_id"))
t5df <- left_join(t4df, sc_reads, by = c("sample_id"))
t6df <- left_join(t5df, no_unassigned_reads, by = c("sample_id"))
final_df <- left_join(t6df, final_reads, by = c("sample_id"))

final_df <- final_df %>% mutate(family = str_extract(sample_id, regex("^\\w{2}\\d{1,2}"))) 
final_df <- final_df %>% mutate(pop = str_extract(family, regex("^\\w{2}"))) %>%
  mutate(amplicon  = str_extract(sample_id, regex("(?<=-).{2}S"))) %>%
  mutate(genotype = str_extract(sample_id, regex("(?<=\\w{2}\\d{1}-).*?(?=\\-T\\d{1}-)"))) %>%
  mutate(timepoint = str_extract(sample_id, regex("(?<=-)T\\d{1}(?=-)")))
final_df$sample_type <- rep("field", times = nrow(final_df))
final_df$sample_type[final_df$pop == "EP"] <- "Control"

View(final_df)

#=================================================================
# RENAMING FOR OUTPUT

final_df_selected <- final_df %>% filter(., sample_type == "field") %>% 
  select(Sample_ID = sample_id,
         Tree_genotype = genotype,
         Tree_provenance_pop = pop,
         Tree_provenance_fam = family,
         Sampling_timepoint = timepoint,
         Raw_reads = raw_reads,
         Filtered = reads_post_qual_filter,
         Denoised = reads_post_denoise,
         Merged = reads_post_merge,
         Non_chimeric = reads_post_chimera_filter,
         Non_mito_chlor = mtchl_depleted_reads,
         Non_syn_comm = sc_depleted_reads,
         No_unassigned = no_unassigned_reads,
         Final_reads_per_sample = final_reads
  )

View(final_df_selected)

final_df_selected <- final_df_selected[order(final_df_selected$Sampling_timepoint, final_df_selected$Sample_ID),]


# Write to table 

write_tsv(final_df_selected, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/Metabarcoding_Reads_per_sample_summary_table_all_steps.tsv")


# Output controls info too 
controls_df_selected <- final_df %>% filter(., sample_type != "field") %>% 
  select(Sample_ID = sample_id,
         Sampling_timepoint = timepoint,
         Raw_reads = raw_reads,
         Filtered = reads_post_qual_filter,
         Denoised = reads_post_denoise,
         Merged = reads_post_merge,
         Non_chimeric = reads_post_chimera_filter,
         Non_mito_chlor = mtchl_depleted_reads,
         Non_syn_comm = sc_depleted_reads,
         Final_reads_per_sample = final_reads
  )

write_tsv(controls_df_selected, "~/Github/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/Metabarcoding_Reads_per_control_summary_table_all_steps.tsv")
#===========================================================================================================================