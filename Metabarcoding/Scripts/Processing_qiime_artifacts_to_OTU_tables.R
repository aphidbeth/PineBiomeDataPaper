# Data paper: Curation of OTU tables

# Author: Beth Moore
# Created: 13/01/2026

#==========
# Set up 
#==========

library("tidyverse"); packageVersion("tidyverse")
library("qiime2R"); packageVersion("qiime2R")
library("phyloseq"); packageVersion("phyloseq")
library("microViz"); packageVersion("microviz")
library("microbiome"); packageVersion("microbiome")
library("ggpubr"); packageVersion("ggpubr")

# > library("tidyverse"); packageVersion("tidyverse")
# [1] ‘2.0.0’
# > library("qiime2R"); packageVersion("qiime2R")
# [1] ‘0.99.6’
# > library("phyloseq"); packageVersion("phyloseq")
# [1] ‘1.50.0’
# > library("microViz"); packageVersion("microviz")
# [1] ‘0.12.5’
# > library("microbiome"); packageVersion("microbiome")
# [1] ‘1.28.0’

#==============================
# IMPORT THE PHYSEQ OBJECTS:
#==============================

T1_ITS_physeq <- qza_to_phyloseq(features="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/table_T1_ITS-nochlor_nomito.qza",                  
                                 taxonomy="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/merged_taxonomy_T1_ITS.qza", 
                                 metadata="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/T1_ITS_metadata.tsv",
                                 tree = "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/tree_T1_ITS.qza"
)
T1_ITS_physeq

T2_ITS_physeq <- qza_to_phyloseq(features="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/table_T2_ITS-nochlor_nomito.qza",                  
                                 taxonomy="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/merged_taxonomy_T2_ITS.qza", 
                                 metadata="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/T2_ITS_metadata.tsv",
                                 tree = "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/tree_T2_ITS.qza"
)
T2_ITS_physeq


T3_ITS_physeq <- qza_to_phyloseq(features="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/table_T3_ITS-nochlor_nomito.qza",                  
                                 taxonomy="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/merged_taxonomy_T3_ITS.qza", 
                                 metadata="~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/T3_ITS_metadata.tsv",
                                 tree = "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Qiime_artifacts/tree_T3_ITS.qza"
)
T3_ITS_physeq



#========================
# Data curation steps
#========================

# Add a new column for contaminant categories
controls <- c("DNAblank", "Amp-Blank", "Idx-Blank", "Syn-Rep", "Syn-Extr")

add_contam_cat <- function(physeq_obj) {
  # Extract sample data as a data frame
  samp_df <- as(sample_data(physeq_obj), "data.frame")
  
  # Add contam_cat column
  samp_df$contam_cat <- ifelse(samp_df$population %in% controls,
                               samp_df$population,
                               "Tree-Sample")
  
  # Reassign updated sample data into the phyloseq object
  sample_data(physeq_obj) <- sample_data(samp_df)
  
  # Return modified phyloseq object
  return(physeq_obj)
}


# Add the categories to the physeq objects: 
T1_ITS_physeq <- add_contam_cat(T1_ITS_physeq)
T2_ITS_physeq <- add_contam_cat(T2_ITS_physeq)
T3_ITS_physeq <- add_contam_cat(T3_ITS_physeq)


#======================================================================
# Identify taxa in synthetic community controls 
#======================================================================

# Content of the synthetic community controls
# - Listeria monocytogenes: 12%
# - Pseudomonas aeruginosa: 12%
# - Bacillus subtilis: 12%
# - Escherichia coli: 12%
# - Salmonella enterica: 12%
# - Lactobacillus fermentum: 12%
# - Enterococcus faecalis: 12%
# - Staphylococcus aureus: 12%
# - Saccharomyces cerevisiae: 2% <===== FUNGAL TAXA
# - Cryptococcus neoformans: 2% <====== FUNGAL TAXA

known_comm <- c("Cryptococcus_neoformans","Cryptococcus Genus", "Saccharomyces_cerevisiae", "Saccharomyces Genus")


# What are the OTU IDs of the known community:
syncomm_T1 <- T1_ITS_physeq %>%
  tax_fix() %>%
  ps_filter(population %in% c("Syn-Rep", "Syn-Extr")) %>%
  ps_melt() %>%
  as_tibble() %>%
  select(OTU, population, Abundance, Kingdom:Species) %>%
  filter(Abundance > 0, (Species %in% known_comm ))

syncomm_T2 <- T2_ITS_physeq %>%
  tax_fix() %>%
  ps_filter(population %in% c("Syn-Rep", "Syn-Extr")) %>%
  ps_melt() %>%
  as_tibble() %>%
  select(OTU, population, Abundance, Kingdom:Species) %>%
  filter(Abundance > 0, (Species %in% known_comm ))

syncomm_T3 <- T3_ITS_physeq %>%
  tax_fix() %>%
  ps_filter(population %in% c("Syn-Rep", "Syn-Extr")) %>%
  ps_melt() %>%
  as_tibble() %>%
  select(OTU, population, Abundance, Kingdom:Species) %>%
  filter(Abundance > 0, (Species %in% known_comm ))

syncomm_T1
syncomm_T2
syncomm_T3
# We can see the Saccharomyces_cerevisiae strain is only classified to genus level here. 


#========================================================================================
# Remove the synthetic community from the samples which are not synthetic community controls
#========================================================================================

# Function to manually zero the syncomm OTU reads for these samples 
set_OTU_to_zero <- function(physeq){
  
  # Extract the synthetic community OTUS:
  syncomm <- physeq %>%
    tax_fix() %>%
    ps_filter(population %in% c("Syn-Rep", "Syn-Extr")) %>%
    ps_melt() %>%
    as_tibble() %>%
    select(OTU, population, Abundance, Kingdom:Species) %>%
    filter(Abundance > 0, (Species %in% known_comm ))
  
  OTU_list <- syncomm$OTU %>% unique()
  
  samples_to_modify <- rownames(sample_data(physeq))[
    !sample_data(physeq)$population %in% c("Syn-Rep", "Syn-Extr")
  ] 
  
  syncomm_samples <-rownames(sample_data(physeq))[
    sample_data(physeq)$population %in% c("Syn-Rep", "Syn-Extr")
  ] 
  
  cat("\nSamples being modified:\n")
  cat(samples_to_modify)
  
  cat("\nOTUs being removed:\n")
  cat(OTU_list)

  # Set these OTUs to zero
  otu_new <- otu_table(physeq)
  otu_new[OTU_list, samples_to_modify] <- 0
  
  
  # Check synthetic controls still have reps for these OTUS:
  cat("\nSynth rep check:\n")
  cat(otu_new[OTU_list, syncomm_samples])
  
  # Check other samples do not
  cat("\nOther sample check:\n")
  cat(otu_new[OTU_list, samples_to_modify])

  otu_table(physeq) <- otu_new
  
  return(physeq)
}

# Run on each dataset
T1_ITS_physeq_sc_depleted <- set_OTU_to_zero(T1_ITS_physeq)
T2_ITS_physeq_sc_depleted <- set_OTU_to_zero(T2_ITS_physeq)
T3_ITS_physeq_sc_depleted <- set_OTU_to_zero(T3_ITS_physeq)

#========================================================================================
# Remove any taxa not assigned to fungal kingdom. 
#========================================================================================

T1_ITS_physeq_no_unassigned <- T1_ITS_physeq_sc_depleted %>% subset_taxa(Kingdom == "Fungi")
T2_ITS_physeq_no_unassigned <- T2_ITS_physeq_sc_depleted %>% subset_taxa(Kingdom == "Fungi")
T3_ITS_physeq_no_unassigned <- T3_ITS_physeq_sc_depleted %>% subset_taxa(Kingdom == "Fungi")


#========================================================================================
# Remove two outlying samples
#========================================================================================
# Two samples are are outlying: one is massively overrepresented in read count and one
# has a taxonomic profile that suggests it is artifactual. These are removed.

outlying_GE1_sample <- T1_ITS_physeq_no_unassigned %>% plot_bar() + theme_bw() + theme(legend.position = "None", axis.text.x = element_text(angle=-90)) 
tiff("~/GitHub/PineBiomeDataPaper/Metabarcoding/Plots/outlying_GE1_sample.tiff", units = "in", height = 7, width = 17, res = 200)
print(outlying_GE1_sample)
dev.off()

outlying_CR3_sample <- T3_ITS_physeq_no_unassigned %>%  transform_sample_counts(., function(x) x / sum(x)) %>%
  plot_bar(fill="Class") + theme_bw() + theme(legend.position = "None", axis.text.x = element_text(angle=-90)) 
tiff("~/GitHub/PineBiomeDataPaper/Metabarcoding/Plots/outlying_CR3_sample.tiff", units = "in", height = 7, width = 17, res = 200)
print(outlying_CR3_sample)
dev.off()


T1_ITS_physeq_final <- T1_ITS_physeq_no_unassigned %>% ps_filter(.,  sample_names(T1_ITS_physeq_no_unassigned) != "GE1-7112-T1-ITS")
T2_ITS_physeq_final <- T2_ITS_physeq_no_unassigned 
T3_ITS_physeq_final <- T3_ITS_physeq_no_unassigned %>% ps_filter(.,  sample_names(T3_ITS_physeq_no_unassigned) != "CR3-7338-T3-ITS")


#========================================================================================
# OUTPUT OTU TABLES
# 
# Note: This represents a minimally filtered version of the dataset. Further decontamination
# and curation of the dataset is needed prior to downstream analysis. 
#========================================================================================

write_OTU_table <- function(timepoint, amplicon, suffix, physeq_obj, output_path) {
  # Output full OTU table
  OTU_df <- as.data.frame(otu_table(physeq_obj))
  OTU_list <- data.frame(OTU_string = rownames(OTU_df))
  final_OTU <- cbind(OTU_list, OTU_df)
  
  # Write to folder
  write.table(final_OTU, 
              file.path(output_path, paste0(timepoint, "_", amplicon, suffix, "_OTU_table.tsv")),
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t")
}

# Run for each dataset
write_OTU_table("T1", "ITS", "_min_clean", T1_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables")
write_OTU_table("T2", "ITS", "_min_clean", T2_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables")
write_OTU_table("T3", "ITS", "_min_clean", T3_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables")

# Note OTU id can be cross referenced with the seqs-and-taxonomy files to find the ITS sequence for a given id.

# I also output the RDS objects for easy use of the datasets in phyloseq
saveRDS(T1_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables/T1_ITS_min_clean_OTU.rds")
saveRDS(T2_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables/T2_ITS_min_clean_OTU.rds")
saveRDS(T3_ITS_physeq_final, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/OTU_tables/T3_ITS_min_clean_OTU.rds")

#=========================================================================================
# Summary statistics of read counts before and after this processing:

# Original counts - This is the raw data after import of the qimme objects, so it has been through
# quality filtering, denoising, merging, chimera removal and removal of mitochondrial reads, but
# it still retains unassigned reads and non-fungal reads. It also still has the cross contamination
# from the synthetic communties into the tree samples. 

# Export sample sums and ASV sums to a tsv for summary stats:
T1_mtchlor_depleted_read_counts <- sample_sums(T1_ITS_physeq)
T1_mtchlor_depleted_asv_counts <- apply((otu_table(T1_ITS_physeq)), 2, function(x) sum(x > 0))

T2_mtchlor_depleted_read_counts <- sample_sums(T2_ITS_physeq)
T2_mtchlor_depleted_asv_counts <- apply((otu_table(T2_ITS_physeq)), 2, function(x) sum(x > 0))

T3_mtchlor_depleted_read_counts <- sample_sums(T3_ITS_physeq)
T3_mtchlor_depleted_asv_counts <- apply((otu_table(T3_ITS_physeq)), 2, function(x) sum(x > 0))


mtchlor_depleted_reads_df <- data.frame(sample_id = c(names(T1_mtchlor_depleted_read_counts), names(T2_mtchlor_depleted_read_counts), names(T3_mtchlor_depleted_read_counts)),
                                        mtchl_depleted_reads = c(T1_mtchlor_depleted_read_counts, T2_mtchlor_depleted_read_counts, T3_mtchlor_depleted_read_counts),
                                        mtchl_depleted_asv = c(T1_mtchlor_depleted_asv_counts, T2_mtchlor_depleted_asv_counts, T3_mtchlor_depleted_asv_counts))

write_tsv(mtchlor_depleted_reads_df, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/mtchlor_depleted_reads_df.tsv")


# After removal of syn-comm from non syn-comm samples
# Export sample sums and ASV sums to a tsv for summary stats:
T1_syncomm_depleted_read_counts <- sample_sums(T1_ITS_physeq_sc_depleted)
T1_syncomm_depleted_asv_counts <- apply((otu_table(T1_ITS_physeq_sc_depleted)), 2, function(x) sum(x > 0))

T2_syncomm_depleted_read_counts <- sample_sums(T2_ITS_physeq_sc_depleted)
T2_syncomm_depleted_asv_counts <- apply((otu_table(T2_ITS_physeq_sc_depleted)), 2, function(x) sum(x > 0))

T3_syncomm_depleted_read_counts <- sample_sums(T3_ITS_physeq_sc_depleted)
T3_syncomm_depleted_asv_counts <- apply((otu_table(T3_ITS_physeq_sc_depleted)), 2, function(x) sum(x > 0))


syncomm_depleted_reads_df <- data.frame(sample_id = c(names(T1_syncomm_depleted_read_counts), names(T2_syncomm_depleted_read_counts), names(T3_syncomm_depleted_read_counts)),
                                           sc_depleted_reads = c(T1_syncomm_depleted_read_counts, T2_syncomm_depleted_read_counts, T3_syncomm_depleted_read_counts),
                                           sc_depleted_asv = c(T1_syncomm_depleted_asv_counts, T2_syncomm_depleted_asv_counts, T3_syncomm_depleted_asv_counts))

write_tsv(syncomm_depleted_reads_df, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/syncomm_depleted_reads_df.tsv")



# After removal of any unassigned and non-fungal reads

T1_no_unassigned_read_counts <- sample_sums(T1_ITS_physeq_no_unassigned)
T1_no_unassigned_asv_counts <- apply((otu_table(T1_ITS_physeq_no_unassigned)), 2, function(x) sum(x > 0))

T2_no_unassigned_read_counts <- sample_sums(T2_ITS_physeq_no_unassigned)
T2_no_unassigned_asv_counts <- apply((otu_table(T2_ITS_physeq_no_unassigned)), 2, function(x) sum(x > 0))

T3_no_unassigned_read_counts <- sample_sums(T3_ITS_physeq_no_unassigned)
T3_no_unassigned_asv_counts <- apply((otu_table(T3_ITS_physeq_no_unassigned)), 2, function(x) sum(x > 0))


no_unassigned_reads_df <- data.frame(sample_id = c(names(T1_no_unassigned_read_counts), names(T2_no_unassigned_read_counts), names(T3_no_unassigned_read_counts)),
                             no_unassigned_reads = c(T1_no_unassigned_read_counts, T2_no_unassigned_read_counts, T3_no_unassigned_read_counts),
                             no_unassigned_asv = c(T1_no_unassigned_asv_counts, T2_no_unassigned_asv_counts, T3_no_unassigned_asv_counts))

write_tsv(no_unassigned_reads_df, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/no_unassigned_reads_df.tsv")

# After removal of any unassigned and non-fungal reads

T1_final_read_counts <- sample_sums(T1_ITS_physeq_final)
T1_final_asv_counts <- apply((otu_table(T1_ITS_physeq_final)), 2, function(x) sum(x > 0))

T2_final_read_counts <- sample_sums(T2_ITS_physeq_final)
T2_final_asv_counts <- apply((otu_table(T2_ITS_physeq_final)), 2, function(x) sum(x > 0))

T3_final_read_counts <- sample_sums(T3_ITS_physeq_final)
T3_final_asv_counts <- apply((otu_table(T3_ITS_physeq_final)), 2, function(x) sum(x > 0))


final_reads_df <- data.frame(sample_id = c(names(T1_final_read_counts), names(T2_final_read_counts), names(T3_final_read_counts)),
                                   final_reads = c(T1_final_read_counts, T2_final_read_counts, T3_final_read_counts),
                                   final_asv = c(T1_final_asv_counts, T2_final_asv_counts, T3_final_asv_counts))

write_tsv(final_reads_df, "~/GitHub/PineBiomeDataPaper/Metabarcoding/Data/Read_stats/final_reads_df.tsv")

#================================================================================================================================
# Comparison of technical replicates

# There are several types of control in this study: negative controls, positive controls
# and techinical replicates of a single pine tree. .

subset_physeq_controls <- function(physeq_obj) {
  all_controls <- physeq_obj %>%
    ps_filter(experimental_factor != "Pinus sylvestris needle phyllosphere")
  
  all_trees <- physeq_obj %>%
    ps_filter(experimental_factor == "Pinus sylvestris needle phyllosphere")
  
  neg_controls <- all_controls %>%
    ps_filter(experimental_factor == "control [EFO_0001461]")
  
  pos_controls <- all_controls %>%
    ps_filter(experimental_factor == "replicate [EFO_0000683]")
  
  tech_reps <- all_controls %>%
    ps_filter(experimental_factor == "Replicate Pinus sylvestris needle phyllosphere")
  
  list(
    all_controls = all_controls,
    all_trees = all_trees,
    neg_controls = neg_controls,
    pos_controls = pos_controls,
    tech_reps = tech_reps
  )
}

T1_ITS_physeq_subsets <- subset_physeq_controls(T1_ITS_physeq_final)
T2_ITS_physeq_subsets <- subset_physeq_controls(T2_ITS_physeq_final)
T3_ITS_physeq_subsets <- subset_physeq_controls(T3_ITS_physeq_final)
# Access the different subsets using the $ index


# Create colour map:
T1_classes <- unique(as.character(tax_table(T1_ITS_physeq_subsets$tech_reps)[, "Class"]))
T2_classes <- unique(as.character(tax_table(T2_ITS_physeq_subsets$tech_reps)[, "Class"]))
T3_classes <- unique(as.character(tax_table(T3_ITS_physeq_subsets$tech_reps)[, "Class"]))
T1_only <- setdiff(T1_classes, union(T2_classes, T3_classes))
T2_only <- setdiff(T2_classes, union(T1_classes, T3_classes))
T3_only <- setdiff(T3_classes, union(T1_classes, T2_classes))
shared_all <- Reduce(intersect, list(T1_classes, T2_classes, T3_classes)) %>% sort()# 21 shared classes

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "darkgrey")

# Create named vector for shared classes
T1_cols <- c(setNames(my_colors, shared_all), setNames(rep("black", times = length(T1_only)),T1_only))
T2_cols <- c(setNames(my_colors, shared_all), setNames(rep("black", times = length(T2_only)),T2_only))
T3_cols <- c(setNames(my_colors, shared_all), setNames(rep("black", times = length(T3_only)),T3_only))

T1_tech_rep_bar_rel <- T1_ITS_physeq_subsets$tech_reps %>%  transform_sample_counts(., function(x) x / sum(x)) %>%  plot_bar(fill = "Class") + ggtitle("T1") + ylab("Relative Abundance")
T2_tech_rep_bar_rel <- T2_ITS_physeq_subsets$tech_reps %>%  transform_sample_counts(., function(x) x / sum(x)) %>%  plot_bar(fill = "Class") + ggtitle("T2") + ylab("Relative Abundance")
T3_tech_rep_bar_rel <- T3_ITS_physeq_subsets$tech_reps %>%  transform_sample_counts(., function(x) x / sum(x)) %>%  plot_bar(fill = "Class") + ggtitle("T3")+ ylab("Relative Abundance")


datapaper_A <- T1_tech_rep_bar_rel  + geom_bar(aes(color=Class, fill=Class), stat="identity") + theme(legend.position = "None") + scale_color_manual(values = T1_cols) + scale_fill_manual(values = T1_cols) 

datapaper_B <-  T2_tech_rep_bar_rel  + geom_bar(aes(color=Class, fill=Class), stat="identity") + theme(legend.position = "None") + scale_color_manual(values = T2_cols) + scale_fill_manual(values = T2_cols) 

datapaper_C <- T3_tech_rep_bar_rel  + geom_bar(aes(color=Class, fill=Class), stat="identity") + theme(legend.position = "None") + scale_color_manual(values = T3_cols) + scale_fill_manual(values = T3_cols) 

tech_reps_datapaper <- ggarrange(datapaper_A,
                                 datapaper_B,
                                 datapaper_C,
                                 nrow = 1, ncol = 3, common.legend = T, 
                                 legend = "bottom", labels = "AUTO" )

# Write to file for data paper

tiff("~/GitHub/PineBiomeDataPaper/Metabarcoding/Plots/Technical_Replicates.tiff", units = "in", height = 7, width = 9, res = 200)
print(tech_reps_datapaper)
dev.off()
