# Count reads for the RNA files

# Get the raw reads and the cleaned reads from the fastp output files

# -------------------------------------------------------------------------
# Parse fastp log files into a QC table => assisted by microsoft co-pilot
# -------------------------------------------------------------------------

library(stringr)
library(dplyr)
library(purrr)
library(readr)

# ---- CONFIG ----
log_dir <- "~/GitHub/PineBiomeDataPaper/RNA_data/logs"

# ---- HELPERS ----

extract_value <- function(text, pattern) {
  match <- str_match(text, pattern)
  if (is.na(match[1,2])) return(NA_real_)
  as.numeric(match[1,2])
}

extract_percent <- function(text, pattern) {
  match <- str_match(text, pattern)
  if (is.na(match[1,2])) return(NA_real_)
  as.numeric(match[1,2])
}

extract_string <- function(txt, pattern) {
  stringr::str_match(txt, pattern)[,2]
}

parse_fastp_log <- function(file) {
  
  sample_name <- str_remove(basename(file), "_fastp\\.log$")
  message("  Parsing sample: ", sample_name)
  
  txt <- read_file(file)
  
  tibble(
    sample = sample_name,
    
    # Reads / bases BEFORE
    reads_before = extract_value(
      txt,
      "Read1 before filtering:[\\s\\S]*?total reads:\\s*(\\d+)"
    ),
    bases_before = extract_value(
      txt,
      "Read1 before filtering:[\\s\\S]*?total bases:\\s*(\\d+)"
    ),
    
    # Reads / bases AFTER
    reads_after = extract_value(
      txt,
      "Read1 after filtering:[\\s\\S]*?total reads:\\s*(\\d+)"
    ),
    bases_after = extract_value(
      txt,
      "Read1 after filtering:[\\s\\S]*?total bases:\\s*(\\d+)"
    ),
    
    # Q30 percentages
    q30_before_read1 = extract_percent(
      txt,
      "Read1 before filtering:[\\s\\S]*?Q30 bases:.*?\\(([^%]+)%\\)"
    ),
    q30_after_read1 = extract_percent(
      txt,
      "Read1 after filtering:[\\s\\S]*?Q30 bases:.*?\\(([^%]+)%\\)"
    ),
    
    # Filtering results
    reads_passed_filter = extract_value(
      txt,
      "reads passed filter:\\s*(\\d+)"
    ),
    reads_failed_low_quality = extract_value(
      txt,
      "reads failed due to low quality:\\s*(\\d+)"
    ),
    reads_failed_too_short = extract_value(
      txt,
      "reads failed due to too short:\\s*(\\d+)"
    ),
    
    # Duplication and insert size
    duplication_rate = extract_percent(
      txt,
      "Duplication rate:\\s*([^%]+)%"
    ),
    insert_size_peak = extract_value(
      txt,
      "Insert size peak.*?:\\s*(\\d+)"
    )
  )
}

# ---- RUN ----

message("Searching for fastp log files...")

log_files <- list.files(
  log_dir,
  pattern = "_fastp\\.log$",
  full.names = TRUE
)

message("Found ", length(log_files), " fastp log file(s)")

message("Beginning parsing of log files...")

qc_table <- map_dfr(log_files, parse_fastp_log)

message("Finished parsing all samples")
message("Total samples parsed: ", nrow(qc_table))

# ---- PREVIEW ----
print(qc_table)



#=============================
# Import rRNA_filtered counts
#=============================
# Reads remaining after filtering rRNA 

filtered <- read.table("~/GitHub/PineBiomeDataPaper/RNA_data/rRNA_read_counts.txt", header = T)
filtered <- filtered  %>% dplyr::rename(sample = Sample)
# Merge with the qc table

t1 <- merge(qc_table, filtered, by = "sample")
str(t1)
write.csv(t1, "~/GitHub/PineBiomeDataPaper/RNA_data/temp_table.csv")
#========================
# Parse HISAT2 logs
#==================================

parse_hisat2_double_log <- function(file) {
  
  txt <- read_file(file)
  
  # ---- Extract sample name ----
  sample_name <- extract_string(
    txt,
    "Aligning reads from sample\\s+(\\S+?)\\s+to"
  )
  
  message("  Parsing sample: ", sample_name)
  
  # ---- Split into two HISAT blocks ----
  blocks <- str_extract_all(
    txt,
    "Aligning reads from sample[\\s\\S]+?overall alignment rate"
  )[[1]]
  
  if (length(blocks) < 2) {
    stop("Expected two HISAT logs in file")
  }
  
  # Normalize encoding
  host_txt <- str_replace_all(blocks[1], "&gt;", ">")
  nonhost_txt <- str_replace_all(blocks[2], "&gt;", ">")
  
  # ========================
  # STEP 1 — HOST FILTERING
  # ========================
  
  step1_input_pairs <- 
    extract_value(
      host_txt,
      "(\\d+) reads;"
    )
  
  
  host_concordant_0 <- extract_value(
    host_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned concordantly 0 times"
  )
  
  host_concordant_1 <- extract_value(
    host_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned concordantly exactly 1 time"
  )
  
  host_concordant_gt1 <- extract_value(
    host_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned concordantly >1 times"
  )
  
  host_discordant <- extract_value(
    host_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned discordantly 1 time"
  )
  
  host_unaligned <- extract_value(
    host_txt,
    "(\\d+) pairs aligned 0 times concordantly or discordantly"
  )
  
  host_total_pairs <- 
    host_concordant_0 +
    host_concordant_1 +
    host_concordant_gt1
  
  host_aligned <- host_concordant_1 + host_concordant_gt1
  host_disconcordant_and_unaligned <- host_unaligned + host_discordant
    
  # ========================
  # STEP 2 — NON-TARGET FILTERING
  # ========================
  
  step2_input_pairs <- 
    extract_value(
      nonhost_txt,
      "(\\d+) reads;"
    )
  
  nontarget_concordant_1 <- extract_value(
    nonhost_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned concordantly exactly 1 time"
  )
  
  nontarget_concordant_gt1 <- extract_value(
    nonhost_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned concordantly >1 times"
  )
  
  nontarget_discordant <- extract_value(
    nonhost_txt,
    "(\\d+) \\([0-9\\.]+%\\) aligned discordantly 1 time"
  )
  
  
  microbial_pairs <- extract_value(
    nonhost_txt,
    "(\\d+) pairs aligned 0 times concordantly or discordantly"
  ) # These are ones that had no kind of alignment to phix or our chloroplast/mito genomes
  
  # Your terminology
  Phix_and_organelle_reads <-
    step2_input_pairs - microbial_pairs
  # So to get ones with some kind of alignment to phix or our chloroplast/mito genomes we can do input - unaligned
  
  
  
  # ========================
  # QC CHECKS
  # ========================
  
  message("  Step1 input pairs: ", step1_input_pairs)
  message("  Host pairs: ", host_aligned)
  message("  Non-aligned pairs (concordant only)", host_unaligned)
  message("  Non Aligned pairs + host discordant reads", host_disconcordant_and_unaligned)
  message("  Step2 input pairs: ", step2_input_pairs)
  message("  Phix/organelle: ", Phix_and_organelle_reads)
  message("  Microbial: ", microbial_pairs)
  
  # ========================
  # OUTPUT
  # ========================
  
  tibble(
    sample = sample_name,
    
    # ---- STEP 1 ----
    step1_input_pairs = step1_input_pairs,
    host_aligned = host_aligned,
    host_unaligned = host_unaligned,
    host_disconcordant_and_unaligned = host_disconcordant_and_unaligned,
    step2_input_pairs = step2_input_pairs,
    Phix_and_organelle_concordant_reads = Phix_and_organelle_reads,
    Possible_microbial_pairs = microbial_pairs
  )
}


# ---- RUN ----

message("Searching for HISAT2 log files...")

log_files <- list.files(
  log_dir,
  pattern = "^Hisat.*\\.out$",
  full.names = TRUE
)

message("Found ", length(log_files), " HISAT2 log file(s)")

if (length(log_files) == 0) {
  stop("No HISAT2 log files found — check log_dir and filename pattern")
}

message("Beginning parsing of log files...")

align_table <- map_dfr(log_files, parse_hisat2_double_log)

message("Finished parsing all samples")
message("Total samples parsed: ", nrow(align_table))

str(align_table)

t2 <- merge(t1, align_table, by = "sample")

#========================
# Merge in tree metadata
#========================

# Split the sample names into timepoint and genotype

qc_table$timepoint <- str_extract(qc_table$sample, "^[^_]+")
qc_table$genotype  <- str_extract(qc_table$sample, "[^_]+$")
str(qc_table)                         

t3 <- merge(t2, qc_table)

genotype_metadata <- read.csv("~/GitHub/PineBiome-Metabarcoding/R/Data/genotype_classification_info.csv")                                  
str(genotype_metadata)

t4 <- left_join(t3, genotype_metadata)
str(t4)


#======================================
# Select columns and order for SI table
#======================================

SI_table2 <- data.frame(
    Sample_ID = t4$sample,
    Tree_genotype = t4$genotype,
    Tree_provenance_pop	= t4$population , 
    Tree_provenance_fam	= t4$family,
    Raw_reads	= t4$reads_before,
    Cleaned_reads = t4$reads_after,
    Filtered = t4$rRNAfiltered_R1_Reads,
    Host_reads = t4$host_aligned, 
    Nonhost_reads = t4$host_disconcordant_and_unaligned )

SI_table2$Percentage_host <- (SI_table2$Host_reads/SI_table2$Filtered)*100

write.csv(SI_table2, "~/GitHub/PineBiomeDataPaper/RNA_data/SI_table2.csv")

#=====================================
# Basic stats
#=====================================

# All reads
max(SI_table2$Raw_reads)
min(SI_table2$Raw_reads)
mean(SI_table2$Raw_reads)
median(SI_table2$Raw_reads)


# Filtered
max(SI_table2$Filtered)
min(SI_table2$Filtered)
mean(SI_table2$Filtered)
median(SI_table2$Filtered)

# Host reads
max(SI_table2$Host_reads)
min(SI_table2$Host_reads)
mean(SI_table2$Host_reads)
median(SI_table2$Host_reads)

# Non-host reads
max(SI_table2$Nonhost_reads)
min(SI_table2$Nonhost_reads)
mean(SI_table2$Nonhost_reads)
median(SI_table2$Nonhost_reads)

# Host perc
min(SI_table2$Percentage_host)
max(SI_table2$Percentage_host)
