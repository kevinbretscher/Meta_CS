# Install the rentrez package if it's not already installed
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

# Load the rentrez package
library(rentrez)


# Function to fetch specific regions of nucleotide sequences based on a data frame of GenBank IDs and coordinates
fetch_sequence_regions <- function(region_df, file_output = NULL) {
  # Check if input data frame has the required columns
  if (!all(c("GenBank_ID", "Start", "End") %in% colnames(region_df))) {
    stop("Input data frame must contain 'GenBank_ID', 'Start', and 'End' columns.")
  }
  
  # Initialize an empty list to store the sequences
  sequences <- list()
  
  # Loop through each row of the data frame to fetch the sequence regions
  for (i in 1:nrow(region_df)) {
    # Extract GenBank ID, start, and end positions for the current row
    genbank_id <- region_df$GenBank_ID[i]
    start <- region_df$Start[i]
    end <- region_df$End[i]
    
    cat("Fetching sequence for GenBank ID:", genbank_id, 
        "from position", start, "to", end, "\n")
    
    # Fetch the specific region of the nucleotide sequence in fasta format
    sequence_data <- entrez_fetch(db = "nucleotide", 
                                  id = genbank_id, 
                                  rettype = "fasta", 
                                  seq_start = start, 
                                  seq_stop = end)
    
    # Store the fetched sequence with an identifier including the range
    sequences[[paste0(genbank_id, "_", start, "-", end)]] <- sequence_data
  }
  
  # Optionally write all sequences to a single file
  if (!is.null(file_output)) {
    cat("Saving sequences to file:", file_output, "\n")
    
    # Combine all sequences into one string
    combined_sequences <- paste(unlist(sequences), collapse = "\n")
    
    # Write the combined sequences to a single file
    write(combined_sequences, file = file_output)
  }
  
  # Return the fetched sequences as a named list
  return(sequences)
}

install.packages("tidyr")
install.packages("dplyr")

library(tidyr)
library(dplyr)


Silva_to_Genebank <- function(L) {
  df <- data.frame(ID = L)
  df2 <- df %>% separate(ID, c("GenBank_ID", "Start","End"))
  return(df2)
}

remove_after_first_dot <- function(input_string) {
  # Use sub function to replace everything after the first dot with an empty string
  result <- sub("\\..*", "", input_string)
  return(result)
}


# # Example usage:
# # Create a data frame with GenBank IDs and specific regions to download
# region_df <- data.frame(
#   GenBank_ID = c("AUOQ02000015", "FPLS01049422"),
#   Start = c(13719, 16),
#   End = c(15183, 1485)
# )
# 
# # Call the function to fetch sequence regions and save to a single FASTA file
# sequences <- fetch_sequence_regions(region_df, file_output = "specified_regions_sequences.fasta")
# 
# 
# #Download all sphingo sequences
# all_sphingo <- readRDS("all_sphingo.rds")
# all_sphingo <- Silva_to_Genebank(all_sphingo)
# sequences <- fetch_sequence_regions(all_sphingo, file_output = "all_sphino.fasta")
# 
# # Print the fetched sequences
# print(sequences)

#Get IDS from enriched OTUs

enriched_OTU <- read.csv("../count_OTUs/input/Enriched_OTUs_Metagenomeseq.csv")
enriched_OTU_ANCOM <- read.csv("../count_OTUs/input/Enriched_OTUs_ANCOM.csv")
enriched_OTU_ANCOM <- enriched_OTU_ANCOM %>% filter(q_C_SS < 0.1)

depleted_OTU <- read.csv("../count_OTUs/input/depleted_OTUs_Metagenomeseq.csv")
depleted_OTU_ANCOM <- read.csv("../count_OTUs/input/depleted_OTUs_ANCOM.csv")
depleted_OTU_ANCOM <- depleted_OTU_ANCOM %>% filter(q_C_SS < 0.1)

#filter

Filter <- "Chitinophagaceae"

#filter results

enriched_OTU_filtered <- enriched_OTU %>% filter(Family == Filter)
enriched_OTU_ANCOM_filtered <- enriched_OTU_ANCOM %>% filter(Family == Filter)
depleted_OTU_filtered <- depleted_OTU %>% filter(Family == Filter)
depleted_OTU_ANCOM_filtered <- depleted_OTU_ANCOM %>% filter(Family == Filter)

full_join(table(enriched_OTU_filtered$Genus), table(enriched_OTU_ANCOM_filtered$Genus),
      table(depleted_OTU_filtered$Genus), table(depleted_OTU_ANCOM_filtered$Genus))
# Convert tables to data frames
df1 <- as.data.frame(table(enriched_OTU_filtered$Genus))
df2 <- as.data.frame(table(enriched_OTU_ANCOM_filtered$Genus))
df3 <- as.data.frame(table(depleted_OTU_filtered$Genus))
df4 <- as.data.frame(table(depleted_OTU_ANCOM_filtered$Genus))

# Rename columns
colnames(df1) <- c("Genus", "Enriched_OTU_metagenomeSeq")
colnames(df2) <- c("Genus", "Enriched_OTU_ANCOM")
colnames(df3) <- c("Genus", "Depleted_OTU_metagenomeSeq")
colnames(df4) <- c("Genus", "Depleted_OTU_ANCOM")

# Merge all data frames using full_join
merged_df <- df1 %>%
  full_join(df2, by = "Genus") %>%
  full_join(df3, by = "Genus") %>%
  full_join(df4, by = "Genus") %>%
  replace(is.na(.), 0)  # Replace NAs with 0

#convert id to silva db id

enriched_OTU_filtered_ID <- Silva_to_Genebank(unique(enriched_OTU_filtered$Taxa))
enriched_OTU_ANCOM_filtered_ID <- Silva_to_Genebank(unique(enriched_OTU_ANCOM_filtered$taxon))
depleted_OTU_ID <- Silva_to_Genebank(unique(depleted_OTU_filtered$Taxa))
depleted_OTU_ANCOM_ID <- Silva_to_Genebank(unique(depleted_OTU_ANCOM_filtered$taxon))

  ### extract all ids
enriched_OTU_ID_all <- Silva_to_Genebank(unique(c(enriched_OTU$Taxa, enriched_OTU_ANCOM$taxon)))
depleted_OTU_ID_all <- Silva_to_Genebank(unique(c(depleted_OTU$Taxa, depleted_OTU_ANCOM$taxon)))

#combine

ID_comb <- rbind(enriched_OTU_filtered_ID, enriched_OTU_ANCOM_filtered_ID, depleted_OTU_ID,depleted_OTU_ANCOM_ID)
ID_comb <- ID_comb %>% distinct()

  ### combine all ids
ID_comb_all <- rbind(enriched_OTU_ID_all, depleted_OTU_ID_all)
ID_comb_all <- ID_comb_all %>% distinct()


#fetch

fetch_sequence_regions(ID_comb, file_output = "Chitinopagacea_new.fasta")

  ### fetch all
fetch_sequence_regions(ID_comb_all, file_output = "all_sig_q0.1/all_sig_OTUs.fasta")


# Convert tables to data frames
df_all_1 <- as.data.frame(table(enriched_OTU$Genus))
df_all_2 <- as.data.frame(table(enriched_OTU_ANCOM$Genus))
df_all_3 <- as.data.frame(table(depleted_OTU$Genus))
df_all_4 <- as.data.frame(table(depleted_OTU_ANCOM$Genus))

# Rename columns
colnames(df_all_1) <- c("Genus", "Enriched_OTU_metagenomeSeq")
colnames(df_all_2) <- c("Genus", "Enriched_OTU_ANCOM")
colnames(df_all_3) <- c("Genus", "Depleted_OTU_metagenomeSeq")
colnames(df_all_4) <- c("Genus", "Depleted_OTU_ANCOM")

# Merge all data frames using full_join
merged_df_all_ <- df_all_1 %>%
  full_join(df_all_2, by = "Genus") %>%
  full_join(df_all_3, by = "Genus") %>%
  full_join(df_all_4, by = "Genus") %>%
  replace(is.na(.), 0)  # Replace NAs with 0

Filter <- "Pseudomonas"

#filter results

enriched_OTU_filtered <- enriched_OTU %>% filter(Genus == Filter)
enriched_OTU_ANCOM_filtered <- enriched_OTU_ANCOM %>% filter(Genus == Filter)
depleted_OTU_filtered <- depleted_OTU %>% filter(Genus == Filter)
depleted_OTU_ANCOM_filtered <- depleted_OTU_ANCOM %>% filter(Genus == Filter)

#convert id to silva db id

enriched_OTU_filtered_ID <- Silva_to_Genebank(unique(enriched_OTU_filtered$Taxa))
enriched_OTU_ANCOM_filtered_ID <- Silva_to_Genebank(unique(enriched_OTU_ANCOM_filtered$taxon))
depleted_OTU_ID <- Silva_to_Genebank(unique(depleted_OTU_filtered$Taxa))
depleted_OTU_ANCOM_ID <- Silva_to_Genebank(unique(depleted_OTU_ANCOM_filtered$taxon))

setdiff(enriched_OTU_ANCOM_filtered_ID$GenBank_ID, depleted_OTU_ANCOM_ID$GenBank_ID)
setdiff(depleted_OTU_ANCOM_ID$GenBank_ID, enriched_OTU_ANCOM_filtered_ID$GenBank_ID)

#combine
ID_comb <- rbind(enriched_OTU_ANCOM_filtered_ID, depleted_OTU_ANCOM_ID)
ID_comb <- ID_comb %>% distinct()

fetch_sequence_regions(ID_comb, file_output = "Pseudomonas/pseudo_ancom_0.1.fasta")

## Final analysis for the manuscript
# only extract the sequences from OTUs belonging to core genera (ANCOM-BC2, FC > 2, FDR adjusted P < 0.1)
priorized_taxa  <- c("Bacillus", "Sphingomonas", "Chitinophaga", "Pseudomonas", "Streptomyces", "Flavobacterium")
for (taxa in priorized_taxa){
  print(taxa)
  OTU_enriched <- core_enriched_ANCOM$taxon[core_enriched_ANCOM$Genus == taxa]
  print(paste("#enriched: ", length(unique(OTU_enriched))))
  OTU_depleted <- core_depleted_ANCOM$taxon[core_depleted_ANCOM$Genus == taxa]
  print(paste("#depleted: ", length(unique(OTU_depleted))))
  OTU_IDs <- Silva_to_Genebank(unique(c(OTU_enriched, OTU_depleted)))
  fetch_sequence_regions(OTU_IDs, file_output = paste("final/", taxa, "_sig_OTU.fasta", sep = ""))
}

priorized_taxa  <- c("Novosphingobium")
for (taxa in priorized_taxa){
  print(taxa)
  OTU_enriched <- core_enriched_ANCOM$taxon[core_enriched_ANCOM$Genus == taxa]
  print(paste("#enriched: ", length(unique(OTU_enriched))))
  OTU_depleted <- core_depleted_ANCOM$taxon[core_depleted_ANCOM$Genus == taxa]
  print(paste("#depleted: ", length(unique(OTU_depleted))))
  OTU_IDs <- Silva_to_Genebank(unique(c(OTU_enriched, OTU_depleted)))
  fetch_sequence_regions(OTU_IDs, file_output = paste("final/", taxa, "_sig_OTU.fasta", sep = ""))
}

vsearch_all <- read.delim("final/blast6_all.txt", header = F)
colnames(vsearch_all)[1] <- "OTU"
colnames(vsearch_all)[2] <- "genome"
colnames(vsearch_all)[3] <- "similarity"
vsearch_all$label_renamed <- sub("\\..*", "", vsearch_all$OTU)
write.table(unique(vsearch_all$genome),"final/closest_ref_genomes/genome_id.txt", quote = F, row.names = F, col.names = F)

vsearch_all <- vsearch_all[,c(1,2,3,13)]
#vsearch_all_annot_enriched <- left_join(vsearch_all, core_enriched_ANCOM, by = "label_renamed")
#vsearch_all_annot_depleted <- left_join(vsearch_all, core_depleted_ANCOM, by = "label_renamed")

# add info if enriched/depleted
vsearch_all$ancom <- "unknown"
vsearch_all$ancom [vsearch_all$label_renamed %in% core_enriched_ANCOM$label_renamed] <- "enriched"
vsearch_all$ancom [vsearch_all$label_renamed %in% core_depleted_ANCOM$label_renamed] <- "depleted"
vsearch_all$genome <- sub("RS_", "", vsearch_all$genome)
vsearch_all$genome <- sub("GB_", "", vsearch_all$genome)

# Select only relevant columns from both ANCOM dataframes
depleted_genus <- core_depleted_ANCOM[, c("label_renamed", "Genus")]
enriched_genus <- core_enriched_ANCOM[, c("label_renamed", "Genus")]

# Combine genus information from both enriched and depleted ANCOM dataframes
all_genus_info <- unique(rbind(depleted_genus, enriched_genus))

# Merge genus information into vsearch_all
vsearch_all <- merge(vsearch_all, all_genus_info, by = "label_renamed")

# Count the number of unique 'ancom' values per genome
genome_status <- aggregate(ancom ~ genome, data = vsearch_all, function(x) length(unique(x)))

# Identify genomes that are both enriched and depleted
conflict_genomes <- genome_status$genome[genome_status$ancom > 1]

# Add a new column to mark conflicts
vsearch_all$conflict <- ifelse(vsearch_all$genome %in% conflict_genomes, "conflict", "unique")

# label where genomes are conflicting
vsearch_all$ancom [vsearch_all$conflict == "conflict"] <- "conflict"

# Drop the conflict column if not needed
#vsearch_filtered <- vsearch_filtered[, !colnames(vsearch_filtered) %in% "conflict"]
vsearch_all_annot <- unique(vsearch_all[,c(3,5,6)])
rownames(vsearch_all_annot) <- vsearch_all_annot$genome
vsearch_all_annot <- vsearch_all_annot[,-1]
#vsearch_all_annot_genus <- unique(vsearch_all[,c(3,5,6)])
#rownames(vsearch_all_annot_genus) <- vsearch_all_annot_genus$genome


