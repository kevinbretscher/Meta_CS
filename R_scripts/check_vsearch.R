Bacillus_vsearch <- read.delim("Bacillus/vsearch/blast6_all.txt", header = F)
colnames(Bacillus_vsearch)[1] <- "OTU"
colnames(Bacillus_vsearch)[2] <- "genome"
colnames(Bacillus_vsearch)[3] <- "similarity"
Bacillus_vsearch_endo <- read.delim("Bacillus/vsearch/blast6_endo.txt", header = F)
colnames(Bacillus_vsearch_endo)[1] <- "OTU"
colnames(Bacillus_vsearch_endo)[2] <- "genome"
colnames(Bacillus_vsearch_endo)[3] <- "similarity"
Bacillus_vsearch_dsmz <- read.delim("Bacillus/vsearch/blast6.txt", header = F)
colnames(Bacillus_vsearch_dsmz)[1] <- "OTU"
colnames(Bacillus_vsearch_dsmz)[2] <- "genome"
colnames(Bacillus_vsearch_dsmz)[3] <- "similarity"

# extract enriched and depleted Sphingomonas OTUs
enriched_bacillus_ancom <- unique(na.omit(enriched_OTU_ANCOM$taxon[enriched_OTU_ANCOM$Genus == "Bacillus"]))
depleted_bacillus_ancom <- unique(na.omit(depleted_OTU_ANCOM$taxon[depleted_OTU_ANCOM$Genus == "Bacillus"]))

# remove the conflicted OTUs (identified as both enriched and depleted)
enriched_bacillus_ancom <- setdiff(enriched_bacillus_ancom, depleted_bacillus_ancom)
depleted_bacillus_ancom <- setdiff(depleted_bacillus_ancom, enriched_bacillus_ancom)

# rename the IDs to match properly
enriched_bacillus_ancom <- sub("\\..*", "", enriched_bacillus_ancom)
depleted_bacillus_ancom <- sub("\\..*", "", depleted_bacillus_ancom)
Bacillus_vsearch$OTU <- sub("\\..*", "", Bacillus_vsearch$OTU)

# add info if enriched/depleted
Bacillus_vsearch$ancom <- "unknown"
Bacillus_vsearch$ancom [Bacillus_vsearch$OTU %in% enriched_bacillus_ancom] <- "enriched"
Bacillus_vsearch$ancom [Bacillus_vsearch$OTU %in% depleted_bacillus_ancom] <- "depleted"
Bacillus_conflicting_genomes <- Bacillus_vsearch %>%
  group_by(genome) %>%
  summarise(unique_ancom = n_distinct(ancom), 
            ancom_classes = paste(unique(ancom), collapse = ", ")) %>%
  filter(unique_ancom > 1)  # Keep only genomes that have both "enriched" & "depleted"

write.table(vsearch_all, file= "final/vsearch_results_annotated.txt",  row.names = F, col.names = T, quote = F)

