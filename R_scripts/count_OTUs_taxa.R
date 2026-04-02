install.packages("dplyr")
library(dplyr)
install.packages("reshape2")
library(reshape2)
install.packages("ggplot2")
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
library(scales)  # Load scales for formatting
install.packages("ggVennDiagram")
library(ggVennDiagram)
enriched_OTUs_ancom <- read.csv("input/Enriched_OTUs_ANCOM.csv")
enriched_OTUs_metagenomeseq <- read.csv("input/Enriched_OTUs_Metagenomeseq.csv")

depleted_OTUs_ancom <- read.csv("input/depleted_OTUs_ANCOM.csv")
depleted_OTUs_metagenomeseq <- read.csv("input/depleted_OTUs_Metagenomeseq.csv")

union(unique(enriched_OTUs_ancom$study),unique(depleted_OTUs_ancom$study))
union(unique(enriched_OTUs_metagenomeseq$study), unique(depleted_OTUs_metagenomeseq$study))


count_enriched_ancom <- enriched_OTUs_ancom %>%
  group_by(study, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(OTU_count = n(), .groups = 'drop')

count_depleted_ancom <- depleted_OTUs_ancom %>%
  group_by(study, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(OTU_count = n(), .groups = 'drop')

count_enriched_metagenomeseq <- enriched_OTUs_metagenomeseq %>%
  group_by(study, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(OTU_count = n(), .groups = 'drop')

count_depleted_metagenomeseq <- depleted_OTUs_metagenomeseq %>%
  group_by(study, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(OTU_count = n(), .groups = 'drop')


unique(count_enriched_ancom$study)
unique(count_enriched_metagenomeseq$study)
unique(count_depleted_ancom$study)
unique(count_depleted_metagenomeseq$study)

sort(table(count_enriched_ancom$Family), decreasing = TRUE)
sort(table(count_enriched_ancom$Genus), decreasing = TRUE)

count_depleted_ancom[is.na(count_depleted_ancom)] <- "unknown"
count_depleted_metagenomeseq[is.na(count_depleted_metagenomeseq)] <- "unknown"
count_enriched_ancom[is.na(count_enriched_ancom)] <- "unknown"
count_enriched_metagenomeseq[is.na(count_enriched_metagenomeseq)] <- "unknown"

#method <- c("ancom", "metagenomeseq")
#type <- c("enriched", "depleted")
# get(paste("count","enriched", "ancom", sep = "_"))

### family level
## enriched
shared_otus_enriched_ancom <- count_enriched_ancom
print(paste("Number of unique family:", length(unique(shared_otus_enriched_ancom$Family))))
df_enriched_ancom <- data.frame(OTUs = numeric(), Family = character(), studies = numeric(), stringsAsFactors = FALSE)
for (family in unique(shared_otus_enriched_ancom$Family)){
      num_studies <- length(unique(shared_otus_enriched_ancom$study[shared_otus_enriched_ancom$Family==family]))
      total_num_OTUs <- sum(shared_otus_enriched_ancom$OTU_count[shared_otus_enriched_ancom$Family==family])
      print(paste(total_num_OTUs, "OTUs belonging to", family, "found in", num_studies, "studies"))
      new_row <- data.frame(OTUs = total_num_OTUs, Family = family, studies = num_studies  )
      df_enriched_ancom <- rbind(df_enriched_ancom, new_row)
      }


shared_otus_enriched_metagenomeseq <- count_enriched_metagenomeseq
print(paste("Number of unique family:", length(unique(shared_otus_enriched_metagenomeseq$Family))))
df_enriched_metagenomeseq <- data.frame(OTUs = numeric(), Family = character(), studies = numeric(), stringsAsFactors = FALSE)
for (family in unique(shared_otus_enriched_metagenomeseq$Family)){
  num_studies <- length(unique(shared_otus_enriched_metagenomeseq$study[shared_otus_enriched_metagenomeseq$Family==family]))
  total_num_OTUs <- sum(shared_otus_enriched_metagenomeseq$OTU_count[shared_otus_enriched_metagenomeseq$Family==family])
  print(paste(total_num_OTUs, "OTUs belonging to", family, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Family = family, studies = num_studies  )
  df_enriched_metagenomeseq <- rbind(df_enriched_metagenomeseq, new_row)
}

## depleted
shared_otus_depleted_ancom <- count_depleted_ancom
print(paste("Number of unique family:", length(unique(shared_otus_depleted_ancom$Family))))
df_depleted_ancom <- data.frame(OTUs = numeric(), Family = character(), studies = numeric(), stringsAsFactors = FALSE)
for (family in unique(shared_otus_depleted_ancom$Family)){
  num_studies <- length(unique(shared_otus_depleted_ancom$study[shared_otus_depleted_ancom$Family==family]))
  total_num_OTUs <- sum(shared_otus_depleted_ancom$OTU_count[shared_otus_depleted_ancom$Family==family])
  print(paste(total_num_OTUs, "OTUs belonging to", family, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Family = family, studies = num_studies  )
  df_depleted_ancom <- rbind(df_depleted_ancom, new_row)
}


shared_otus_depleted_metagenomeseq <- count_depleted_metagenomeseq
print(paste("Number of unique family:", length(unique(shared_otus_depleted_metagenomeseq$Family))))
df_depleted_metagenomeseq <- data.frame(OTUs = numeric(), Family = character(), studies = numeric(), stringsAsFactors = FALSE)
for (family in unique(shared_otus_depleted_metagenomeseq$Family)){
  num_studies <- length(unique(shared_otus_depleted_metagenomeseq$study[shared_otus_depleted_metagenomeseq$Family==family]))
  total_num_OTUs <- sum(shared_otus_depleted_metagenomeseq$OTU_count[shared_otus_depleted_metagenomeseq$Family==family])
  print(paste(total_num_OTUs, "OTUs belonging to", family, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Family = family, studies = num_studies  )
  df_depleted_metagenomeseq <- rbind(df_depleted_metagenomeseq, new_row)
}


### genus level
## enriched
print(paste("Number of unique genus:", length(unique(shared_otus_enriched_ancom$Genus))))
df_enriched_ancom_genus <- data.frame(OTUs = numeric(), Genus = character(), studies = numeric(), stringsAsFactors = FALSE)
for (genus in unique(shared_otus_enriched_ancom$Genus)){
  num_studies <- length(unique(shared_otus_enriched_ancom$study[shared_otus_enriched_ancom$Genus==genus]))
  total_num_OTUs <- sum(shared_otus_enriched_ancom$OTU_count[shared_otus_enriched_ancom$Genus==genus])
  print(paste(total_num_OTUs, "OTUs belonging to", genus, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Genus = genus, studies = num_studies  )
  df_enriched_ancom_genus <- rbind(df_enriched_ancom_genus, new_row)
}

print(paste("Number of unique genus:", length(unique(shared_otus_enriched_metagenomeseq$Genus))))
df_enriched_metagenomeseq_genus <- data.frame(OTUs = numeric(), Genus = character(), studies = numeric(), stringsAsFactors = FALSE)
for (genus in unique(shared_otus_enriched_metagenomeseq$Genus)){
  num_studies <- length(unique(shared_otus_enriched_metagenomeseq$study[shared_otus_enriched_metagenomeseq$Genus==genus]))
  total_num_OTUs <- sum(shared_otus_enriched_metagenomeseq$OTU_count[shared_otus_enriched_metagenomeseq$Genus==genus])
  print(paste(total_num_OTUs, "OTUs belonging to", genus, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Genus = genus, studies = num_studies  )
  df_enriched_metagenomeseq_genus <- rbind(df_enriched_metagenomeseq_genus, new_row)
}

## depleted
print(paste("Number of unique genus:", length(unique(shared_otus_depleted_ancom$Genus))))
df_depleted_ancom_genus <- data.frame(OTUs = numeric(), Genus = character(), studies = numeric(), stringsAsFactors = FALSE)
for (genus in unique(shared_otus_depleted_ancom$Genus)){
  num_studies <- length(unique(shared_otus_depleted_ancom$study[shared_otus_depleted_ancom$Genus==genus]))
  total_num_OTUs <- sum(shared_otus_depleted_ancom$OTU_count[shared_otus_depleted_ancom$Genus==genus])
  print(paste(total_num_OTUs, "OTUs belonging to", genus, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Genus = genus, studies = num_studies  )
  df_depleted_ancom_genus <- rbind(df_depleted_ancom_genus, new_row)
}

print(paste("Number of unique genus:", length(unique(shared_otus_depleted_metagenomeseq$Genus))))
df_depleted_metagenomeseq_genus <- data.frame(OTUs = numeric(), Genus = character(), studies = numeric(), stringsAsFactors = FALSE)
for (genus in unique(shared_otus_depleted_metagenomeseq$Genus)){
  num_studies <- length(unique(shared_otus_depleted_metagenomeseq$study[shared_otus_depleted_metagenomeseq$Genus==genus]))
  total_num_OTUs <- sum(shared_otus_depleted_metagenomeseq$OTU_count[shared_otus_depleted_metagenomeseq$Genus==genus])
  print(paste(total_num_OTUs, "OTUs belonging to", genus, "found in", num_studies, "studies"))
  new_row <- data.frame(OTUs = total_num_OTUs, Genus = genus, studies = num_studies  )
  df_depleted_metagenomeseq_genus <- rbind(df_depleted_metagenomeseq_genus, new_row)
}


df_depleted_ancom$method <- "ancom"
df_depleted_ancom_genus$method <- "ancom"
df_enriched_ancom$method <- "ancom"
df_enriched_ancom_genus$method <- "ancom"

df_depleted_metagenomeseq$method <- "metagenomeseq"
df_depleted_metagenomeseq_genus$method <- "metagenomeseq"
df_enriched_metagenomeseq$method <- "metagenomeseq"
df_enriched_metagenomeseq_genus$method <- "metagenomeseq"

df_depleted_ancom$type <- "depleted"
df_depleted_ancom_genus$type <- "depleted"
df_enriched_ancom$type <- "enriched"
df_enriched_ancom_genus$type <- "enriched"

df_depleted_metagenomeseq$type <- "depleted"
df_depleted_metagenomeseq_genus$type <- "depleted"
df_enriched_metagenomeseq$type <- "enriched"
df_enriched_metagenomeseq_genus$type <- "enriched"


## merge the tables
### count per genus
df_ancom_genus <- merge(df_enriched_ancom_genus, df_depleted_ancom_genus, by = "Genus", all = TRUE,  suffixes = c(".enriched",".depleted"))
df_ancom_genus$OTUs.enriched[is.na(df_ancom_genus$OTUs.enriched)] <- 0
df_ancom_genus$studies.enriched[is.na(df_ancom_genus$studies.enriched)] <- 0
df_ancom_genus$OTUs.depleted[is.na(df_ancom_genus$OTUs.depleted)] <- 0
df_ancom_genus$studies.depleted[is.na(df_ancom_genus$studies.depleted)] <- 0
df_ancom_genus$total_OTUs <- df_ancom_genus$OTUs.enriched + df_ancom_genus$OTUs.depleted
df_ancom_genus <- df_ancom_genus[,c(1,2,3,6,7,10)]

df_metagenomeseq_genus <- merge(df_enriched_metagenomeseq_genus, df_depleted_metagenomeseq_genus, by = "Genus", all = TRUE,  suffixes = c(".enriched",".depleted"))
df_metagenomeseq_genus$OTUs.enriched[is.na(df_metagenomeseq_genus$OTUs.enriched)] <- 0
df_metagenomeseq_genus$studies.enriched[is.na(df_metagenomeseq_genus$studies.enriched)] <- 0
df_metagenomeseq_genus$OTUs.depleted[is.na(df_metagenomeseq_genus$OTUs.depleted)] <- 0
df_metagenomeseq_genus$studies.depleted[is.na(df_metagenomeseq_genus$studies.depleted)] <- 0
df_metagenomeseq_genus$total_OTUs <- df_metagenomeseq_genus$OTUs.enriched + df_metagenomeseq_genus$OTUs.depleted
df_metagenomeseq_genus <- df_metagenomeseq_genus[,c(1,2,3,6,7,10)]

df_genus <- merge(df_ancom_genus, df_metagenomeseq_genus, by = "Genus", all = TRUE,  suffixes = c(".ancom",".metagenomeseq"))
df_genus[is.na(df_genus)] <- 0
write.csv(df_genus, "count_shared_OTUs_genus_new.csv", row.names = FALSE)

### count per family
df_ancom_family <- merge(df_enriched_ancom, df_depleted_ancom, by = "Family", all = TRUE,  suffixes = c(".enriched",".depleted"))
df_ancom_family$OTUs.enriched[is.na(df_ancom_family$OTUs.enriched)] <- 0
df_ancom_family$studies.enriched[is.na(df_ancom_family$studies.enriched)] <- 0
df_ancom_family$OTUs.depleted[is.na(df_ancom_family$OTUs.depleted)] <- 0
df_ancom_family$studies.depleted[is.na(df_ancom_family$studies.depleted)] <- 0
df_ancom_family$total_OTUs <- df_ancom_family$OTUs.enriched + df_ancom_family$OTUs.depleted
df_ancom_family <- df_ancom_family[,c(1,2,3,6,7,10)]

df_metagenomeseq_family <- merge(df_enriched_metagenomeseq, df_depleted_metagenomeseq, by = "Family", all = TRUE,  suffixes = c(".enriched",".depleted"))
df_metagenomeseq_family$OTUs.enriched[is.na(df_metagenomeseq_family$OTUs.enriched)] <- 0
df_metagenomeseq_family$studies.enriched[is.na(df_metagenomeseq_family$studies.enriched)] <- 0
df_metagenomeseq_family$OTUs.depleted[is.na(df_metagenomeseq_family$OTUs.depleted)] <- 0
df_metagenomeseq_family$studies.depleted[is.na(df_metagenomeseq_family$studies.depleted)] <- 0
df_metagenomeseq_family$total_OTUs <- df_metagenomeseq_family$OTUs.enriched + df_metagenomeseq_family$OTUs.depleted
df_metagenomeseq_family <- df_metagenomeseq_family[,c(1,2,3,6,7,10)]

df_family <- merge(df_ancom_family, df_metagenomeseq_family, by = "Family", all = TRUE,  suffixes = c(".ancom",".metagenomeseq"))
df_family[is.na(df_family)] <- 0
write.csv(df_family, "count_shared_OTUs_family_new.csv", row.names = FALSE)
# family_genus_enriched <- merge(enriched_OTUs_ancom[, c("Genus", "Family")], enriched_OTUs_metagenomeseq[, c("Genus", "Family")], by = "Genus", all = TRUE)
# family_genus_enriched <- family_genus_enriched %>% mutate(Family = coalesce(Family.x, Family.y)) %>% distinct()
# family_genus_depleted <- merge(depleted_OTUs_ancom[, c("Genus", "Family")], depleted_OTUs_metagenomeseq[, c("Genus", "Family")], by = "Genus", all = TRUE)
# family_genus_depleted <- family_genus_depleted %>% mutate(Family = coalesce(Family.x, Family.y)) %>% distinct()
# family_genus <- merge(family_genus_enriched[, c("Genus", "Family")], family_genus_depleted[, c("Genus", "Family")], by = "Genus", all = TRUE)
# family_genus <- family_genus %>% mutate(Family = coalesce(Family.x, Family.y)) %>% distinct()
# df_genus <- merge(df_genus, family_genus, by = "Genus", all.x = TRUE)
 

## check genera overlap
enriched_genera_ancom <- read.csv("../shared_genera/enriched_genera_overlap_ANCOM.csv")
enriched_genera_metagenomeseq <- read.csv("../shared_genera/enriched_genera_overlap_CSS.csv")
depleted_genera_ancom <- read.csv("../shared_genera/depleted_genera_overlap_ANCOM (1).csv")
depleted_genera_metagenomeseq <- read.csv("../shared_genera/depleted_genera_overlap_CSS (1).csv")
enriched_genera_ancom$Genus <- enriched_genera_ancom$taxon
depleted_genera_ancom$Genus <- depleted_genera_ancom$taxon
genera_ancom <- merge(enriched_genera_ancom[, c("Genus", "n")], depleted_genera_ancom[, c("Genus", "n")], by = "Genus", all = TRUE, suffixes = c(".enriched", ".depleted"))
genera_metagenomeseq<- merge(enriched_genera_metagenomeseq[, c("Genus", "n")], depleted_genera_metagenomeseq[, c("Genus", "n")], by = "Genus", all = TRUE, suffixes = c(".enriched", ".depleted"))


## compare ancom and metagenomeseq results
enriched_OTUs <- list(ANCOM = enriched_OTUs_ancom$taxon,
                      metagenomeSeq = enriched_OTUs_metagenomeseq$Taxa)

enriched_OTUs_gg_venn <- ggVennDiagram(enriched_OTUs) +
  theme_minimal() + 
  labs(title = "Enriched OTUs") + coord_flip() +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")


depleted_OTUs <- list(ANCOM = depleted_OTUs_ancom$taxon,
                      metagenomeSeq = depleted_OTUs_metagenomeseq$Taxa)

depleted_OTUs_gg_venn <- ggVennDiagram(depleted_OTUs) +
  theme_minimal() + 
  labs(title = "Depleted OTUs") + coord_flip() +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")


enriched_OTUs_genus <- list(ANCOM = enriched_OTUs_ancom$Genus,
                            metagenomeSeq = enriched_OTUs_metagenomeseq$Genus)

enriched_OTUs_genus_gg_venn <- ggVennDiagram(enriched_OTUs_genus) +
  theme_minimal() + 
  labs(title = "Enriched OTUs (genus)") + coord_flip() +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")

depleted_OTUs_genus <- list(ANCOM = depleted_OTUs_ancom$Genus,
                            metagenomeSeq = depleted_OTUs_metagenomeseq$Genus)

depleted_OTUs_genus_gg_venn <- ggVennDiagram(depleted_OTUs_genus) +
  theme_minimal() + 
  labs(title = "Depleted OTUs (genus)") + coord_flip() +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")


OTUs_list <- list(ANCOM_enriched = enriched_OTUs_ancom$taxon, 
                  metagenomeSeq_enriched = enriched_OTUs_metagenomeseq$Taxa,
                  metagenomeSeq_depleted = depleted_OTUs_metagenomeseq$Taxa,
                  ANCOM_depleted = depleted_OTUs_ancom$taxon)
OTUs_gg_venn <- ggVennDiagram(OTUs_list,set_color = c("blue","black","red","yellow")) +
  theme_minimal() + 
  labs(title = "OTUs") +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")

OTUs_genus_list <- list(ANCOM_enriched = enriched_OTUs_ancom$Genus, 
                  metagenomeSeq_enriched = enriched_OTUs_metagenomeseq$Genus,
                  metagenomeSeq_depleted = depleted_OTUs_metagenomeseq$Genus,
                  ANCOM_depleted = depleted_OTUs_ancom$Genus)
OTUs_genus_gg_venn <- ggVennDiagram(OTUs_genus_list,label = "count",label_alpha = 0) +
  theme_minimal() + 
  labs(title = "OTUs") +
  theme(
    axis.title.x = element_blank(),  # Remove X axis title
    axis.title.y = element_blank(),  # Remove Y axis title
    axis.ticks = element_blank(),     # Remove axis ticks
    axis.text.x = element_blank(),    # Remove X axis text
    axis.text.y = element_blank(),      # Remove Y axis text
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank()   # Remove background panel
  ) +
  scale_fill_distiller(palette = "RdBu")




df_combined_ancom <- rbind(df_enriched_ancom, df_depleted_ancom)
df_transformed_ancom <- df_combined_ancom %>%
  mutate(OTUs = ifelse(type == "depleted", -OTUs, OTUs)) %>%
  arrange(desc(OTUs)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))


ggplot(df_transformed_ancom, aes(x = Family, y = OTUs, fill = as.integer(studies))) +
  geom_bar(stat = "identity") +
  labs(title = "ANCOM: OTU Counts by Family (Enriched vs. Depleted)",
       x = "Family",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies",
                      labels = c("2", "4", "6", "8", "10"),
                      breaks = c(2, 4, 6, 8, 10)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  


df_combined_metagenomeseq <- rbind(df_enriched_metagenomeseq, df_depleted_metagenomeseq)
df_transformed_metagenomeseq <- df_combined_metagenomeseq %>%
  mutate(OTUs = ifelse(type == "depleted", -OTUs, OTUs)) %>%
  arrange(desc(OTUs)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))


ggplot(df_transformed_metagenomeseq, aes(x = Family, y = OTUs, fill = as.integer(studies))) +
  geom_bar(stat = "identity") +
  labs(title = "metagenomeSeq: OTU Counts by Family (Enriched vs. Depleted)",
       x = "Family",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies",
                      labels = c("2", "4", "6", "8", "10"),
                      breaks = c(2, 4, 6, 8, 10)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  


df_combined_ancom_genus <- rbind(df_enriched_ancom_genus, df_depleted_ancom_genus)
df_transformed_ancom_genus <- df_combined_ancom_genus %>%
  mutate(OTUs = ifelse(type == "depleted", -OTUs, OTUs)) %>%
  arrange(desc(OTUs)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))


ggplot(df_transformed_ancom_genus, aes(x = Genus, y = OTUs, fill = as.integer(studies))) +
  geom_bar(stat = "identity") +
  labs(title = "ANCOM: OTU Counts by Genus (Enriched vs. Depleted)",
       x = "Genus",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies",
                      labels = c("2", "4", "6", "8", "10"),
                      breaks = c(2, 4, 6, 8, 10)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  


df_combined_metagenomeseq_genus <- rbind(df_enriched_metagenomeseq_genus, df_depleted_metagenomeseq_genus)
df_transformed_metagenomeseq_genus <- df_combined_metagenomeseq_genus %>%
  mutate(OTUs = ifelse(type == "depleted", -OTUs, OTUs)) %>%
  arrange(desc(OTUs)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))


ggplot(df_transformed_metagenomeseq_genus, aes(x = Genus, y = OTUs, fill = as.integer(studies))) +
  geom_bar(stat = "identity") +
  labs(title = "metagenomeSeq: OTU Counts by Genus (Enriched vs. Depleted)",
       x = "Genus",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies",
                      labels = c("2", "4", "6", "8", "10"),
                      breaks = c(2, 4, 6, 8, 10)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  

# Select top 20 genus by #enriched OTUs 
df_top_genera_ancom <- df_transformed_ancom_genus[seq(1,20),]

# Plot only top families
ggplot(df_top_genera_ancom, aes(x = Genus, y = OTUs, fill = studies)) +
  geom_bar(stat = "identity") +
  labs(title = "ANCOM: Top 20 Genera by OTU Counts (Enriched)",
       x = "Genus",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")

# Select top 20 genus by #enriched OTUs 
df_top_genera_metagenomeseq <- df_transformed_metagenomeseq_genus [seq(1,20),]

# Plot only top families
ggplot(df_top_genera_metagenomeseq, aes(x = Genus, y = OTUs, fill = studies)) +
  geom_bar(stat = "identity") +
  labs(title = "metagenomeSeq: Top 20 Genera by OTU Counts (Enriched)",
       x = "Genus",
       y = "Number of OTUs") +
  theme_pubclean() + theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Number of Studies") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")

# intersect(df_enriched_ancom$Family, RF_family$Family)
# intersect(df_depleted_ancom$Family, RF_family$Family)
# 
# intersect(intersect(df_enriched_ancom$Family, RF_family$Family),
#           intersect(df_depleted_ancom$Family, RF_family$Family))
# 
# intersect(df_enriched_ancom$Family, RF_family$Family)[!
#   (intersect(df_enriched_ancom$Family, RF_family$Family) %in% intersect(df_depleted_ancom$Family, RF_family$Family))]

RF_family <- read.csv("../RF/RF_Output_family_renamed.csv")

table(RF_family$Phylum) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
RF_family$Phylum_plot <- "Other"
RF_family$Phylum_plot[RF_family$Phylum=="Proteobacteria"] <- "Proteobacteria"
RF_family$Phylum_plot[RF_family$Phylum=="Actinobacteriota"] <- "Actinobacteriota"
RF_family$Phylum_plot[RF_family$Phylum=="Acidobacteriota"] <- "Acidobacteriota"
RF_family$Phylum_plot[RF_family$Phylum=="Bacteroidota"] <- "Bacteroidota"
RF_family$Phylum_plot[RF_family$Phylum=="Firmicutes"] <- "Firmicutes"

df_RF_family <- RF_family %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#999999",  "#D55E00","#CC79A7", "#0072B2")
ggplot(df_RF_family[seq(1:20),], aes(x = Family, y = MeanDecreaseAccuracy, fill = Phylum_plot)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 most important features",
       x = "Family",
       y = "Mean Decrease in Accuracy") +
  theme_pubclean() + theme_pubr() +
  scale_fill_manual(values = colors, name = "Phylum") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  )


RF_genus <- read.csv("../RF/RF_Output_genus_renamed.csv")

table(RF_genus$Phylum) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
RF_genus$Phylum_plot <- "Other"
RF_genus$Phylum_plot[RF_genus$Phylum=="Proteobacteria"] <- "Proteobacteria"
RF_genus$Phylum_plot[RF_genus$Phylum=="Actinobacteriota"] <- "Actinobacteriota"
RF_genus$Phylum_plot[RF_genus$Phylum=="Acidobacteriota"] <- "Acidobacteriota"
RF_genus$Phylum_plot[RF_genus$Phylum=="Bacteroidota"] <- "Bacteroidota"
RF_genus$Phylum_plot[RF_genus$Phylum=="Firmicutes"] <- "Firmicutes"

df_RF_genus <- RF_genus %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#999999",  "#D55E00","#CC79A7", "#0072B2")
ggplot(df_RF_genus[seq(1:20),], aes(x = Genus, y = MeanDecreaseAccuracy, fill = Phylum_plot)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 most important features",
       x = "Genus",
       y = "Mean Decrease in Accuracy") +
  theme_pubclean() + theme_pubr() +
  scale_fill_manual(values = colors, name = "Phylum") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

RF_OTU <- read.csv("../RF/RF_Output_OTU_Without_SHAP_filter_0001_25.csv")

table(RF_OTU$Phylum) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
RF_OTU$Phylum_plot <- "Other"
RF_OTU$Phylum_plot[RF_OTU$Phylum=="Proteobacteria"] <- "Proteobacteria"
RF_OTU$Phylum_plot[RF_OTU$Phylum=="Actinobacteriota"] <- "Actinobacteriota"
RF_OTU$Phylum_plot[RF_OTU$Phylum=="Acidobacteriota"] <- "Acidobacteriota"
RF_OTU$Phylum_plot[RF_OTU$Phylum=="Bacteroidota"] <- "Bacteroidota"
RF_OTU$Phylum_plot[RF_OTU$Phylum=="Firmicutes"] <- "Firmicutes"

df_RF_OTU <- RF_OTU %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(taxa = factor(taxa, levels = unique(taxa)))

colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#999999",  "#D55E00","#CC79A7", "#0072B2")
ggplot(df_RF_OTU[seq(1,20),], aes(x = taxa, y = MeanDecreaseAccuracy, fill = Phylum_plot)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 most important features",
       x = "OTU",
       y = "Mean Decrease in Accuracy") +
  theme_pubclean() + theme_pubr() +
  scale_fill_manual(values = colors, name = "Phylum") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),  # Rotate and size adjustments
    plot.title = element_text(hjust = 0.5)  # Center the title
  )


table(df_RF_OTU$Genus[seq(1,50)])
