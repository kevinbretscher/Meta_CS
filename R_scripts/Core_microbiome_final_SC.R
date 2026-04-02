library(qiime2R)
library(phyloseq)
library(metagenomeSeq)
library(tidyverse)
library(microbiome)
#library(microViz)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggpubr)
#library(ggside)
library("doParallel")
library(speedyseq)



CSS_norm <- function(physeq1){
  metseq <- phyloseq_to_metagenomeSeq(physeq1)
  #calculate the proper percentile by which to normalize counts
  p <- cumNormStatFast(metseq)
  #calculate the scaling factors 
  metseq <- cumNorm(metseq, p = p)
  #create dataframe with normalized counts
  metseq.counts <- MRcounts(metseq, norm = TRUE, log = FALSE)
  #create phyloseq object with transformed data
  otu_table(physeq1) <- otu_table(metseq.counts, taxa_are_rows = TRUE)
  return(physeq1)
}

Rename_taxa <- function(physeq,taxon) {
  physeq2 <- physeq
  if (taxon == "Genus") {
    taxa_names(physeq2) <- paste(tax_table(physeq2)[,"Phylum"],
                                 tax_table(physeq2)[,"Class"],
                                 tax_table(physeq2)[,"Order"],
                                 tax_table(physeq2)[,"Family"],
                                 tax_table(physeq2)[,"Genus"],sep="+")
  } else if (taxon == "Family") {
    taxa_names(physeq2) <- paste(tax_table(physeq2)[,"Phylum"],
                                 tax_table(physeq2)[,"Class"],
                                 tax_table(physeq2)[,"Order"],
                                 tax_table(physeq2)[,"Family"],sep="+")
  } else if (taxon == "Order") {
    taxa_names(physeq2) <- paste(tax_table(physeq2)[,"Phylum"],
                                 tax_table(physeq2)[,"Class"],
                                 tax_table(physeq2)[,"Order"],sep="+")
  } else if (taxon == "Class") {
    taxa_names(physeq2) <- paste(tax_table(physeq2)[,"Phylum"],
                                 tax_table(physeq2)[,"Class"],sep="+")
  } else if (taxon == "Phylum") {
    taxa_names(physeq2) <- paste(tax_table(physeq2)[,"Phylum"])
  }
  return(physeq2)
}


ps.import <- readRDS("phyloseq.RDS")

physeq.Rhizo <- subset_samples(ps.import, B_R=="R")

#### Family


ps.family <- speedyseq::tax_glom(physeq.Rhizo, taxrank="Family")
ps.family <- CSS_norm(ps.family)
ps.family <- Rename_taxa(ps.family,"Family")

ps.family <- transform_sample_counts(ps.family, function(x) x / sum(x) )

df.family <- psmelt(ps.family)

# how many samples per category

unique_sample_count.family <- df.family  %>%
  group_by(C_S) %>%
  summarize(unique_sample_count = n_distinct(Sample), .groups = 'drop')

df.family.reduced <- df.family %>% 
  filter(Abundance > 0) # drop useless rows to make processing faster and to not count absent taxa

#group by sample and family and count rows. Divide by unqiue samples

df.family.core <- df.family.reduced %>%
  group_by(OTU, C_S) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    PR_samples = ifelse(mean(Abundance, na.rm = TRUE) > 0, n(), 0),
    PR_studies = ifelse(mean(Abundance, na.rm = TRUE) > 0, n_distinct(Study), 0),
    .groups = 'drop'
  ) %>%
  left_join(unique_sample_count.family, by = "C_S") %>%
  mutate(PR_samples_rel = PR_samples / unique_sample_count) %>%
  separate(OTU, into = c("Phylum", "Class", "Order", "Family"), sep = "\\+", remove = FALSE)

 C.plot.family <- ggscatter(df.family.core[df.family.core$C_S == "C",], x = "log10mean_abundance", y = "PR_samples_rel", color = "PR_studies", alpha = 0.5, title = " Con") + 
     scale_colour_gradient(low = "dimgray", high = "#FEE08B", na.value = NA, name = "Conducive\n#studies", breaks = c(0, 4, 8, 12)) +
     scale_y_continuous(limits = c(0, 1)) + 
     #scale_x_continuous(trans='log10') +
     geom_hline(yintercept=c(0.9), linetype="dashed") +
     geom_vline(xintercept=c(-3), linetype="dashed")  + font("xy.text", size = 6) + font("xylab", size = 8, face = "bold") + theme(plot.title = element_text(hjust = 0.5)) +
     #geom_text_repel(aes(label = ifelse(PR_samples_rel ==1 & mean_abundance > 0.01 ,as.character(Family),"")),max.overlaps = Inf, force = 10)+
     labs(x="Log 10(mean relative abundance)", y="Prevalence", color = "Conducive\n#studies", title = "Family")+
     theme_pubr(base_size = 18)
 S.plot.family <- ggscatter(df.family.core[df.family.core$C_S == "S",], x = "log10mean_abundance", y = "PR_samples_rel", color = "PR_studies", alpha = 0.5, title = "Sup") + 
     scale_colour_gradient(low = "dimgray", high = "#91BFDB", na.value = NA, name = "Suppressive\n#studies", breaks = c(0, 4, 8, 12)) +
     scale_y_continuous(limits = c(0, 1)) + 
     #scale_x_continuous(trans='log10') +
     geom_hline(yintercept=c(0.9), linetype="dashed") +
     geom_vline(xintercept=c(-3), linetype="dashed")  + font("xy.text", size = 6) + font("xylab", size = 8, face = "bold") + theme(plot.title = element_text(hjust = 0.5))+
     #geom_text_repel(aes(label = ifelse(PR_samples_rel==1 & mean_abundance > 0.01 ,as.character(Family),"")),max.overlaps = Inf, force = 10) +
     labs(x="Log 10(mean relative abundance)", y="Prevalence", color = "# Studies", title = "Family") +
     theme_pubr(base_size = 18)


ps.Genus <- speedyseq::tax_glom(physeq.Rhizo, taxrank="Genus")
ps.Genus <- CSS_norm(ps.Genus)
ps.Genus <- Rename_taxa(ps.Genus,"Genus")

ps.Genus <- transform_sample_counts(ps.Genus, function(x) x / sum(x) )

df.Genus <- psmelt(ps.Genus)

# how many samples per category

unique_sample_count.Genus <- df.Genus  %>%
  group_by(C_S) %>%
  summarize(unique_sample_count = n_distinct(Sample), .groups = 'drop')

df.Genus.reduced <- df.Genus %>% 
  filter(Abundance > 0) # drop useless rows to make processing faster and to not count absent taxa

#group by sample and Genus and count rows. Divide by unqiue samples

df.Genus.core <- df.Genus.reduced %>%
  group_by(OTU, C_S) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    PR_samples = ifelse(mean(Abundance, na.rm = TRUE) > 0, n(), 0),
    PR_studies = ifelse(mean(Abundance, na.rm = TRUE) > 0, n_distinct(Study), 0),
    .groups = 'drop'
  ) %>%
  left_join(unique_sample_count.Genus, by = "C_S") %>%
  mutate(PR_samples_rel = PR_samples / unique_sample_count) %>%
  separate(OTU, into = c("Phylum", "Class", "Order", "Family","Genus"), sep = "\\+", remove = FALSE)

C.plot.Genus <- ggscatter(df.Genus.core[df.Genus.core$C_S == "C",], x = "log10mean_abundance", y = "PR_samples_rel", color = "PR_studies", alpha = 0.5, title = " Con") + 
  scale_colour_gradient(low = "dimgray", high = "#FEE08B", na.value = NA, name = "Conducive\n#studies", breaks = c(0, 4, 8, 12)) +
  scale_y_continuous(limits = c(0, 1)) + 
  #scale_x_continuous(trans='log10') +
  geom_hline(yintercept=c(0.9), linetype="dashed") +
  geom_vline(xintercept=c(-3), linetype="dashed")  + font("xy.text", size = 6) + font("xylab", size = 8, face = "bold") + theme(plot.title = element_text(hjust = 0.5)) +
  #geom_text_repel(aes(label = ifelse(PR_samples_rel ==1 & mean_abundance > 0.01 ,as.character(Genus),"")),max.overlaps = Inf, force = 10)+
  labs(x="Log 10(mean relative abundance)", y="Prevalence", color = "Conducive\n#studies", title = "Genus")+
  theme_pubr(base_size = 18)
S.plot.Genus <- ggscatter(df.Genus.core[df.Genus.core$C_S == "S",], x = "log10mean_abundance", y = "PR_samples_rel", color = "PR_studies", alpha = 0.5, title = "Sup") + 
  scale_colour_gradient(low = "dimgray", high = "#91BFDB", na.value = NA, name = "Suppressive\n#studies", breaks = c(0, 4, 8, 12)) +
  scale_y_continuous(limits = c(0, 1)) + 
  #scale_x_continuous(trans='log10') +
  geom_hline(yintercept=c(0.9), linetype="dashed") +
  geom_vline(xintercept=c(-3), linetype="dashed")  + font("xy.text", size = 6) + font("xylab", size = 8, face = "bold") + theme(plot.title = element_text(hjust = 0.5))+
  #geom_text_repel(aes(label = ifelse(PR_samples_rel==1 & mean_abundance > 0.01 ,as.character(Genus),"")),max.overlaps = Inf, force = 10) +
  labs(x="Log 10(mean relative abundance)", y="Prevalence", color = "Suppressive\n#studies", title = "Genus") +
  theme_pubr(base_size = 18)
ggarrange(C.plot.family, S.plot.family, C.plot.Genus, S.plot.Genus,
          ncol = 2, nrow = 2,
          common.legend = F) +
  theme(rect = element_rect(fill = "transparent"))

# plot % of total microbiome

C.core.family <- df.family.core %>%
  filter(C_S == "C") %>%
  filter(mean_abundance > 0.001) %>%
  filter(PR_samples_rel > 0.9)

S.core.family <- df.family.core %>%
  filter(C_S == "S") %>%
  filter(mean_abundance > 0.001) %>%
  filter(PR_samples_rel > 0.9) 

C.core.Genus <- df.Genus.core %>%
  filter(C_S == "C") %>%
  filter(mean_abundance > 0.001) %>%
  filter(PR_samples_rel > 0.9)

S.core.Genus <- df.Genus.core %>%
  filter(C_S == "S") %>%
  filter(mean_abundance > 0.001) %>%
  filter(PR_samples_rel > 0.9) 


df.core.acc.C.family <- df.family %>%
  filter(C_S == "C") %>%
  mutate(core.taxa = case_when(OTU %in% C.core.family$OTU ~ "Core", TRUE ~ "Accessory" )) %>%
  group_by(Sample,core.taxa) %>%
  summarise(Abundance = sum(Abundance))

df.core.acc.S.family <- df.family %>%
  filter(C_S == "S") %>%
  mutate(core.taxa = case_when(OTU %in% S.core.family$OTU ~ "Core", TRUE ~ "Accessory" )) %>%
  group_by(Sample,core.taxa) %>%
  summarise(Abundance = sum(Abundance))

df.core.acc.C.Genus <- df.Genus %>%
  filter(C_S == "C") %>%
  mutate(core.taxa = case_when(OTU %in% C.core.Genus$OTU ~ "Core", TRUE ~ "Accessory" )) %>%
  group_by(Sample,core.taxa) %>%
  summarise(Abundance = sum(Abundance))

df.core.acc.S.Genus <- df.Genus %>%
  filter(C_S == "S") %>%
  mutate(core.taxa = case_when(OTU %in% S.core.Genus$OTU ~ "Core", TRUE ~ "Accessory" )) %>%
  group_by(Sample,core.taxa) %>%
  summarise(Abundance = sum(Abundance))


df.core.acc.C.family <- ggboxplot(df.core.acc.C.family, "core.taxa", "Abundance", ylim = c(0,1), size = 1, color = "#E69F00", fill = "#E69F00", alpha = 0.6, title = "Family con") + theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")
df.core.acc.S.family <- ggboxplot(df.core.acc.S.family, "core.taxa", "Abundance", ylim = c(0,1), size = 1, color = "#56B4E9", fill = "#56B4E9", alpha = 0.6, title = "Family sup") + theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")
df.core.acc.C.Genus <- ggboxplot(df.core.acc.C.Genus, "core.taxa", "Abundance", ylim = c(0,1), size = 1, color = "#009E73", fill = "#009E73", alpha = 0.6, title = "Genus con") + theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")
df.core.acc.S.Genus <- ggboxplot(df.core.acc.S.Genus, "core.taxa", "Abundance", ylim = c(0,1), size = 1, color = "#F0E442", fill = "#F0E442", alpha = 0.6, title = "Genus sup") + theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")


ggarrange(df.core.acc.C.family, df.core.acc.S.family, df.core.acc.C.Genus, 
          df.core.acc.S.Genus,
          ncol = 2, nrow = 2) + 
  theme(rect = element_rect(fill = "transparent")) 



##### VENN

venn.list <- list(Conducive_family = C.core.family$OTU, suppressive_family = S.core.family$OTU,
                  Conducive_Genus = C.core.Genus$OTU, suppressive_Genus = S.core.Genus$OTU)

venn.list_genus <- list(Conducive = C.core.Genus$OTU, Suppressive = S.core.Genus$OTU)
venn.list_family <- list(Conducive = C.core.family$OTU, Suppressive = S.core.family$OTU)

# library(ggVennDiagram)
# ggVennDiagram(venn.list[c("Conducive_family","suppressive_family")], label = "count") + scale_fill_distiller(palette = "RdBu")+ coord_flip()
# ggVennDiagram(venn.list[c("Conducive_Genus","suppressive_Genus")], label = "count") + scale_fill_distiller(palette = "RdBu") + coord_flip()

install.packages("eulerr")
library("eulerr")

p_venn_core_genus <- plot(euler(venn.list_genus), 
     fills = list(fill = c("Conducive" = "#FEE08B",
                           "Suppressive" =  "#91BFDB"), 
                  alpha = 0.8),
     labels = list(col = "black", font = 2, fontsize = 10),
     quantities = list(col = "black", font = 1, fontsize = 10),
     main = "")
ggsave("../../map/venn_core_genus.pdf", plot = p_venn_core_genus, width = 4, height = 3, units = "in", dpi = 300)



p_venn_core_family <- plot(euler(venn.list_family), 
     fills = list(fill = c("Conducive" = "#FEE08B",
                           "Suppressive" =  "#91BFDB"), 
                  alpha = 0.8),
     labels = list(col = "black", font = 2, fontsize = 10),
     quantities = list(col = "black", font = 1, fontsize = 10),
     main = "")
ggsave("../../map/venn_core_family.pdf", plot = p_venn_core_family, width = 4, height = 3, units = "in", dpi = 300)

venn <- Venn(venn.list) 


#unique in suppressive
discern(venn, slice1 = "suppressive_family", slice2 = "Conducive_family")
discern(venn, slice1 = "Conducive_family", slice2 = "suppressive_family")

discern(venn, slice1 = "suppressive_Genus", slice2 = "Conducive_Genus")
discern(venn, slice1 = "Conducive_Genus", slice2 = "suppressive_Genus")

# overlap
overlap(venn, slice = c("Conducive_family","suppressive_family"))
overlap(venn, slice = c("Conducive_Genus","suppressive_Genus"))

# combine
union(venn.list_family$Conducive, venn.list_family$Suppressive)
union(venn.list_genus$Conducive, venn.list_genus$Suppressive)

# intersect(df.Genus.core$OTU[df.Genus.core$C_S=="S" & 
#                     df.Genus.core$mean_abundance > 0.001 & 
#                     df.Genus.core$PR_samples_rel > 0.9],
# df.Genus.core$OTU[df.Genus.core$C_S=="C" & 
#                     df.Genus.core$mean_abundance > 0.001 & 
#                     df.Genus.core$PR_samples_rel < 0.9])

table(df.family.core$Family[df.family.core$mean_abundance > 0.001 & df.family.core$PR_samples_rel > 0.9])
table(df.Genus.core$Genus[df.Genus.core$mean_abundance > 0.001 & df.Genus.core$PR_samples_rel > 0.9])

union(df.Genus.core$OTU[df.Genus.core$C_S=="S" & 
                                         df.Genus.core$mean_abundance > 0.001 & 
                                         df.Genus.core$PR_samples_rel > 0.9],
df.Genus.core$OTU[df.Genus.core$C_S=="C" & 
                    df.Genus.core$mean_abundance > 0.001 & 
                    df.Genus.core$PR_samples_rel > 0.9])

df.Genus.core_filtered <- df.Genus.core[df.Genus.core$mean_abundance>0.001 & df.Genus.core$PR_samples_rel > 0.9,]

core_family <- unique(df.family.core$Family[df.family.core$mean_abundance > 0.001 & df.family.core$PR_samples_rel > 0.9])
core_genus <- unique(df.Genus.core$Genus[df.Genus.core$mean_abundance > 0.001 & df.Genus.core$PR_samples_rel > 0.9])
core_family_S <- unique(df.family.core$Family[df.family.core$C_S=="S" & df.family.core$mean_abundance > 0.001 & df.family.core$PR_samples_rel > 0.9])
core_genus_S <- unique(df.Genus.core$Genus[df.Genus.core$C_S=="S" & df.Genus.core$mean_abundance > 0.001 & df.Genus.core$PR_samples_rel > 0.9])
core_family_C <- unique(df.family.core$Family[df.family.core$C_S=="C" & df.family.core$mean_abundance > 0.001 & df.family.core$PR_samples_rel > 0.9])
core_genus_C <- unique(df.Genus.core$Genus[df.Genus.core$C_S=="C" & df.Genus.core$mean_abundance > 0.001 & df.Genus.core$PR_samples_rel > 0.9])

uniq_core_genus_C <- setdiff(core_genus_C, core_genus_S)
uniq_core_genus_S <- setdiff(core_genus_S, core_genus_C)
uniq_core_genus <- c(uniq_core_genus_C, uniq_core_genus_S)

uniq_core_family_C <- setdiff(core_family_C, core_family_S)
uniq_core_family_S <- setdiff(core_family_S, core_family_C)
uniq_core_family <- c(uniq_core_family_C, uniq_core_family_S)

count_shared_family <- read.csv("../../count_OTUs/count_shared_OTUs_family_new.csv")
count_shared_genus <- read.csv("../../count_OTUs/count_shared_OTUs_genus_new.csv")

count_shared_family_core <- count_shared_family[count_shared_family$Family %in% core_family,]
count_shared_genus_core <- count_shared_genus[count_shared_genus$Genus %in% core_genus,]

count_shared_family_core_S <- count_shared_family[count_shared_family$Family %in% core_family_S,]
count_shared_genus_core_S <- count_shared_genus[count_shared_genus$Genus %in% core_genus_S,]

write.csv(df.family.core, "data_core_family_plot.csv")
write.csv(df.Genus.core, "data_core_genus_plot.csv")

write.csv(C.core.family,"C_core_family.csv")
write.csv(S.core.family,"S_core_family.csv")
write.csv(C.core.Genus,"C_core_genus.csv")
write.csv(S.core.Genus,"S_core_genus.csv")
