library(qiime2R)
library(phyloseq)
install.packages("phyloseq")

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("metagenomeSeq")

# install.packages("devtools")
# library("devtools")
# install_github("Bioconductor-mirror/metagenomeSeq")

library(metagenomeSeq)
library(tidyverse)
#library(microbiome)
#library(microViz)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggside)
library(ggsci)
library("doParallel")
# install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
library(speedyseq)

#source("Custom_functions.R")

ps.import <- readRDS("../Scripts/phyloseq.RDS")
ps.phylum <- speedyseq::tax_glom(ps.import, taxrank="Phylum")
ps.order <- speedyseq::tax_glom(ps.import, taxrank="Order")
ps.family <- speedyseq::tax_glom(ps.import, taxrank="Family")
ps.genus <- speedyseq::tax_glom(ps.import, taxrank="Genus")

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

sample_ps <- function(ps, FUN = sample, ...){
  ids <- sample_names(ps)
  sampled_ids <- FUN(ids, ...)
  ps <- prune_samples(sampled_ids, ps)
  return(ps)
}

rel_abundance <- function(physeq1){
  physeq1 <- transform_sample_counts(physeq1, function(x) x/sum(x))
  
}

physeq.Rhizo <- subset_samples(ps.import, B_R=="R")

physeq.Rhizo_order <- subset_samples(ps.order, B_R=="R")
physeq.Rhizo_family <- subset_samples(ps.family, B_R=="R")
physeq.Rhizo_genus <- subset_samples(ps.genus, B_R=="R")

physeq.Rhizo.norm <- CSS_norm(physeq.Rhizo)
physeq.Rhizo.norm_order <- CSS_norm(physeq.Rhizo_order)
physeq.Rhizo.norm_family <- CSS_norm(physeq.Rhizo_family)
physeq.Rhizo.norm_genus <- CSS_norm(physeq.Rhizo_genus)

sample_data(physeq.Rhizo.norm)$Study <- as.character(sample_data(physeq.Rhizo.norm)$Study)
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="Adam"] <- "Ossowicki 2020, Wheat"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="B2R-Lilian"] <- "Abreu 2024, Wheat"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="Kwak"] <- "Cha 2015, Strawberry"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJCA010283"] <- "Wen 2023, Cucumber"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA432446"] <- "Zhou 2019, Banana"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA578725"] <- "Yin 2021, Wheat"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA603213"] <- "Ding 2021, Tobacco"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA627608"] <- "Shen 2022, Banana"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA649668"] <- "Goh 2020, Oil palm"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA921884"] <- "Yang 2023, Cucumber"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="PRJNA989386"] <- "Luo 2024, Peanut"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="Ruth-zwaagdijk"] <- "Expósito 2017, Sugar beet"
sample_data(physeq.Rhizo.norm)$Study[sample_data(physeq.Rhizo.norm)$Study=="Xiaogang"] <- "Li 2015, Cotton"
sample_data(physeq.Rhizo.norm)$Study <- factor(sample_data(physeq.Rhizo.norm)$Study, levels = c(
  "Shen 2022, Banana","Zhou 2019, Banana",
  "Li 2015, Cotton", 
  "Wen 2023, Cucumber", "Yang 2023, Cucumber",
  "Goh 2020, Oil palm",
  "Luo 2024, Peanut",
  "Cha 2015, Strawberry",
  "Expósito 2017, Sugar beet",
  "Ding 2021, Tobacco",
  "Abreu 2024, Wheat",  "Ossowicki 2020, Wheat",  "Yin 2021, Wheat"
))

table(sample_data(physeq.Rhizo.norm)$Region)
table(sample_data(physeq.Rhizo.norm)$Pathogen_Genus)
sample_data(physeq.Rhizo.norm)$Pathogen_Genus[sample_data(physeq.Rhizo.norm)$Pathogen_Genus=="Rizoctonia Solani AG2-2IIIB"] <- "Rhizoctonia"
sample_data(physeq.Rhizo.norm)$Region[sample_data(physeq.Rhizo.norm)$Region=="v3-v4"] <- "V3-V4"

r_OTU_CS <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ C_S,
                    data = as(sample_data(physeq.Rhizo.norm), "data.frame"))
r_OTU_Study <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ Study,
                       data = as(sample_data(physeq.Rhizo.norm), "data.frame"))
r_OTU_crop <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ Crop,
                       data = as(sample_data(physeq.Rhizo.norm), "data.frame"))
r_OTU_pathogen <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ Pathogen_Genus,
                       data = as(sample_data(physeq.Rhizo.norm), "data.frame"))
r_OTU_region <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ Region,
                        data = as(sample_data(physeq.Rhizo.norm), "data.frame"))

print(Suppressiveness = r_OTU_CS$R2[1],
           Study = r_OTU_Study$R2[1],
           Crop = r_OTU_crop$R2[1],
           Pathogen = r_OTU_pathogen$R2[1])
# Create the dataframe
df_R2 <- data.frame(
  Category = c("Soil", "Study", "Crop", "Pathogen", "Primer"),
  Value = c(2.083352, 54.67261, 35.01257, 22.72436, 15.905)
)

library(ggplot2)
library(ggpubr)

# Create the bar plot
ggplot(df_R2, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "PERMANOVA R2 (%)", title = "") +
  scale_fill_brewer(palette = "Set3") +
  theme_pubr(border = T) +
  theme(legend.position = "none") # Remove legend since categories are on the x-axis
  


#Cap analysis

GP.cap <- ordinate(physeq.Rhizo.norm, "CAP", "bray", ~C_S)
GP.cap_study <- ordinate(physeq.Rhizo.norm, "CAP", "bray", ~Study)

GP.cap_genus <- ordinate(physeq.Rhizo.norm_genus, "CAP", "bray", ~C_S)
GP.cap_family <- ordinate(physeq.Rhizo.norm_family, "CAP", "bray", ~C_S)
GP.cap_order <- ordinate(physeq.Rhizo.norm_order, "CAP", "bray", ~C_S)

p <- plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color="C_S", title="CAP ~C_S") +
  scale_color_aaas() +
  scale_fill_aaas() +
  geom_point(size = 3) + 
  geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr()

library(ggsci)

# Manually combine colors from multiple palettes
palette_extended <- c(pal_npg("nrc")(10), pal_jco("default")(10))

p_cap_CS <- plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color="C_S", title="") +
  geom_point(size = 2, alpha = 0.65) + scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"), name = "", labels =  c( "C"="Conducive", "S" = "Suppressive ")) +
  theme_pubr(border=T) + theme(legend.position = "none")
ggsave("../../map/cap_CS.pdf", plot = p_cap_CS, width = 3, height = 3, units = "in", dpi = 300)

p_cap_study <- plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color="Study", title="") +
  geom_point(size = 2, alpha = 0.5) + scale_color_manual(values =palette_extended ) +
  theme_pubr(border=T) + theme(legend.position = "none", legend.text = element_text(size = 10))
ggsave("../../map/cap_study.pdf", plot = p_cap_study, width = 3, height = 3, units = "in", dpi = 300)


print(p)
print(p2)

plot_ordination(physeq.Rhizo.norm, GP.cap_study, type="samples", color="Study", shape = "C_S") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color = "C_S", title="CAP ~C_S") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  geom_point(size = 2) + 
  # geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  # geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr(border=T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color = "Study", title="CAP ~C_S") +
  #scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  geom_point(size = 2) + 
  # geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  # geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr(border=T) + theme(legend.position = "none")




plot_ordination(physeq.Rhizo.norm_genus, GP.cap_genus, type="samples", color = "C_S", title="CAP ~C_S") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  geom_point(size = 3) + 
  geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr()

plot_ordination(physeq.Rhizo.norm_family, GP.cap_family, type="samples", color = "C_S", title="CAP ~C_S") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  geom_point(size = 3) + 
  geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr()

plot_ordination(physeq.Rhizo.norm_order, GP.cap_order, type="samples", color = "C_S", title="CAP ~C_S") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  geom_point(size = 3) + 
  geom_xsidedensity(aes(y = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  geom_ysidedensity(aes(x = after_stat(density), fill = C_S, color = C_S, alpha = 0.5)) +
  theme_pubr()

### statistics bray

metadata <- as(sample_data(physeq.Rhizo.norm), "data.frame")
r_OTU_CS <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ C_S,
             data = as(sample_data(physeq.Rhizo.norm), "data.frame"))
r_OTU_Study <- adonis2(phyloseq::distance(physeq.Rhizo.norm, method="bray") ~ Study,
        data = as(sample_data(physeq.Rhizo.norm), "data.frame"))

r_genus_CS <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus, method="bray") ~ C_S,
        data = as(sample_data(physeq.Rhizo.norm_genus), "data.frame"))
r_genus_Study <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus, method="bray") ~ Study,
        data = as(sample_data(physeq.Rhizo.norm_genus), "data.frame"))
r_genus_crop <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus, method="bray") ~ Crop,
                         data = as(sample_data(physeq.Rhizo.norm_genus), "data.frame"))
r_genus_pathogen <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus, method="bray") ~ Crop,
                        data = as(sample_data(physeq.Rhizo.norm_genus), "data.frame"))

r_family_CS <- adonis2(phyloseq::distance(physeq.Rhizo.norm_family, method="bray") ~ C_S,
        data = as(sample_data(physeq.Rhizo.norm_family), "data.frame"))
r_family_Study <- adonis2(phyloseq::distance(physeq.Rhizo.norm_family, method="bray") ~ Study,
                       data = as(sample_data(physeq.Rhizo.norm_family), "data.frame"))

r_order_CS <- adonis2(phyloseq::distance(physeq.Rhizo.norm_order, method="bray") ~ C_S,
        data = as(sample_data(physeq.Rhizo.norm_order), "data.frame"))
r_order_Study <- adonis2(phyloseq::distance(physeq.Rhizo.norm_order, method="bray") ~ Study,
                      data = as(sample_data(physeq.Rhizo.norm_order), "data.frame"))

###########PCOA

registerDoParallel(cores=20)
GP.ord <- ordinate(physeq.Rhizo.norm, "PCoA", "bray")
GP.ord_order <- ordinate(physeq.Rhizo.norm_order, "PCoA", "bray")
GP.ord_family <- ordinate(physeq.Rhizo.norm_family, "PCoA", "bray")
GP.ord_genus <- ordinate(physeq.Rhizo.norm_genus, "PCoA", "bray")

# OTU level
r_OTU_CS
r_OTU_Study
plot_ordination(physeq.Rhizo.norm, GP.ord, type="samples", color= "C_S") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 72462 OTUs")+
  annotate("text", x = -0.2, y = 0.5, hjust = 0, label = "Soil suppressiveness: R2=0.02083, P=0.001") +
  annotate("text", x = -0.2, y = 0.45, hjust = 0, label = "Study: R2=0.54673, P=0.001") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm, GP.ord, type="samples", color="Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 72462 OTUs")+
  annotate("text", x = -0.2, y = 0.5, hjust = 0, label = "Soil suppressiveness: R2=0.02083, P=0.001") +
  annotate("text", x = -0.2, y = 0.45, hjust = 0, label = "Study: R2=0.54673, P=0.001") +
  theme_pubr(border = T) + theme(legend.position = "none")

# Genus level
r_genus_CS
r_genus_Study
plot_ordination(physeq.Rhizo.norm_genus, GP.ord_genus, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 2752 Genera")+
  annotate("text", x = -0.1, y = 0.3, hjust = 0, label = "Soil suppressiveness: R2=0.01587, P=0.001") +
  annotate("text", x = -0.1, y = 0.25, hjust = 0, label = "Study: R2=0.66373, P=0.001") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) +  theme(legend.position = "none")
plot_ordination(physeq.Rhizo.norm_genus, GP.ord_genus, type="samples", color="Study", shape = "C_S") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 2752 Genera")+
  annotate("text", x = -0.1, y = 0.3, hjust = 0, label = "Soil suppressiveness: R2=0.01587, P=0.001") +
  annotate("text", x = -0.1, y = 0.25, hjust = 0, label = "Study: R2=0.66373, P=0.001") +
  theme_pubr(border = T) +  theme(legend.position = "none")

# Family level
r_family_CS
r_family_Study
plot_ordination(physeq.Rhizo.norm_family, GP.ord_family, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 907 Families")+
  annotate("text", x = -0.1, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01348, P=0.001") +
  annotate("text", x = -0.1, y = 0.3, hjust = 0, label = "Study: R2=0.6729, P=0.001") +
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) +  theme(legend.position = "none")
# Order level
r_order_CS
r_order_Study
plot_ordination(physeq.Rhizo.norm_order, GP.ord_order, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = "PCoA, 512 Orders")+
  annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

# Drop 1 taxa >> drop core taxa

# physeq.Rhizo.norm_test <- subset_taxa(physeq.Rhizo.norm, Genus!="Sphingomonas")
# physeq.Rhizo.norm_genus_test <- subset_taxa(physeq.Rhizo.norm_genus, Genus!="Sphingomonas")
# physeq.Rhizo.norm_family_test <- subset_taxa(physeq.Rhizo.norm_family, Family!="Sphingomonadaceae")
# physeq.Rhizo.norm_order_test <- subset_taxa(physeq.Rhizo.norm_order, Genus!="Sphingomonas")

physeq.Rhizo.norm_ex_core <- subset_taxa(physeq.Rhizo.norm, !Genus %in% core_genus )
physeq.Rhizo.norm_genus_ex_core <- subset_taxa(physeq.Rhizo.norm_genus, !Genus %in% core_genus )
physeq.Rhizo.norm_family_ex_core <- subset_taxa(physeq.Rhizo.norm_family, !Family %in% core_family )

# GP.ord_test <- ordinate(physeq.Rhizo.norm_test, "PCoA", "bray")
GP.ord_ex_core <- ordinate(physeq.Rhizo.norm_ex_core, "PCoA", "bray")
GP.ord_genus_ex_core <- ordinate(physeq.Rhizo.norm_genus_ex_core, "PCoA", "bray")
GP.ord_family_ex_core <- ordinate(physeq.Rhizo.norm_family_ex_core, "PCoA", "bray")


# r_CS_test <- adonis2(phyloseq::distance(physeq.Rhizo.norm_test, method="bray") ~ C_S,
#                       data = as(sample_data(physeq.Rhizo.norm_test), "data.frame"))
# r_Study_test <- adonis2(phyloseq::distance(physeq.Rhizo.norm_test, method="bray") ~ Study,
#                          data = as(sample_data(physeq.Rhizo.norm_test), "data.frame"))

r_CS_test_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_ex_core, method="bray") ~ C_S,
                     data = as(sample_data(physeq.Rhizo.norm_ex_core), "data.frame"))
r_Study_test_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_ex_core, method="bray") ~ Study,
                        data = as(sample_data(physeq.Rhizo.norm_ex_core), "data.frame"))

r_CS_test_genus_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus_ex_core, method="bray") ~ C_S,
                             data = as(sample_data(physeq.Rhizo.norm_genus_ex_core), "data.frame"))
r_Study_test_genus_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_genus_ex_core, method="bray") ~ Study,
                                data = as(sample_data(physeq.Rhizo.norm_genus_ex_core), "data.frame"))

r_CS_test_family_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_family_ex_core, method="bray") ~ C_S,
                             data = as(sample_data(physeq.Rhizo.norm_family_ex_core), "data.frame"))
r_Study_test_family_ex_core <- adonis2(phyloseq::distance(physeq.Rhizo.norm_family_ex_core, method="bray") ~ Study,
                                data = as(sample_data(physeq.Rhizo.norm_family_ex_core), "data.frame"))
# plot_ordination(physeq.Rhizo.norm_test, GP.ord_test, type="samples", color="C_S", shape = "Study") +
#   scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
#   # labs(title = "PCoA, 512 Orders")+
#   # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
#   # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
#   scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm_ex_core, GP.ord_ex_core, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm, GP.ord, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

# CAP drop 1
# GP.cap_test <- ordinate(physeq.Rhizo.norm_test, "CAP", "bray", ~C_S)
# GP.cap_test <- ordinate(physeq.Rhizo.norm_genus_test, "CAP", "bray", ~C_S)
# 
# GP.cap_family_test <- ordinate(physeq.Rhizo.norm_family_test, "CAP", "bray", ~C_S)
GP.cap
GP.cap_test_ex_core <- ordinate(physeq.Rhizo.norm_ex_core, "CAP", "bray", ~C_S)
GP.cap_genus_test_ex_core <- ordinate(physeq.Rhizo.norm_genus_ex_core, "CAP", "bray", ~C_S)
GP.cap_family_test_ex_core <- ordinate(physeq.Rhizo.norm_family_ex_core, "CAP", "bray", ~C_S)


plot_ordination(physeq.Rhizo.norm, GP.cap, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm_ex_core, GP.cap_test_ex_core, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")


plot_ordination(physeq.Rhizo.norm_genus, GP.cap_genus, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm_genus_ex_core, GP.cap_genus_test_ex_core, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")


plot_ordination(physeq.Rhizo.norm_family, GP.cap_family, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

plot_ordination(physeq.Rhizo.norm_family_ex_core, GP.cap_family_test_ex_core, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  # labs(title = "PCoA, 512 Orders")+
  # annotate("text", x = -0.2, y = 0.35, hjust = 0, label = "Soil suppressiveness: R2=0.01309, P=0.001") +
  # annotate("text", x = -0.2, y = 0.3, hjust = 0, label = "Study: R2=0.65893, P=0.001") + 
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")

# loop
physeq_list <- list()
physeq_list$all <- list()
physeq_list$all$OTU <- physeq.Rhizo
physeq_list$all$genus <- physeq.Rhizo_genus
physeq_list$all$family <- physeq.Rhizo_family
physeq_list$all$order <- physeq.Rhizo_order

physeq_list$ex <- list()
physeq_list$ex$OTU <- subset_samples(physeq.Rhizo, !Study %in% c("Ruth-zwaagdijk", "PRJNA578725"))
physeq_list$ex$genus <- subset_samples(physeq.Rhizo_genus, !Study %in% c("Ruth-zwaagdijk", "PRJNA578725"))
physeq_list$ex$family <- subset_samples(physeq.Rhizo_family, !Study %in% c("Ruth-zwaagdijk", "PRJNA578725"))
physeq_list$ex$order <- subset_samples(physeq.Rhizo_order, !Study %in% c("Ruth-zwaagdijk", "PRJNA578725"))

pcoa_list <- list()
pcoa_list$all <- list()
pcoa_list$ex <- list()

cap_list <- list()
cap_list$all <- list()
cap_list$ex <- list()

pcoa_sub_list <- list()
pcoa_sub_list$all <- list()
pcoa_sub_list$ex <- list()

cap_sub_list <- list()
cap_sub_list$all <- list()
cap_sub_list$ex <- list()

adonis_list <- list()
adonis_list$all <- list()
adonis_list$ex <- list()

adonis_list_sub <- list()
adonis_list_sub$all <- list()
adonis_list_sub$ex <- list()

for (s in c("all", "ex")){
for (t in c("OTU", "genus", "family", "order")){
p <- physeq_list[[s]][[t]]
physeq_norm <- CSS_norm(p)
adonis_CS <- adonis2(phyloseq::distance(physeq_norm, method="bray") ~ C_S,
                      data = as(sample_data(physeq_norm), "data.frame"))
adonis_Study <- adonis2(phyloseq::distance(physeq_norm, method="bray") ~ Study,
                         data = as(sample_data(physeq_norm), "data.frame"))
adonis_list[[s]][[t]][["CS"]] <- adonis_CS
adonis_list[[s]][[t]][["Study"]] <- adonis_Study

GP.pcoa <- ordinate(physeq_norm, "PCoA", "bray")
GP.cap <- ordinate(physeq_norm, "CAP", "bray", ~ C_S)

pcoa <- plot_ordination(physeq_norm, GP.pcoa, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = paste("PCoA", ntaxa(physeq_norm), t, s))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")
pcoa_list[[s]][[t]] <- pcoa

cap <- plot_ordination(physeq_norm, GP.cap, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = paste("CAP", ntaxa(physeq_norm), t, s))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")
cap_list[[s]][[t]] <- cap

# Remove Sphingomonadales
p_sub <- subset_taxa(p, Order!="Sphingomonadales")
physeq_norm_sub <- CSS_norm(p_sub)
adonis_CS_sub <- adonis2(phyloseq::distance(physeq_norm_sub, method="bray") ~ C_S,
                     data = as(sample_data(physeq_norm_sub), "data.frame"))
adonis_Study_sub <- adonis2(phyloseq::distance(physeq_norm_sub, method="bray") ~ Study,
                        data = as(sample_data(physeq_norm_sub), "data.frame"))
adonis_list_sub[[s]][[t]][["CS"]] <- adonis_CS_sub
adonis_list_sub[[s]][[t]][["Study"]] <- adonis_Study_sub

GP.pcoa_sub <- ordinate(physeq_norm_sub, "PCoA", "bray")
GP.cap_sub <- ordinate(physeq_norm_sub, "CAP", "bray", ~ C_S)

pcoa_sub <- plot_ordination(physeq_norm_sub, GP.pcoa_sub, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = paste("PCoA", ntaxa(physeq_norm_sub), t, s, "Sphingomonadales excluded"))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")
pcoa_sub_list[[s]][[t]] <- pcoa_sub

cap_sub <- plot_ordination(physeq_norm_sub, GP.cap_sub, type="samples", color="C_S", shape = "Study") +
  scale_shape_manual(values = 0:20) + geom_point(alpha=0.65) +
  labs(title = paste("CAP", ntaxa(physeq_norm_sub), t, s, "Sphingomonadales excluded"))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B")) + theme_pubr(border = T) + theme(legend.position = "none")
cap_sub_list[[s]][[t]] <- cap_sub
}}


  for (t in c("OTU", "genus", "family", "order")){
    for (s in c("all", "ex")){
    print(ggarrange(pcoa_list[[s]][[t]], cap_list[[s]][[t]],
              pcoa_sub_list[[s]][[t]], cap_sub_list[[s]][[t]]))
    }}


ggarrange(pcoa_list[["all"]][["OTU"]], cap_list[["all"]][["OTU"]],
          pcoa_sub_list[["all"]][["OTU"]], cap_sub_list[["all"]][["OTU"]])

ggarrange(pcoa_list[["all"]][["genus"]], cap_list[["all"]][["genus"]],
          pcoa_sub_list[["all"]][["genus"]], cap_sub_list[["all"]][["genus"]])

ggarrange(pcoa_list[["all"]][["family"]], cap_list[["all"]][["family"]],
          pcoa_sub_list[["all"]][["family"]], cap_sub_list[["all"]][["family"]])

ggarrange(pcoa_list[["all"]][["order"]], cap_list[["all"]][["order"]],
          pcoa_sub_list[["all"]][["order"]], cap_sub_list[["all"]][["order"]])


ggarrange(pcoa_list[["ex"]][["OTU"]], cap_list[["ex"]][["OTU"]],
          pcoa_sub_list[["ex"]][["OTU"]], cap_sub_list[["ex"]][["OTU"]])

ggarrange(pcoa_list[["ex"]][["genus"]], cap_list[["ex"]][["genus"]],
          pcoa_sub_list[["ex"]][["genus"]], cap_sub_list[["ex"]][["genus"]])

ggarrange(pcoa_list[["ex"]][["family"]], cap_list[["ex"]][["family"]],
          pcoa_sub_list[["ex"]][["family"]], cap_sub_list[["ex"]][["family"]])

ggarrange(pcoa_list[["ex"]][["order"]], cap_list[["ex"]][["order"]],
          pcoa_sub_list[["ex"]][["order"]], cap_sub_list[["ex"]][["order"]])


adonis_list$all$OTU
adonis_list_sub$all$OTU




