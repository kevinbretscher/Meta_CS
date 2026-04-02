library(phytools)
library(ape)
ID_comb_all$GenBank_ID
trait_enriched_OTU_ID_all <- enriched_OTU_ID_all$GenBank_ID

# Create a named vector with default value 0 (not enriched/depleted)
trait_vector <- setNames(rep(0, length(ID_comb_all$GenBank_ID)), ID_comb_all$GenBank_ID)

# Assign '1' to enriched OTUs
trait_vector[enriched_OTU_ID_all$GenBank_ID] <- 1

# Assign '0' to depleted OTUs (already default but ensures explicit labeling)
trait_vector[depleted_OTU_ID_all$GenBank_ID] <- 0
print(trait_vector)

enriched_OTU_ID_all$GenBank_ID[enriched_OTU_ID_all$GenBank_ID %in% depleted_OTU_ID_all$GenBank_ID]
## Some OTUs are both enriched and depleted in different studies
sum(depleted_OTU_ID_all$GenBank_ID %in% enriched_OTU_ID_all$GenBank_ID)
duplicated_IDs_all <- depleted_OTU_ID_all$GenBank_ID[depleted_OTU_ID_all$GenBank_ID %in% enriched_OTU_ID_all$GenBank_ID]

trait_vector[duplicated_IDs_all]


## filter also based on foldchange
sum(unique(c(enriched_OTU$Taxa[enriched_OTU$logFC>1],enriched_OTU_ANCOM$taxon[enriched_OTU_ANCOM$lfc_C_SS>1])) %in%
unique(c(depleted_OTU$Taxa[depleted_OTU$logFC < (-1)],depleted_OTU_ANCOM$taxon[depleted_OTU_ANCOM$lfc_C_SS < (-1)])))

filter_lfc <- 1
combined_enriched_IDs_filtered <- unique(c(enriched_OTU$Taxa[enriched_OTU$logFC > filter_lfc],
                                           enriched_OTU_ANCOM$taxon[enriched_OTU_ANCOM$lfc_C_SS > filter_lfc]))
combined_depleted_IDs_filtered <- unique(c(depleted_OTU$Taxa[depleted_OTU$logFC < (-filter_lfc)],
                                           depleted_OTU_ANCOM$taxon[depleted_OTU_ANCOM$lfc_C_SS < (-filter_lfc)]))

duplicated_IDs <- combined_depleted_IDs_filtered[combined_depleted_IDs_filtered%in%combined_enriched_IDs_filtered]
test <- enriched_OTU_ANCOM[enriched_OTU_ANCOM$taxon %in% duplicated_IDs,]
test_1 <- depleted_OTU_ANCOM[depleted_OTU_ANCOM$taxon %in% duplicated_IDs,]

### mixing ancom and metagenomeseq is not a good idea, can check it separately for the assessment of phylogenetic signals
# ANCOM
enriched_OTU_ANCOM_IDs <- enriched_OTU_ANCOM$taxon
depleted_OTU_ANCOM_IDs <- depleted_OTU_ANCOM$taxon
enriched_ANCOM_trait_vector <- setNames(rep(1, length(enriched_OTU_ANCOM_IDs)), enriched_OTU_ANCOM_IDs)
depleted_ANCOM_trait_vector <- setNames(rep(0, length(depleted_OTU_ANCOM_IDs)), depleted_OTU_ANCOM_IDs)
ANCOM_trait_vector <- c(enriched_ANCOM_trait_vector,depleted_ANCOM_trait_vector)
sum(enriched_OTU_ANCOM_IDs%in%depleted_OTU_ANCOM_IDs)
dup_ANCOM <- enriched_OTU_ANCOM_IDs[enriched_OTU_ANCOM_IDs%in%depleted_OTU_ANCOM_IDs]
names(ANCOM_trait_vector)[names(ANCOM_trait_vector) %in% dup_ANCOM]

# remove IDs duplicated
ANCOM_trait_vector <- ANCOM_trait_vector[!names(ANCOM_trait_vector) %in% dup_ANCOM]

# trait_df <- data.frame(OTU = names(ANCOM_trait_vector), Value = ANCOM_trait_vector)
# 
# trait_df <- trait_df %>%
#   group_by(OTU) %>%
#   summarise(Final_Value = as.integer(mean(Value) > 0.5), .groups = "drop")

#metagenomeseq

enriched_OTU_metagenomeseq_IDs <- enriched_OTU$Taxa
depleted_OTU_metagenomeseq_IDs <- depleted_OTU$Taxa
enriched_metagenomeseq_trait_vector <- setNames(rep(1, length(enriched_OTU_metagenomeseq_IDs)), enriched_OTU_metagenomeseq_IDs)
depleted_metagenomeseq_trait_vector <- setNames(rep(0, length(depleted_OTU_metagenomeseq_IDs)), depleted_OTU_metagenomeseq_IDs)
metagenomeseq_trait_vector <- c(enriched_metagenomeseq_trait_vector, depleted_metagenomeseq_trait_vector)
sum(enriched_OTU_metagenomeseq_IDs%in%depleted_OTU_metagenomeseq_IDs)
dup_metagenomeseq <- enriched_OTU_metagenomeseq_IDs[enriched_OTU_metagenomeseq_IDs%in%depleted_OTU_metagenomeseq_IDs]
metagenomeseq_trait_vector <- metagenomeseq_trait_vector[!names(metagenomeseq_trait_vector) %in% dup_metagenomeseq]

### only take the sig OTUs belonging to the core genera, and filter based on foldchange
core_genus_C <- read.csv("../diversity/diversity/C_core_genus.csv")
core_genus_S <- read.csv("../diversity/diversity/S_core_genus.csv")

# handle uncultured
core_genus_C$Genus[core_genus_C$Genus == "uncultured"] <- paste("uncultured", 
                                                                core_genus_C$Family[core_genus_C$Genus == "uncultured"], sep = "_")
core_genus_C$Genus[core_genus_C$Genus == "uncultured_uncultured"] <- paste("uncultured", 
                                                                core_genus_C$Order[core_genus_C$Genus == "uncultured_uncultured"], sep = "_")
core_genus_C$Genus[core_genus_C$Genus == "uncultured_uncultured"] <- paste("uncultured", 
                                                                           core_genus_C$Class[core_genus_C$Genus == "uncultured_uncultured"], sep = "_")

core_genus_S$Genus[core_genus_S$Genus == "uncultured"] <- paste("uncultured", 
                                                                core_genus_S$Family[core_genus_S$Genus == "uncultured"], sep = "_")
core_genus_S$Genus[core_genus_S$Genus == "uncultured_uncultured"] <- paste("uncultured", 
                                                                           core_genus_S$Order[core_genus_S$Genus == "uncultured_uncultured"], sep = "_")
core_genus_S$Genus[core_genus_S$Genus == "uncultured_uncultured"] <- paste("uncultured", 
                                                                           core_genus_S$Class[core_genus_S$Genus == "uncultured_uncultured"], sep = "_")

core_genus <- unique(c(core_genus_C$Genus, core_genus_S$Genus))

find_uncultured <- enriched_OTU$Genus == "uncultured" & !is.na(enriched_OTU$Genus)
enriched_OTU$Genus[find_uncultured] <- paste("uncultured", 
                                             enriched_OTU$Family[find_uncultured], sep = "_")

# write it into a function
update_uncultured_genus <- function(df) {
  find_uncultured <- df$Genus == "uncultured" & !is.na(df$Genus)
  df$Genus[find_uncultured] <- paste("uncultured", df$Family[find_uncultured], sep = "_")
  return(df)
}

depleted_OTU <- update_uncultured_genus(depleted_OTU)
enriched_OTU_ANCOM <- update_uncultured_genus(enriched_OTU_ANCOM)
depleted_OTU_ANCOM <- update_uncultured_genus(depleted_OTU_ANCOM)

sum(enriched_OTU$Genus %in% core_genus)
sum(depleted_OTU$Genus %in% core_genus)

sum(enriched_OTU_ANCOM$Genus %in% core_genus)
sum(depleted_OTU_ANCOM$Genus %in% core_genus)


core_enriched_metagenomeseq <- enriched_OTU[enriched_OTU$Genus %in% core_genus,]
core_enriched_metagenomeseq <- core_enriched_metagenomeseq[core_enriched_metagenomeseq$logFC > 1,]
core_depleted_metagenomeseq <- depleted_OTU[depleted_OTU$Genus %in% core_genus,]
core_depleted_metagenomeseq <- core_depleted_metagenomeseq[core_depleted_metagenomeseq$logFC < (-1),]
# check duplicated
sum(core_enriched_metagenomeseq$Taxa %in% core_depleted_metagenomeseq$Taxa)
sum(core_depleted_metagenomeseq$Taxa %in% core_enriched_metagenomeseq$Taxa)
dup_metagenomeseq_core <- unique(core_enriched_metagenomeseq$Taxa
                                 [core_enriched_metagenomeseq$Taxa %in% core_depleted_metagenomeseq$Taxa])

core_enriched_ANCOM <- enriched_OTU_ANCOM[enriched_OTU_ANCOM$Genus %in% core_genus,]
core_enriched_ANCOM <- core_enriched_ANCOM[core_enriched_ANCOM$lfc_C_SS > 1,]
core_depleted_ANCOM <- depleted_OTU_ANCOM[depleted_OTU_ANCOM$Genus %in% core_genus,]
core_depleted_ANCOM <- core_depleted_ANCOM[core_depleted_ANCOM$lfc_C_SS < (-1),]
# check duplicated
sum(core_enriched_ANCOM$taxon %in% core_depleted_ANCOM$taxon)
sum(core_depleted_ANCOM$taxon %in% core_enriched_ANCOM$taxon)
dup_ANCOM_core <- unique(core_enriched_ANCOM$taxon[core_enriched_ANCOM$taxon %in% core_depleted_ANCOM$taxon])

ID_core_metagenomeseq <- Silva_to_Genebank(unique(c(core_enriched_metagenomeseq$Taxa, core_depleted_metagenomeseq$Taxa)))
ID_core_ANCOM <- Silva_to_Genebank(unique(c(core_enriched_ANCOM$taxon, core_depleted_ANCOM$taxon)))

#vec <- c(core_enriched_metagenomeseq$Taxa, core_depleted_metagenomeseq$Taxa)
#vec[!duplicated(vec) & !duplicated(vec, fromLast = TRUE)]

fetch_sequence_regions(ID_core_metagenomeseq, file_output = "all_sig_q0.1_fc2/metagenomeseq/metagenomeseq_core_sig_seqs.fasta")
fetch_sequence_regions(ID_core_ANCOM, file_output = "all_sig_q0.1_fc2/ancom/ancom_core_sig.fasta")

# make alighment with muscle and tree with iqtree2

# define trait vector for the core sig OTUs
# ANCOM
trait_vector_core_ANCOM_enriched <- setNames(rep(1, length(core_enriched_ANCOM$taxon)), core_enriched_ANCOM$taxon)
trait_vector_core_ANCOM_depleted <- setNames(rep(0, length(core_depleted_ANCOM$taxon)), core_depleted_ANCOM$taxon)
trait_vector_core_ANCOM <- c(trait_vector_core_ANCOM_enriched, trait_vector_core_ANCOM_depleted)
length(trait_vector_core_ANCOM[!names(trait_vector_core_ANCOM) %in% 
                          names(trait_vector_core_ANCOM)[duplicated(names(trait_vector_core_ANCOM)) | duplicated(names(trait_vector_core_ANCOM), fromLast = TRUE)]])
# metagenomeseq
trait_vector_core_metagenomeseq_enriched <- setNames(rep(1, length(core_enriched_metagenomeseq$Taxa)), core_enriched_metagenomeseq$Taxa)
trait_vector_core_metagenomeseq_depleted <- setNames(rep(0, length(core_depleted_metagenomeseq$Taxa)), core_depleted_metagenomeseq$Taxa)
trait_vector_core_metagenomeseq<- c(trait_vector_core_metagenomeseq_enriched, trait_vector_core_metagenomeseq_depleted)

library(phytools)
library(ape)

# Read tree
tree_ancom <- read.tree("all_sig_q0.1_fc2/ancom/ancom_core_sig.afa.treefile")
tip.label(tree_ancom) <- sub("\\..*", "", tip.label(tree_ancom))

tree_metagenomeseq <- read.tree("all_sig_q0.1_fc2/metagenomeseq/metagenomeseq_core_sig.afa.treefile")
tip.label(tree_metagenomeseq) <- sub("\\..*", "", tip.label(tree_metagenomeseq))

# Extract the trait (e.g., enrichment as continuous or binary)
names(trait_vector_core_ANCOM) <- sub("\\..*", "", names(trait_vector_core_ANCOM))
names(trait_vector_core_metagenomeseq) <- sub("\\..*", "", names(trait_vector_core_metagenomeseq))

# Compute Pagel’s λ
lambda_result <- phylosig(tree_ancom, trait_vector_core_ANCOM, method="lambda", test=TRUE)
lambda_result_ancom <- lambda_result
print(lambda_result_ancom)

lambda_result_metagenomeseq <- phylosig(tree_metagenomeseq, trait_vector_core_metagenomeseq, method="lambda", test=TRUE)
print(lambda_result_metagenomeseq)

# fun part: let's visualize the tree of sig OTUs
iqtree_ancom <- read.iqtree("all_sig_q0.1_fc2/ancom/ancom_core_sig.afa.treefile")

tree_data_ancom$label_renamed <- sub("\\..*", "", tree_data_ancom$label)
names(trait_vector_core_ANCOM_enriched) <- sub("\\..*", "", names(trait_vector_core_ANCOM_enriched))
names(trait_vector_core_ANCOM_depleted) <- sub("\\..*", "", names(trait_vector_core_ANCOM_depleted))
tree_data_ancom$group <- "Other"
tree_data_ancom$group[tree_data_ancom$label_renamed %in% names(trait_vector_core_ANCOM_enriched)] <- "Enriched"
tree_data_ancom$group[tree_data_ancom$label_renamed %in% names(trait_vector_core_ANCOM_depleted)] <- "Depleted"
tree_data_ancom$group[tree_data_ancom$label_renamed =="NR_102775"] <- "Outgroup Thermotoga"

library(ape)
library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)
library(treedataverse)
library(ggnewscale)
library(tidytree)
library(ggpubr)
library(treeplyr)
library(dplyr)

p <- ggtree(iqtree_ancom,layout= "fan") %<+% tree_data_ancom +  # Attach metadata to the tree
  geom_tiplab(aes(label = label, color = group), show.legend = FALSE, size= 2, align=TRUE, linesize=.5) +  # Color labels
  geom_tippoint(aes(color = group), size = 1) +  # Color tips
  scale_color_manual(values = c("Enriched" = "red","Depleted" = "blue", "Other" = "black"))

core_enriched_ANCOM$label_renamed <- sub("\\..*", "", core_enriched_ANCOM$taxon)
tree_data_ancom_enriched <- data.frame(label_renamed = core_enriched_ANCOM$label_renamed,
                                       Phylum = core_enriched_ANCOM$Phylum,
                                       Order = core_enriched_ANCOM$Order,
                                       Family = core_enriched_ANCOM$Family,
                                       Genus = core_enriched_ANCOM$Genus)
core_depleted_ANCOM$label_renamed <- sub("\\..*", "", core_depleted_ANCOM$taxon)
tree_data_ancom_depleted <- data.frame(label_renamed = core_depleted_ANCOM$label_renamed,
                                       Phylum = core_depleted_ANCOM$Phylum,
                                       Order = core_depleted_ANCOM$Order,
                                       Family = core_depleted_ANCOM$Family,
                                       Genus = core_depleted_ANCOM$Genus)
tree_data_ancom_phylum <- unique(rbind(tree_data_ancom_enriched, tree_data_ancom_depleted))
rownames(tree_data_ancom_phylum) <- tree_data_ancom_phylum$label_renamed
tree_data_ancom_phylum[2]

tip.label(iqtree_ancom) <- sub("\\..*", "", tip.label(iqtree_ancom))
tree_data_ancom <- iqtree_ancom %>% fortify()

gheatmap(p, tree_data_ancom_phylum[2], offset=.8, width=.2,
         colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")


to_drop_ancom <- c("HE583116",
"FR872466",
"FPLS01007673",
"FPLS01015191",
"FR687516",
"CP000927",
"FPLS01038872",
"FPLS01007672",
"FPLS01005249",
"FPLS01005733",
"FPLL01001173",
"FPLS01021518",
"FR872477",
"NEJH01000025",
"FPLS01036208",
"HE583140",
"CP033019",
"AJLS01000166",
"AY281355",
"FPLS01051226",
"FPLS01029794",
"FPLS01057109",
"JN566182",
"LMEZ01000006",
"FPLS01009538",
"FPLS01029226",
"FPLS01027161",
"FPLS01062878",
"FPLS01000104")

iqtree.reduced_ancom <- drop.tip(iqtree_ancom, to_drop_ancom)
iqtree.reduced_ancom <- root(iqtree.reduced_ancom, outgroup = "NR_102775", edgelabel = TRUE)
#tree_data_ancom_phylum <- left_join(tree_data_ancom_phylum, tree_data_ancom, by = "label_renamed")
cb_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#999999", "#882255", "#44AA99", "#AA4499", "#DDCC77", "#117733", "#332288",
  "#88CCEE", "#661100", "#6699CC", "#AA4466"
)
tree_data_ancom_phylum$Family_highlight <- "Other"
to_highlight <- c("Chitinophagaceae",
                  "Sphingomonadaceae",
                  "Bacillaceae",
                  "Comamonadaceae",
                  "Rhizobiaceae",
                  "Sphingobacteriaceae",
                  "Oxalobacteraceae",
                  "Pseudomonadaceae",
                  "Nitrospiraceae")
tree_data_ancom_phylum$Family_highlight[tree_data_ancom_phylum$Family%in%to_highlight] <- tree_data_ancom_phylum$Family[tree_data_ancom_phylum$Family%in%to_highlight]

selected_genera <- c("Bacillus",
                     "Flavobacterium",
                     "Streptomyces",
                     "Chitinophaga",
                     "Pseudomonas",
                     "Novosphingobium",
                     "Sphingomonas")
selected_labels <- tree_data_ancom_phylum$label_renamed[tree_data_ancom_phylum$Phylum %in%
                                                          c("Proteobacteria", "Firmicutes", "Baceroidota")]
selected_genus_info <- all_genus_info[all_genus_info$Genus %in% selected_genera,]
rownames(selected_genus_info) <- selected_genus_info$label_renamed

tree_data_ancom_phylum$Genus_highlight <- "Other"
tree_data_ancom_phylum$Genus_highlight[tree_data_ancom_phylum$Genus%in%selected_genera] <- tree_data_ancom_phylum$Genus[tree_data_ancom_phylum$Genus%in%selected_genera]
# Define the desired order
desired_order <- c("Bacillus", "Streptomyces", "Flavobacterium", "Chitinophaga",
                   "Pseudomonas", "Novosphingobium", "Sphingomonas", "Other")

# Convert Genus_highlight to a factor with the desired order
tree_data_ancom_phylum <- tree_data_ancom_phylum %>%
  mutate(Genus_highlight = factor(Genus_highlight, levels = desired_order))

ggtree(iqtree.reduced_ancom,layout= "circular") %<+% tree_data_ancom_phylum +  # Attach metadata to the tree
  geom_tippoint(aes(color = Family, na.translate = FALSE), size = 0.75)

p <- ggtree(iqtree.reduced_ancom,layout= "circular") %<+% tree_data_ancom_phylum +  # Attach metadata to the tree
  geom_tippoint(aes(color = Genus_highlight, na.translate = FALSE), size = 0.75) +  # Color tips
  scale_color_manual(name = "Genus",
                    values = c("Bacillus" = "#E69F00",
                               "Flavobacterium" = "#56B4E9",
                               "Streptomyces" = "#009E73",
                               "Chitinophaga"= "#F0E442",
                               "Pseudomonas" = "#0072B2",
                               "Novosphingobium" = "#D55E00",
                               "Sphingomonas" = "#CC79A7",
                               "Other" = "#999999"
                               ),
                    na.translate = FALSE)

p1 <- gheatmap(p, tree_data_ancom_phylum[2], offset=0.1, width=.1, color = NA, colnames = F) +
  scale_fill_manual(name = "Phylum",
                    values = c("Proteobacteria" = "#E64B35FF",
                               "Actinobacteriota" = "#4DBBD5FF",
                               "Acidobacteriota" = "#00A087FF",
                               "Bacteroidota" = "#3C5488FF",
                               "Firmicutes" = "#F39B7FFF",
                               "Verrucomicrobiota"= "#8491B4FF",
                               "Chloroflexi" = "#91D1C2FF",
                               "Planctomycetota" = "#DC0000FF",
                               "Myxococcota"="#8c564b",
                               "Nitrospirota"="#bcbd22",
                               "Gemmatimonadota"= "#1f77b4",
                               "Bdellovibrionota"="#ff7f0e"),
                    na.translate = FALSE)
# pdf("all_sig_q0.1_fc2/ancom/tree.pdf", paper="a4")
# print(p1)
# dev.off()
# pal_npg("nrc")(10)
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF"
# [9] "#7E6148FF" "#B09C85FF"
# 
# "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#fdae61", "#5e3c99"

# Count unique occurrences of "study" for each "taxon"
count_core_enriched_ANCOM <- core_enriched_ANCOM %>%
  group_by(taxon) %>%
  summarise(n_studies = n_distinct(study))
count_core_depleted_ANCOM <- core_depleted_ANCOM %>%
  group_by(taxon) %>%
  summarise(n_studies = n_distinct(study))
count_core_enriched_ANCOM$taxon <- sub("\\..*", "", count_core_enriched_ANCOM$taxon)
count_core_depleted_ANCOM$taxon <- sub("\\..*", "", count_core_depleted_ANCOM$taxon)
tree_data_ancom_enriched_studies <- data.frame(label_renamed = count_core_enriched_ANCOM$taxon,
                                               n_studies = count_core_enriched_ANCOM$n_studies,
                                               suppressive_soil = "Enriched")
tree_data_ancom_depleted_studies <- data.frame(label_renamed = count_core_depleted_ANCOM$taxon,
                                               n_studies = count_core_depleted_ANCOM$n_studies,
                                               suppressive_soil = "Depleted")
tree_data_ancom_studies <- unique(rbind(tree_data_ancom_enriched_studies, tree_data_ancom_depleted_studies))
rownames(tree_data_ancom_phylum) <- tree_data_ancom_phylum$label_renamed

ancom_lfc_enriched <- data.frame(label_renamed = core_enriched_ANCOM$label_renamed,
                                 lfc = core_enriched_ANCOM$lfc_C_SS,
                                 condition = "Enriched")
ancom_lfc_depleted <- data.frame(label_renamed = core_depleted_ANCOM$label_renamed,
                                 lfc = core_depleted_ANCOM$lfc_C_SS,
                                 condition = "Depleted")
tree_data_ancom_lfc <- unique(rbind(ancom_lfc_enriched, ancom_lfc_depleted))
#rownames(tree_data_ancom_lfc) <- tree_data_ancom_lfc$label_renamed


library(ggtreeExtra)
p2 <- p1 + new_scale_fill() +  geom_fruit(
  data = tree_data_ancom_studies, 
  geom = geom_bar, 
  mapping = aes(y = label_renamed, x = n_studies, fill = suppressive_soil), 
  stat = "identity",
  orientation = "y", 
  offset = 0.4,
  axis.params=list(
    axis       = "x",
    text.size  = 2,
    hjust      = 1,
    vjust      = 0.5,
    nbreak     = 2,
  ),
  grid.params=list()
) +
  scale_fill_manual(name = "Suppressive soils",
                    values =  c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
  


p3 <- p2 +  new_scale_fill() +
  geom_fruit(
    data = tree_data_ancom_lfc, 
    geom = geom_bar, 
    mapping = aes(y = label_renamed, x = lfc, fill = condition),
    stat = "identity",
    offset = 0.3,
    pwidth = 0.2,
    axis.params=list(
      axis       = "x",
      text.size  = 2,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 2,
    ),
    grid.params=list()) +
  scale_fill_manual(name="Log2FC",values=c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
pdf("all_sig_q0.1_fc2/ancom/tree_annot.pdf", paper="a4")
print(p3)
dev.off()

library(RColorBrewer)
palette_21 <- colorRampPalette(brewer.pal(12, "Set3"))(21) # Expanding Set1 to 10 colors

samples_to_drop <- tree_data_ancom_phylum$label_renamed[tree_data_ancom_phylum$Phylum!="Proteobacteria"]
sub_tree <- drop.tip(iqtree.reduced_ancom, samples_to_drop)
#sub_tree <- root(sub_tree, outgroup = "NR_102775", edgelabel = TRUE)
p_proteobacteria <- ggtree(sub_tree,layout= "circular")
tree_data_ancom_proteobacteria <- tree_data_ancom_phylum[tree_data_ancom_phylum$Phylum=="Proteobacteria",]
tree_data_ancom_proteobacteria$Family_highlight <- "Other"
is_top10_proteobacteria <- tree_data_ancom_proteobacteria$Family %in% names(sort(table(tree_data_ancom_proteobacteria$Family), decreasing = T)[1:10])
tree_data_ancom_proteobacteria$Family_highlight[is_top10_proteobacteria] <- tree_data_ancom_proteobacteria$Family[is_top10_proteobacteria]
tree_data_ancom_proteobacteria$Family_highlight <- factor(tree_data_ancom_proteobacteria$Family_highlight, levels = 
                                                            c(names(sort(table(tree_data_ancom_proteobacteria$Family), decreasing = T))[1:10], "Other"))

p_proteobacteria1 <- gheatmap(p_proteobacteria, tree_data_ancom_proteobacteria[6], offset=0.2, width=.1, color = NA, colnames = F) +
 scale_fill_manual(name = "Family", 
                   values = c(palette_20[c(1:3,5:11)], "#D4D2CC"), 
                   breaks =c(names(sort(table(tree_data_ancom_proteobacteria$Family), decreasing = T))[1:10], "Other"),
                   na.translate = F)

p_proteobacteria2 <- p_proteobacteria1 + new_scale_fill() +
  geom_fruit(
    data = tree_data_ancom_lfc, 
    geom = geom_bar, 
    mapping = aes(y = label_renamed, x = lfc, fill = condition),
    stat = "identity",
    offset = 0.1,
    pwidth = 0.2) +
  scale_fill_manual(name="Suppressive soils",values=c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
p_proteobacteria2
pdf("all_sig_q0.1_fc2/ancom/tree_annot_proteobacteria.pdf", paper="a4")
print(p_proteobacteria2)
dev.off()


samples_to_drop <- tree_data_ancom_phylum$label_renamed[tree_data_ancom_phylum$Phylum!="Bacteroidota"]
sub_tree <- drop.tip(iqtree.reduced_ancom, samples_to_drop)
#sub_tree <- root(sub_tree, outgroup = "NR_102775", edgelabel = TRUE)
p_bacteroidetes <- ggtree(sub_tree,layout= "circular") %<+% tree_data_ancom_bacteroidetes +
geom_tippoint(aes(color = Genus, na.translate = FALSE), size = 1) +  # Color tips
  scale_color_manual(name = "Genus",
                     values = c("Chitinophaga" = cb_palette[1],
                                "Flavobacterium" = cb_palette[2],
                                "Flavisolibacter" = cb_palette[3],
                                "Pedobacter" = cb_palette[4] ,
                                "Ferruginibacter"= cb_palette[5],
                                "Edaphobaculum"= cb_palette[6],
                                "Terrimonas"= cb_palette[7],
                                "Niastella"= cb_palette[9],
                                "uncultured_Chitinophagaceae"= cb_palette[8]),
                     na.translate = F)
tree_data_ancom_bacteroidetes <- tree_data_ancom_phylum[tree_data_ancom_phylum$Phylum=="Bacteroidota",]
tree_data_ancom_bacteroidetes$Family_highlight <- "Other"
is_top10_bacteroidetes <- tree_data_ancom_bacteroidetes$Family %in% names(sort(table(tree_data_ancom_bacteroidetes$Family), decreasing = T)[1:10])
tree_data_ancom_bacteroidetes$Family_highlight[is_top10_bacteroidetes] <- tree_data_ancom_bacteroidetes$Family[is_top10_bacteroidetes]
tree_data_ancom_bacteroidetes$Family_highlight <- factor(tree_data_ancom_bacteroidetes$Family_highlight, levels = 
                                                            c(names(sort(table(tree_data_ancom_bacteroidetes$Family), decreasing = T))[1:10], "Other"))

p_bacteroidetes1 <- gheatmap(p_bacteroidetes, tree_data_ancom_bacteroidetes[6], offset=0.2, width=.1, color = NA, colnames = F) +
  scale_fill_brewer(name = "Family", 
                    palette = "Set3", 
                    na.translate = F)
p_bacteroidetes2 <- p_bacteroidetes1 + new_scale_fill() +
  geom_fruit(
    data = tree_data_ancom_lfc, 
    geom = geom_bar, 
    mapping = aes(y = label_renamed, x = lfc, fill = condition),
    stat = "identity",
    offset = 0.1,
    pwidth = 0.2) +
  scale_fill_manual(name="Suppressive soils",values=c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
p_bacteroidetes2
pdf("all_sig_q0.1_fc2/ancom/tree_annot_bacteroidetes.pdf", paper="a4")
print(p_bacteroidetes2)
dev.off()

samples_to_drop <- tree_data_ancom_phylum$label_renamed[tree_data_ancom_phylum$Phylum!="Firmicutes"]
sub_tree <- drop.tip(iqtree.reduced_ancom, samples_to_drop)
tree_data_ancom_Firmicutes <- tree_data_ancom_phylum[tree_data_ancom_phylum$Phylum=="Firmicutes",]

#sub_tree <- root(sub_tree, outgroup = "NR_102775", edgelabel = TRUE)
p_Firmicutes <- ggtree(sub_tree,layout= "circular") %<+% tree_data_ancom_Firmicutes +
  geom_tippoint(aes(color = Genus, na.translate = FALSE), size = 1) +  # Color tips
  scale_color_manual(name = "Genus",
                     values = c("Bacillus" = cb_palette[1],
                                "Paenibacillus" = cb_palette[2],
                                "Tumebacillus" = cb_palette[3]),
                     na.translate = F)

p_Firmicutes1 <- gheatmap(p_Firmicutes, tree_data_ancom_Firmicutes[4], offset=0.2, width=.1, color = NA, colnames = F) +
  scale_fill_brewer(name = "Family", 
                    palette = "Set3", 
                    na.translate = F)
p_Firmicutes2 <- p_Firmicutes1 + new_scale_fill() +
  geom_fruit(
    data = tree_data_ancom_lfc, 
    geom = geom_bar, 
    mapping = aes(y = label_renamed, x = lfc, fill = condition),
    stat = "identity",
    offset = 0.1,
    pwidth = 0.2) +
  scale_fill_manual(name="Suppressive soils",values=c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
p_Firmicutes2
pdf("all_sig_q0.1_fc2/ancom/tree_annot_Firmicutes.pdf", paper="a4")
print(p_Firmicutes2)
dev.off()

samples_to_drop <- tree_data_ancom_phylum$label_renamed[tree_data_ancom_phylum$Phylum!="Actinobacteriota"]
sub_tree <- drop.tip(iqtree.reduced_ancom, samples_to_drop)
tree_data_ancom_Actinobacteriota <- tree_data_ancom_phylum[tree_data_ancom_phylum$Phylum=="Actinobacteriota",]

#sub_tree <- root(sub_tree, outgroup = "NR_102775", edgelabel = TRUE)
genus_Actinobacteriota <- names(sort(table(tree_data_ancom_Actinobacteriota$Genus), decreasing = T))
p_Actinobacteriota <- ggtree(sub_tree,layout= "circular") %<+% tree_data_ancom_Actinobacteriota +
  geom_tippoint(aes(color = Genus, na.translate = FALSE), size = 1) +  # Color tips
  scale_color_manual(name = "Genus",
                     values = cb_palette,
                     na.translate = F)

p_Actinobacteriota1 <- gheatmap(p_Actinobacteriota, tree_data_ancom_Actinobacteriota[4], offset=0.2, width=.1, color = NA, colnames = F) +
  scale_fill_brewer(name = "Family", 
                    palette = "Set3", 
                    na.translate = F)
p_Actinobacteriota2 <- p_Actinobacteriota1 + new_scale_fill() +
  geom_fruit(
    data = tree_data_ancom_lfc, 
    geom = geom_bar, 
    mapping = aes(y = label_renamed, x = lfc, fill = condition),
    stat = "identity",
    offset = 0.1,
    pwidth = 0.2) +
  scale_fill_manual(name="Suppressive soils",values=c("Enriched"="#91BFDB", "Depleted"="#FEE08B"))
p_Actinobacteriota2
pdf("all_sig_q0.1_fc2/ancom/tree_annot_Actinobacteriota.pdf", paper="a4")
print(p_Actinobacteriota2)
dev.off()
