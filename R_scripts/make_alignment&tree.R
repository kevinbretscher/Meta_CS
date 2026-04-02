##Install packages if not installed
listOfPackages <- c("ape","Biostrings","msa","seqinr","phangorn",
                    "treedataverse","ggnewscale","tidytree","ggpubr",
                    "treeplyr")
for (i in listOfPackages){
  if(! i %in% installed.packages()){
    install.packages(i, dependencies = TRUE)
  }
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("remotes")
BiocManager::install("YuLab-SMU/treedataverse")
remotes::install_github("uyedaj/treeplyr")

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


#read fasta file with sequences that are erniched or depleted
seqs <- readDNAStringSet("Bacillus/Bacillus_new.fasta")
seqs
#rename Silva sequences
names(seqs) <- remove_after_first_dot(names(seqs))
seqs #sanity check for name

#readt representative sequences and outgroup
seqsrep <- readDNAStringSet(c("Bacillus/Bacillus_ref_dsmz.fasta"))
outgroup <- readDNAStringSet(c("Bacillus/outgroup_Paenibacillus_16S.fasta"))

#combine to one faste. Instead of this you could also do cat in de commandline I guess

comb_seqs <- singleseq <- c(seqs,seqsrep,outgroup)
writeXStringSet(comb_seqs, "Bacillus/Bacillus_all.fasta" ,format ="fasta")
readDNAStringSet(c("Bacillus/Bacillus_all.fasta"))
## same for chitino
#read fasta file with sequences that are erniched or depleted
seqs <- readDNAStringSet("Chitino/Chitinophagaceae_new.fasta")
seqs
#rename Silva sequences
names(seqs) <- remove_after_first_dot(names(seqs))
seqs #sanity check for name

seqsrep <- readDNAStringSet(c("Chitino/Chitino_ref_dsmz.fasta"))
outgroup <- readDNAStringSet(c("Chitino/outgroup_Cytophaga_16S"))

#combine to one faste. Instead of this you could also do cat in de commandline I guess

comb_seqs <- singleseq <- c(seqs,seqsrep,outgroup)
writeXStringSet(comb_seqs, "Chitino/Chitino_all.fasta" ,format ="fasta")
#Run following command in comand line

# muscle -align Sphingo_comb_.fasta -output Sphingo_comb_.afa
# Or in case of many sequences (hundreds) muscle -super5 Sphingo_comb_.fasta -output Sphingo_comb.afa\
#then run iqtree, make sure to give the outgroup (! the taxa id in the header, not the file)
# iqtree2 -s Sphingo_comb_with_all.afa -B 1000 -alrt 1000 -nt auto -o Novosphingobium_capsulatum__GIFU_11526__D16147  


library(treedataverse)
library(ggnewscale)
library(tidytree)
library(treeplyr)
library(ggtreeExtra)
library(ggsci)
#install.packages("ggstar")
library(ggstar)
library(treeio)

iqtree <- read.iqtree("Bacillus/Bacillus_all.afa.treefile")
iqtree_chitino <- read.iqtree("Chitino/Chitino_all.afa.treefile")

iqtree.all <- read.iqtree("Sphingo/Sphingo_comb_with_all.afa.treefile")
ggtree(iqtree, layout='fan')
ggtree(iqtree_chitino, layout='fan')

iqtree <- root(iqtree, outgroup = "AJ320493.1", edgelabel = TRUE)


to_drop <- c("AJLS01000166",
             "U20384")
to_drop_chitino <- c("FPLS01052098","FPLS01051226","FPLS01006454", "FPLS01049069")


iqtree.reduced <- drop.tip(iqtree, to_drop)
iqtree.reduced <- root(iqtree.reduced, outgroup = "AJ320493.1", edgelabel = TRUE)

iqtree.reduced_chitino <- drop.tip(iqtree_chitino, to_drop_chitino)
iqtree.reduced_chitino <- root(iqtree.reduced_chitino, outgroup = "M58768.1", edgelabel = TRUE)
ggtree(iqtree.reduced_chitino)
tree_data_chitino <- iqtree.reduced_chitino %>% fortify()
tree_data_chitino$group <- "Other"
tree_data_chitino$group[tree_data_chitino$label %in% annotation$label] <- "Enriched"
tree_data_chitino$group[tree_data_chitino$label %in% annotation_depleted$label] <- "Depleted"

ggtree(iqtree.reduced_chitino, layout= "circular") %<+% tree_data_chitino +  # Attach metadata to the tree
  geom_tiplab(aes(label = label, color = group), show.legend = FALSE, size= 1, align=TRUE, linesize=.5) +  # Color labels
  geom_tippoint(aes(color = group), size = 1) +  # Color tips
  scale_color_manual(values = c("Enriched" = "red","Depleted" = "blue", "Other" = "black")) 

ggtree(iqtree.reduced_chitino) %<+% tree_data_chitino +  # Attach metadata to the tree
  geom_tiplab(aes(label = label, color = group), show.legend = FALSE, size= 2, align=TRUE, linesize=.5) +  # Color labels
  geom_tippoint(aes(color = group), size = 1) +  # Color tips
  scale_color_manual(values = c("Enriched" = "red","Depleted" = "blue", "Other" = "black")) 

# ggtree(iqtree.reduced, layout= "circular")+#, aes(color = label %in% rownames(heatmap_matrix))) +
#   geom_tiplab(size=1, align=TRUE, linesize=.5) +
#   geom_point(aes(subset = label %in% annotation$label), color = "red", size = 3)

# Convert tree into a dataframe for joining
tree_data_bacillus <- iqtree.reduced %>% 
  fortify()

tree_data_bacillus$group <- "Other"
tree_data_bacillus$group[tree_data_bacillus$label %in% annotation$label] <- "Enriched"
tree_data_bacillus$group[tree_data_bacillus$label %in% annotation_depleted$label] <- "Depleted"

# Plot with highlighted tips and labels
ggtree(iqtree.reduced, layout= "circular") %<+% tree_data_bacillus +  # Attach metadata to the tree
  geom_tiplab(aes(label = label, color = group), show.legend = FALSE, size=2, align=TRUE, linesize=.5) +  # Color labels
  geom_tippoint(aes(color = group), size = 1) +  # Color tips
  scale_color_manual(values = c("Enriched" = "red","Depleted" = "blue", "Other" = "black")) 

sort(table(enriched_OTU_ANCOM$Genus[enriched_OTU_ANCOM$Family=="Chitinophagaceae"]))
sort(table(depleted_OTU_ANCOM$Genus[depleted_OTU_ANCOM$Family=="Chitinophagaceae"]))

setdiff(enriched_OTU_ANCOM$Genus[enriched_OTU_ANCOM$Family=="Chitinophagaceae"],
        depleted_OTU_ANCOM$Genus[depleted_OTU_ANCOM$Family=="Chitinophagaceae"])
setdiff(depleted_OTU_ANCOM$Genus[depleted_OTU_ANCOM$Family=="Chitinophagaceae"],
        enriched_OTU_ANCOM$Genus[enriched_OTU_ANCOM$Family=="Chitinophagaceae"])

sort(table(enriched_OTU_ANCOM$study[enriched_OTU_ANCOM$Genus=="Sphingomonas"]))
sort(table(enriched_OTU_ANCOM$study[enriched_OTU_ANCOM$Genus=="Flavobacterium"]))
sort(table(depleted_OTU_ANCOM$study[depleted_OTU_ANCOM$Genus=="Flavobacterium"]))
sort(table(enriched_OTU_ANCOM$study[enriched_OTU_ANCOM$Genus=="Chitinophaga"]))
sort(table(depleted_OTU_ANCOM$study[depleted_OTU_ANCOM$Genus=="Chitinophaga"]))
#read enriched otus from Download 16S script and moduft label so that it matches the tree

annotation <- enriched_OTU_ANCOM %>%  group_by(taxon) %>% 
  summarize(studies_up = paste0(study, collapse = "|")) %>%
  mutate(Taxa = remove_after_first_dot(taxon)) %>%
  rename(label = Taxa)
annotation$sig <- "enriched"

annotation_depleted <- depleted_OTU_ANCOM %>%  group_by(taxon) %>% 
  summarize(studies_down = paste0(study, collapse = "|")) %>%
  mutate(Taxa = remove_after_first_dot(taxon)) %>%
  rename(label = Taxa)
annotation_depleted$sig <- "depleted"


#merge the data with our tree

bp2.df.meta <- left_join(x=iqtree.reduced, y=annotation, by="label")
bp2.df.meta <- left_join(x=bp2.df.meta, y=annotation_depleted, by="label")

# #read lifestyle doc and rename because the f*king names never match
# 
# lifestyle.i <- read.delim("metadata_annotated.txt")
# lifestyle.i$label <- paste(lifestyle.i$Assembly.Accession,"_O", sep = "")
# row.names(lifestyle.i) <- lifestyle.i$label
# 
# #rename and shorten genome labels so they match the lifestyle file i got
# 
# extract_gcf_text <- function(input_strings) {
#   sapply(input_strings, function(input_string) {
#     # Find the matching pattern
#     match <- regmatches(input_string, regexpr("GCF.*?_O", input_string))
#     
#     # Check if a match was found; return it if so, else return the original string
#     if (length(match) > 0) {
#       return(match)
#     } else {
#       return(input_string)
#     }
#   })
# }
# 
# 
# 
# bp2.df.meta@phylo$tip.label <- extract_gcf_text(bp2.df.meta@phylo$tip.label)
# 
# #drop some more outliers
# 
# to_drop <- c("GCF.036855915.1_O",
#              "GCF.022601895.1_O",
#              "GCF.000732685.2_O")
# 
# bp2.df.meta <- drop.tip(bp2.df.meta, to_drop)
# 
# #merge with the new data
# 
# bp2.df.meta <- left_join(x=bp2.df.meta, y=lifestyle.i, by="label")
# #bp2.df.meta <- left_join(x=bp2.df.meta, y=study_wide_heatmap, by="label")

##Make a tree with lifestyle
# 
# tree.fan <- ggtree(bp2.df.meta, layout='fan') +
#   geom_tiplab(align=TRUE, size = 0.5) 
# 
# tree.fan.bp <- tree.fan + new_scale_fill() + new_scale_color() + geom_point(aes(fill = UFboot, color = UFboot), shape = 23, size = 1, na.rm = TRUE) + 
#   scale_fill_gsea(name="UFboot", na.value=NA) + scale_color_gsea(name="UFboot", na.value=NA)
# 
# p.life <- tree.fan.bp + new_scale_fill() + new_scale_color() +
#   geom_fruit(geom=geom_tile, mapping=aes(y = label, fill = Lifestyle),  width=0.01 , offset = 0.5, size = 1) +
#   scale_fill_discrete(na.value = 'white') 


# alternative tree with studies 
p.final <- tree.fan.bp + new_scale_fill() + new_scale_color() +
  geom_fruit(geom=geom_point, mapping=aes(y=label, color = studies_down, size=studies_down), offset = 0.5, size = 1) +
  scale_color_discrete(na.value = 'white') 

p.final2 <- p.final + new_scale_fill() + new_scale_color() +
  geom_fruit(geom=geom_point, mapping=aes(y=label, color = studies_up, size=studies_up), offset = 0.1, size = 1) +
  scale_color_discrete(na.value = 'white')

