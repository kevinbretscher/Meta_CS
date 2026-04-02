iqtree.all <- read.iqtree("Sphingo/Sphingo_comb_with_all.afa.treefile")
iqtree.genome_sphingo <- read.iqtree("../baclife_sphingo/baclife/workflow/intermediate_files/phylophlan/output_phylophlan/input_resolved.tre")

metadata_sphingo <- read.delim("Sphingo/metadata_for_baclife.txt")
baclife_names <- read.delim("../baclife_sphingo/baclife/workflow/Shiny_app/input/names_equivalence.txt", sep = " ")

library(dplyr)
library(stringr)

# Extract accession ID from Full_name
baclife_names_1 <- baclife_names %>%
  mutate(Accession_ID = str_extract(Full_name, "GCF\\.\\d+\\.\\d+"))

metadata_sphingo <- metadata_sphingo %>%
  mutate(Accession_ID = gsub("_", ".", Assembly.Accession))

baclife_names$Accession_ID %in% metadata_sphingo$Accession_ID

merged_baclife_mapping <- left_join(baclife_names_1, metadata_sphingo, 
                       by =  "Accession_ID")

merged_baclife_mapping[is.na(merged_baclife_mapping$Tree),]

ggtree(iqtree.genome_sphingo, layout='fan')
tip.label(iqtree.genome_sphingo)  <- gsub("'", "", tip.label(iqtree.genome_sphingo))

merged_baclife_mapping$tip.label <- gsub(".fna", "", merged_baclife_mapping$bacLIFE_name)
tip.label(iqtree.genome_sphingo) %in% merged_baclife_mapping$tip.label
highlighted_labels <- merged_baclife_mapping$tip.label[merged_baclife_mapping$Tree=="supp"]


tree_data <- data.frame(label = tip.label(iqtree.genome_sphingo),
                        group = ifelse(ip.label(iqtree.genome_sphingo)%in% highlighted_labels, "Highlighted", "Other"))


ggtree(iqtree.genome_sphingo)  %<+% tree_data +
  geom_tiplab(aes(label = label, color = group),size = 1, offset = 0.5, show.legend = FALSE) +  # Color labels
  geom_tippoint(aes(color = group), size = 0.1) +  # Color tips
  scale_color_manual(values = c("Highlighted" = "red", "Other" = "black"))  # Define colors

                    