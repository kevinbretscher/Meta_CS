library(dplyr)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")
library(phyloseq)
#library(speedyseq)
library(vegan)
library(ape)
install.packages("ggdist")
library(ggdist)

ps.import <- readRDS("../Scripts/phyloseq.RDS")
ps.import <- subset_samples(ps.import, B_R=="R")

metadata <- as.matrix(sample_data(ps.import))
metadata <- as.data.frame.matrix(metadata)
metadata$SampleID <- rownames(metadata)


max.runs <- 100
min_lib <- min(sample_sums(ps.import))
nsamp = nsamples(ps.import)


observed <- matrix(nrow = nsamp, ncol = max.runs)
row.names(observed) <- sample_names(ps.import)

shannon <- matrix(nrow = nsamp, ncol = max.runs)
row.names(shannon) <- sample_names(ps.import)

set.seed(1998)



for (i in 1:max.runs) {
  # Sub sample
  r <- rarefy_even_depth(ps.import, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  obs <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  observed[ ,i] <- obs
  
  # Calculate evenness
  shann <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
  shannon[ ,i] <- shann
}


# Create a new dataframe to hold the means and standard deviations of obs estimates
SampleID <- row.names(observed)
mean <- apply(observed, 1, mean)
sd <- apply(observed, 1, sd)
measure <- rep("observed", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of shannon estimates
SampleID <- row.names(shannon)
mean <- apply(shannon, 1, mean)
sd <- apply(shannon, 1, sd)
measure <- rep("Shannon", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

#merge with metadata
alpha <- rbind(rich_stats, even_stats)
alpha <- merge(alpha, metadata, by = "SampleID") 

# Reorder factor levels alphabetically
alpha$C_S <- factor(alpha$C_S, levels = sort(levels(factor(alpha$C_S))))

#visualize

library(ggplot2)
library(ggpubr)
library(ggsci)
library(gghalves)

my_comparisons <- list( c("C", "S"))

obs_plot <- ggviolin(alpha[alpha$measure == 'observed', ], x = "C_S", y = "mean", fill = "C_S", size = 0.8 ,draw_quantiles = 0.5, title = "Observed ASVs", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" , label.y = 2200) + scale_fill_npg() + rotate_x_text(angle = 45) +
  theme(
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1)
  ) 

shann_plot <- ggviolin(alpha[alpha$measure == 'Shannon', ], x = "C_S", y = "mean", fill = "C_S", size = 0.8 ,draw_quantiles = 0.5, title = "Shannon", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" , label.y = 7.1) + scale_fill_npg() + rotate_x_text(angle = 45) +
  theme(
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1)
  )

ggarrange(obs_plot,shann_plot, common.legend = TRUE)


ggviolin(alpha, x = "C_S", y = "mean", fill = "C_S", size = 0.8 ,draw_quantiles = 0.5, title = "Alpha diversity", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" ) + rotate_x_text(angle = 45) +
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  facet_wrap(~measure, scales = "free_y")+
  scale_x_discrete(labels = c("S" = "", "C"="")) +
  theme(
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1)
  ) +
  theme_pubclean() + theme_pubr(border=T)

obs_plot <- ggplot(alpha[alpha$measure == 'observed', ], aes(x = C_S, y = mean, fill = C_S)) +
  ggdist::stat_halfeye(
    adjust = 1,        # Increase smoothing for density
    width = 0.6,       # Width of density
    justification = -0.3, # Offset density plot
    point_colour = NA   # Remove points in density
  ) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
  geom_jitter(
    aes(color = C_S),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Observed",
    x = "Group",
    y = "Value"
  )+
  theme_pubclean() + theme_pubr(border=T)

shann_plot <- ggplot(alpha[alpha$measure == 'Shannon', ], aes(x = C_S, y = mean, fill = C_S)) +
  ggdist::stat_halfeye(
    adjust = 1,        # Increase smoothing for density
    width = 0.6,       # Width of density
    justification = -0.3, # Offset density plot
    point_colour = NA   # Remove points in density
  ) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
  geom_jitter(
    aes(color = C_S),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +
  scale_fill_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Shannon",
    x = "Group",
    y = "Value"
  )+
  theme_pubclean() + theme_pubr(border=T)

p_alpha_combined <- ggplot(alpha, aes(x = C_S, y = mean)) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
geom_jitter(
    aes(color = C_S),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +  
  facet_wrap(~measure, scales = "free_y", labeller = labeller(measure = c(observed = "Observed OTUs", Shannon = "Shannon Index")))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"), name = "", labels = c(C="Conducive", S="Suppressive"))+
  #theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" ) +
  labs(
    title = "",
    x = "",
    y = ""
  )+
  scale_x_discrete(labels = c(C="Conducive", S="Suppressive"))+
  theme_pubclean() + theme_pubr(border=T) +   theme(legend.position = "right") 
ggsave("../../map/alpha_combined.pdf", plot = p_alpha_combined, width = 7, height = 4, units = "in", dpi = 300)


ggplot(alpha, aes(x = C_S, y = mean)) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
  geom_jitter(
    aes(color = C_S, shape = Study),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +  
  facet_wrap(~measure, scales = "free_y", labeller = labeller(measure = c(observed = "Observed OTUs", Shannon = "Shannon Index")))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  scale_shape_manual(values = 0:20) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" ) +
  labs(
    title = "Alpha Diversity",
    x = "Group",
    y = "Value"
  )+
  scale_x_discrete(labels = c(C="Conducive", S="Suppressive"))+
  theme_pubclean() + theme_pubr(border=T)


ggplot(alpha[!alpha$Study %in% c("Ruth-zwaagdijk", "PRJNA578725"),], aes(x = C_S, y = mean)) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
  geom_jitter(
    aes(color = C_S),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +  
  facet_wrap(~measure, scales = "free_y", labeller = labeller(measure = c(observed = "Observed OTUs", Shannon = "Shannon Index")))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" ) +
  labs(
    title = "Alpha Diversity",
    x = "Group",
    y = "Value"
  )+
  scale_x_discrete(labels = c(C="Conducive", S="Suppressive"))+
  theme_pubclean() + theme_pubr(border=T)

#per study
ggplot(alpha, aes(x = C_S, y = mean)) +
  geom_boxplot(
    width = 0.15,      # Box width
    outlier.shape = NA, # Remove outliers
    alpha = 0.5        # Transparency
  ) +
  geom_jitter(
    aes(color = C_S),
    width = 0.1,       # Jitter width
    alpha = 0.6        # Transparency
  ) +  
  facet_grid(measure~Study, scales = "free_y", labeller = labeller(measure = c(observed = "Observed OTUs", Shannon = "Shannon Index")))+
  scale_color_manual(values = c(S="#91BFDB", C ="#FEE08B"))+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif" ) +
  labs(
    title = "Alpha Diversity",
    x = "Group",
    y = "Value"
  )+
  theme_pubclean() + theme_pubr(border=T) +
  theme(
        axis.title.x = element_blank(),   # Remove x-axis label
        axis.text.x = element_blank(),    # Remove x-axis text
        axis.ticks.x = element_blank())    # Remove x-axis ticks) +
