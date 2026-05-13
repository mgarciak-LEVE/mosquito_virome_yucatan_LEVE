setwd("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/scripts/czid_analysis")

# Load the libraries and dataset
library(vegan)
library(labdsv)
library(ggplot2)
library(dplyr)
library(barrel)
library(cooccure)
library(igraph)
library(reshape2)

## DATASET LOADING -------

dataset_main <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/df_final.csv")
dataset_metadata <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/metadata.csv")
dataset_diversity <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/diversity_metrics.csv")

## MAIN DATASET MANIPULATION ----

# Filter to keep only virus rows (where category == "viruses")
virus_data <- dataset_main[dataset_main$category == "viruses", ]

# Abundance matrix
abundance_matrix <- reshape2::dcast(virus_data, 
                                    sample_name ~ name, 
                                    value.var = "nt_rpm", 
                                    fun.aggregate = sum, 
                                    fill = 0)

abundance_matrix[is.na(abundance_matrix)] <- 0

# Set sample names as row names
rownames(abundance_matrix) <- abundance_matrix$sample_name
abundance_matrix$sample_name <- NULL

## DIVERSITY ANALYSES -------



# RAREFACTION CURVES ----


# NMDS ----

dist_matrix <- vegdist(abundance_matrix, method = "bray")
set.seed(123)  # for reproducibility
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)

nmds_result$stress

# View convergence details
nmds_result

# Extract NMDS scores for samples
nmds_scores <- as.data.frame(nmds_result$points)
colnames(nmds_scores) <- c("NMDS1", "NMDS2")
nmds_scores$sample <- rownames(nmds_scores)

# Add metadata to the scores
nmds_scores <- merge(nmds_scores, dataset_metadata, by.x = "sample", by.y = "sample_name")

# Create plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = habitat)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = habitat), linetype = "dashed") +
  theme_minimal() +
  labs(title = "NMDS of Viral Communities",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       x = "NMDS1", y = "NMDS2")


