setwd("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE")

# Load the libraries and dataset
library(vegan)
library(labdsv)
library(ggplot2)
library(dplyr)
library(barrel)
library(reshape2)
library(tidydr)
library(iNEXT)
library(fossil)
library(tidyr)
library(gridExtra)

## DATASET LOADING -------

dataset_main <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/df_final.csv")
dataset_metadata <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/metadata.csv")
dataset_diversity <- read.csv("/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/docs/czid_sequence_reports/processed/diversity_metrics.csv")

## MAIN DATASET MANIPULATION ----

# Filter to keep only virus rows (where category == "viruses")
virus_data <- dataset_main[dataset_main$category == "viruses", ]

## ABUNDANCE MATRIX ----


abundance_matrix <- reshape2::dcast(virus_data, 
                                    sample_name ~ name, 
                                    value.var = "nt_rpm", 
                                    fun.aggregate = sum, 
                                    fill = 0)

# NAs into zeros
abundance_matrix[is.na(abundance_matrix)] <- 0

# Set sample names as row names
rownames(abundance_matrix) <- abundance_matrix$sample_name
abundance_matrix$sample_name <- NULL

# Remove habitat rows with zero abundance
abundance_matrix <- abundance_matrix[rowSums(abundance_matrix) > 0, ]
# Remove species columns with zero abundance
abundance_matrix <- abundance_matrix[, colSums(abundance_matrix) > 0, ]


# A. DESCRIPTIVE PLOTS ----

# A.1. VENN DIAGRAM ----

png("results/czid/virome_venn_diagram.png",
    width = 10, height = 8, units = "in", res = 600, bg = "white")

# Create presence/absence by habitat
viral_presence <- virus_data %>%
  dplyr::select(sample_name, name, habitat) %>%
  distinct() %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = habitat, values_from = present, 
                     values_fill = 0, values_fn = max)

# Get viral lists per habitat
habitat_viruses <- list()
for(h in unique(virus_data$habitat)) {
  habitat_viruses[[h]] <- virus_data %>%
    filter(habitat == h) %>%
    pull(name) %>%
    unique()
}


# Create Venn diagram (if you have 2-4 habitats)

library(VennDiagram)

venn.plot <- venn.diagram(
  x = habitat_viruses,
  category.names = names(habitat_viruses),
  filename = NULL,
  output = TRUE,
  # Circles
  lwd = 1,
  lty = 'solid',
  fill = RColorBrewer::brewer.pal(length(habitat_viruses), "Dark2"),
  col = "black",
  # Make circles more circular
  rotation.degree = 0,        # Rotation angle
  scaled = TRUE,              # Scales circles to match overlap sizes
  euler.d = TRUE,
  # Numbers
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.0,
  cat.default.pos = "outer",
  cat.fontface = "bold"
  # Output features
)

grid.newpage()
grid.draw(venn.plot)
dev.off()

## B. DIVERSITY ANALYSES -------

# Calculate diversity indices
diversity <- data.frame(
  sample_name = rownames(abundance_matrix),
  richness = specnumber(abundance_matrix), # Number of detected viruses
  shannon = vegan::diversity(abundance_matrix),           # Default is Shannon
  simpson = vegan::diversity(abundance_matrix, "simpson"),
  invsimpson = vegan::diversity(abundance_matrix, "invsimpson"),
  chao1 = chao1(abundance_matrix, taxa.row = TRUE),
  
  total_rpm = rowSums(abundance_matrix)
)

# Metadata
diversity <- diversity %>%
  left_join(dataset_metadata, by = "sample_name")

# B.1. Alpha Diversity ----

# A.1.1 PLOTS ----

## A.1.1 Boxplots By Habitat ----

png("results/czid/alpha_diversity_habitat.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Shannon diversity
plot1 <- ggplot(diversity, aes(x = habitat, y = shannon, color = habitat)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Shannon Diversity by Habitat",
                      x = "Habitat", y = "Shannon Diversity Index")

# Simpson index
plot2 <- ggplot(diversity, aes(x = habitat, y = simpson, color = habitat)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Simpson Index by Habitat",
       x = "Habitat", y = "Simpson Index")

# Species richness
plot3 <- ggplot(diversity, aes(x = habitat, y = richness, color = habitat)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Species Richness by Habitat",
       x = "Habitat", y = "Species Richness")


grid.arrange(plot1, plot2, plot3, ncol = 3)
dev.off()


# A.1.2 Boxplots By Species ----

# Shannon diversity
png("results/czid/alpha_diversity_species.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

plot4 <- ggplot(diversity, aes(x = species, y = shannon, fill = species)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Shannon Diversity by Mosquito Species",
       x = "Species", y = "Shannon Diversity Index")

# Simpson index
plot5 <- ggplot(diversity, aes(x = species, y = simpson, fill = species)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Simpson Index by Mosquito Species",
       x = "Species", y = "Simpson Index")

# Species richness
plot6 <- ggplot(diversity, aes(x = species, y = richness, fill = species)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(position = position_jitter(width = 0.2)) +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Species Richness by Mosquito Species",
       x = "Species", y = "Species Richness")

grid.arrange(plot4, plot5, plot6, ncol = 3)
dev.off()

# B.2 Beta Diversity ----

## Bray-Curtis

png("results/czid/nmds_plot.png", 
    width = 10, height = 8, units = "in", res = 600, bg = "white")

# Dissimilarity calculation
set.seed(123)
braycurtis_dist <- vegdist(abundance_matrix, method = "bray", binary = T)

# B.2.1 NMDS ----

# Run NMDS
nmds <- metaMDS(braycurtis_dist, k = 2, trymax = 100)

# Extract scores
nmds_scores <- as.data.frame(vegan::scores(nmds, display = "sites"))
nmds_scores <- cbind(nmds_scores, diversity)

# Calculate centroid for each habitat
nmds_centroids <- nmds_scores %>%
  group_by(habitat) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = 'drop')

# Plot NMDS

nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = habitat, shape = species), size = 4, alpha = 0.8) +
  geom_text(aes(label = sample_name), hjust = -0.3, vjust = -0.3, size = 2.5) +
  stat_ellipse(aes(color = habitat), level = 0.95, alpha = 0.2) +
  geom_point(data = nmds_centroids, aes(x = NMDS1, y = NMDS2, fill = habitat),
             size = 5, shape = 23, stroke = 1.5) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Viral Community Composition (NMDS)",
       subtitle = paste("Stress =", round(nmds$stress, 3), 
                        "| Diamonds = habitat centroids"),
       color = "Habitat", shape = "Species", fill = "Habitat Centroid")

print(nmds_plot)
dev.off

# B.2.2.PCOA ----

png("results/czid/pcoa_plot.png", 
    width = 10, height = 8, units = "in", res = 600, bg = "white")

pcoa_result <- cmdscale(braycurtis_dist, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores <- cbind(pcoa_scores, diversity)

# Calculate variance explained
pcoa_var <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100, 1)

pcoa_plot <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = habitat, shape = species), size = 4, alpha = 0.8) +
  stat_ellipse(aes(color = habitat), level = 0.95, alpha = 0.2) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Principal Coordinates Analysis (Bray-Curtis)",
       subtitle = paste0("PCoA1: ", pcoa_var[1], "% | PCoA2: ", pcoa_var[2], "%"),
       x = paste0("PCoA1 (", pcoa_var[1], "%)"),
       y = paste0("PCoA2 (", pcoa_var[2], "%)"),
       color = "Habitat", shape = "Species")

print(pcoa_plot)
dev.off()

