setwd("C:/Users/troop/GitHub/mosquito_virome_yucatan_LEVE")

# Libraries
library(ggcoverage)
library(ggplot2)

#### OSFV #####

png("results/czid/Coverage_Plots_PM2528.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Original genome-read mapping PM2528
PM2528_depth_file <- "data/raw/czid_raw/mapping/OSFV/PM2528_S2_923330/PM2528_S2_923330_depth.txt"
PM2528_before <- read.table(PM2469_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2528_before$stage <- "Before (CZID)"

# Small RNA mapping PM2528
PM2528_smallrna_file <- "results/czid/osfv_complete_genome/PM2528_S2_923330_depth.txt"
PM2528_after <- read.table(PM2528_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2528_after$stage <- "After (small RNA)"

# Combine
combined <- rbind(PM2469_after, PM2469_before)

# Coverage plots for PM2528
ggplot(PM2528_after, aes(x=pos, y=depth, fill=stage)) +
  geom_area(alpha=0.5) +
  facet_wrap(~stage, ncol=1, scales="free_y") +
  labs(x="Genome Position (bp)", y="Coverage Depth (×)",
       title="OSFV PM2528 Coverage After Small RNA Mapping") +
  theme_minimal()

dev.off()


png("results/czid/Coverage_Plots_PM3153.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Original genome-read mapping PM3153
PM3153_depth_file <- "data/raw/czid_raw/mapping/OSFV/PM3153_S5_923338/PM3153_S5_923338_depth.txt"
PM3153_before <- read.table(PM3183_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM3153_before$stage <- "Before (CZID)"

# Small RNA mapping PM3153
PM3153_smallrna_file <- "results/czid/osfv_complete_genome/PM3153_S5_923338_depth.txt"
PM3153_after <- read.table(PM3183_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM3153_after$stage <- "After (small RNA)"

# Combine
combined <- rbind(PM3153_after, PM3153_before)

# Coverage plots for PM3153
ggplot(PM3153_after, aes(x=pos, y=depth, fill=stage)) +
  geom_area(alpha=0.5) +
  facet_wrap(~stage, ncol=1, scales="free_y") +
  labs(x="Genome Position (bp)", y="Coverage Depth (×)",
       title="OSFV PM3153 Coverage After Small RNA Mapping") +
  theme_minimal()

dev.off()
