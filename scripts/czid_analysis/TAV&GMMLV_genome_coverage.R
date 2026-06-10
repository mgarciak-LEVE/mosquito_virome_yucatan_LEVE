setwd("C:/Users/troop/GitHub/mosquito_virome_yucatan_LEVE")

# Libraries
library(ggcoverage)
library(ggplot2)

#### TAV #####

png("results/czid/Coverage_Plots_PM2469.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Original genome-read mapping PM2528
PM2469_depth_file <- "data/raw/czid_raw/mapping/TAV/PM2469_S1_923325/PM2469_S1_923325_depth.txt"
PM2469_before <- read.table(PM2469_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2469_before$stage <- "Before (CZID)"

# Small RNA mapping PM2528
PM2469_smallrna_file <- "results/czid/tav_complete_genome/PM2469_S1_923325_depth.txt"
PM2469_after <- read.table(PM2469_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2469_after$stage <- "After (small RNA)"

# Combine
combined <- rbind(PM2469_before, PM2469_after)

# Coverage plots for PM2528
ggplot(combined, aes(x=pos, y=depth, fill=stage)) +
  geom_area(alpha=0.5) +
  facet_wrap(~stage, ncol=1, scales="free_y") +
  labs(x="Genome Position (bp)", y="Coverage Depth (×)",
       title="TAV PM2469 Coverage After Small RNA Mapping") +
  theme_minimal()

dev.off()

############################################################################################

#### GMMLV #####

png("results/czid/Coverage_Plots_PM2622_GMMLV.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Original genome-read mapping PM3183
PM2622_depth_file <- "data/raw/czid_raw/mapping/GMMLV/PM2622_S32_923335/PM2622_S32_923335_depth.txt"
PM2622_before <- read.table(PM2622_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2622_before$stage <- "Before (CZID)"

# Small RNA mapping PM3183
PM2622_smallrna_file <- "results/czid/gmmlv_complete_genome/PM2622_S32_923335_depth.txt"
PM2622_after <- read.table(PM2622_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2622_after$stage <- "After (small RNA)"

# Combine
combined <- rbind(PM2622_after, PM2622_before)

# Coverage plots for PM3183
ggplot(combined, aes(x=pos, y=depth, fill=stage)) +
  geom_area(alpha=0.5) +
  facet_wrap(~stage, ncol=1, scales="free_y") +
  labs(x="Genome Position (bp)", y="Coverage Depth (×)",
       title="GMMLV PM2622 Coverage After Small RNA Mapping") +
  theme_minimal()

dev.off()

