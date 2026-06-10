setwd("C:/Users/troop/GitHub/mosquito_virome_yucatan_LEVE")

# Libraries
library(ggcoverage)
library(ggplot2)

#### OSFV #####

png("results/czid/Coverage_Plots_PM2528.png", 
    width = 18, height = 6, units = "in", res = 600, bg = "white")

# Original genome-read mapping PM2528
PM2528_depth_file <- "data/raw/czid_raw/mapping/OSFV/PM2528_S2_923330/PM2528_S2_923330_depth.txt"
PM2528_before <- read.table(PM2528_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2528_before$stage <- "Before (CZID)"

# Small RNA mapping PM2528
PM2528_smallrna_file <- "results/czid/osfv_complete_genome/PM2528_S2_923330_depth.txt"
PM2528_after <- read.table(PM2528_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM2528_after$stage <- "After (small RNA)"

# Combine
combined <- rbind(PM2528_before, PM2528_after)

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
PM3153_before <- read.table(PM3153_depth_file, header=FALSE, col.names=c("contig", "pos", "depth"))
PM3153_before$stage <- "Before (CZID)"

# Small RNA mapping PM3153
PM3153_smallrna_file <- "results/czid/osfv_complete_genome/PM3153_S5_923338_depth.txt"
PM3153_after <- read.table(PM3153_smallrna_file, header=FALSE, col.names=c("contig", "pos", "depth"))
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



##############################################################

# Small RNA mapping to PM3153 consensus genome

####========================================####
####  PM2528: Three Coverage Panels         ####
####========================================####

png("results/czid/Coverage_Plots_PM2528_extra.png", 
    width = 18, height = 8, units = "in", res = 600, bg = "white")

# Small RNA mapping (PM2528 reads → PM2528 consensus)
PM2528_own <- read.table("results/czid/osfv_complete_genome/PM2528_S2_923330_depth.txt",
                         col.names=c("contig","pos","depth"))
PM2528_own$stage <- "Small RNA mapped to PM2528 consensus"

# Cross-mapping (PM2528 reads → PM3153 consensus)
PM2528_cross <- read.table("results/czid/osfv_complete_genome/PM2528_on_PM3153_depth.txt",
                           col.names=c("contig","pos","depth"))
PM2528_cross$stage <- "Small RNA cross-mapped to PM3153 consensus"

# Combine
combined_pm2528 <- rbind(PM2528_own, PM2528_cross)

ggplot(combined_pm2528, aes(x=pos, y=depth, fill=stage)) +
  geom_area(alpha=0.6) +
  facet_wrap(~stage, ncol=1, scales="free_y") +
  scale_fill_manual(values=c("Small RNA mapped to PM2528 consensus"="steelblue",
                             "Small RNA cross-mapped to PM3153 consensus"="darkgreen")) +
  labs(x="Genome Position (bp)", y="Coverage Depth (×)",
       title="OSFV PM2528: Small RNA Coverage",
       subtitle="Own consensus vs. Cross-mapped to PM3153 reference genome") +
  theme_minimal() +
  theme(strip.text=element_text(size=10, face="bold"),
        plot.title=element_text(size=14, face="bold"))

dev.off()

