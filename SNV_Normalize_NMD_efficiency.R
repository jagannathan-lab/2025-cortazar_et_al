# Title: Normalization of NMD efficiency values
# Description: Calculates SNV-normalized NMD efficiency values for PTC and SNV variants
#Example shown for how to plot NMD efficiency values for PTC

library(tidyverse)


#Input files with replicate values
data_PTC <- "NMD_efficiency_values_PTC.txt"
data_SNV <- "NMD_efficiency_values_SNV.txt"

#The NMD matrix table contains the header
dfx_data_PTC <- read.table(data_PTC, header = TRUE)

#colnames for the SNV replicates
colnames = c("codon", "mutation", "NMD")
dfx_SNV <- read.table(data_SNV, col.names=colnames)

#Range of codon positions of SNVs in each SGE experiment
#Exon (X), Region (R)
#Exon 1 SGE Experiments codons 1 - 118
#X1R1 1 - 27
#X1R2 28 - 58
#X1R3 59 - 88
#X1R4 89 - 118
#Exon 4 SGE Experiments codons 214 - 270
#X4R1 214 - 246
#X4R2 247 - 270
#X7R3 434 - 460
#X10R1 537 - 566
#Exon 11 567 - 656

#Calculate the median of SNVs in each individual SGE experiment
SNVs_median_X1R1 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon <= 27) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))

#for the range of PTCs from each SGE region, reported in the merge file,
#store the corresponding SNV NMD median value
dfx_SNVs_median_X1R1 <- tibble(codon = 1:23) %>%
  mutate(SNVmedian_NMD = SNVs_median_X1R1$SNVmedian_NMD)

SNVs_median_X1R2 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 28 & codon <= 58) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X1R2 <- tibble(codon = 24:58) %>%
  mutate(SNVmedian_NMD = SNVs_median_X1R2$SNVmedian_NMD)

SNVs_median_X1R3 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 59 & codon <= 88) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X1R3 <- tibble(codon = 59:86) %>%
  mutate(SNVmedian_NMD = SNVs_median_X1R3$SNVmedian_NMD)

SNVs_median_X1R4 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 89 & codon <= 118) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X1R4 <- tibble(codon = 87:118) %>%
  mutate(SNVmedian_NMD = SNVs_median_X1R4$SNVmedian_NMD)
        

SNVs_median_X4R1 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 214 & codon <= 246) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X4R1 <- tibble(codon = 214:239) %>%
  mutate(SNVmedian_NMD = SNVs_median_X4R1$SNVmedian_NMD)

SNVs_median_X4R2 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 247 & codon <= 270) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X4R2 <- tibble(codon = 240:270) %>%
  mutate(SNVmedian_NMD = SNVs_median_X4R2$SNVmedian_NMD)

SNVs_median_X7R3 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 434 & codon <= 460) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X7R3 <- tibble(codon = 434:460) %>%
  mutate(SNVmedian_NMD = SNVs_median_X7R3$SNVmedian_NMD)

SNVs_median_X10R1 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 537 & codon <= 566) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X10R1 <- tibble(codon = 537:566) %>%
  mutate(SNVmedian_NMD = SNVs_median_X10R1$SNVmedian_NMD)

SNVs_median_X11R3 <- dfx_SNV %>% 
    select(codon, mutation, NMD) %>%
    mutate(codon = as.numeric(as.character(codon)), NMD = as.numeric(NMD)) %>%
    filter(codon >= 567 & codon <= 656) %>%
    group_by(codon, mutation) %>%
    summarize(SNVexp_NMD = mean(NMD)) %>%
    ungroup() %>%
    summarize(SNVmedian_NMD = median(SNVexp_NMD))
dfx_SNVs_median_X11R3 <- tibble(codon = 567:656) %>%
  mutate(SNVmedian_NMD = SNVs_median_X11R3$SNVmedian_NMD)
        

##Apend all dfx_SNVs_median values into a single data frame by stacking rows 
dfx_SNVs_medians <- bind_rows(dfx_SNVs_median_X1R1, dfx_SNVs_median_X1R2, dfx_SNVs_median_X1R3,
        dfx_SNVs_median_X1R4, dfx_SNVs_median_X4R1, dfx_SNVs_median_X4R2, dfx_SNVs_median_X7R3,
        dfx_SNVs_median_X10R1, dfx_SNVs_median_X11R3)


######### Prepare subset data
PTC_data <- dfx_data_PTC %>%
  select(codon, mutation, mean_val, std_val) %>%
  mutate(codon = as.numeric(as.character(codon)), mean_val = as.numeric(mean_val), std_val = as.numeric(std_val)) %>% 
  filter(std_val <= 0.2)

#Adding a new column for the median NMD value from SNVs and normalizing to SNV median
DMSO_PTC_data_and_medians <- PTC_data %>%
    left_join(dfx_SNVs_medians, by = "codon") %>%
    mutate(norm_mean = mean_val - SNVmedian_NMD)

#Create unique codon mapping for evenly spaced x positions
#Create codon map with gaps between exons
#Create codon map with gaps between exons
codon_map <- dfx_data_PTC %>%
  distinct(codon) %>%
  arrange(codon) %>%
  mutate(
    #Create a numeric x position with artificial gaps
    row = case_when(
      codon <= 118 ~ row_number(),                     #Exon 1
      codon > 118 & codon <= 270 ~ row_number() + 3,  #+5 gap after exon 1
      codon > 270 & codon <= 460 ~ row_number() + 6,  #+5 more after exon 4
      codon > 460 & codon <= 566 ~ row_number() + 9,  #+5 more after exon 7
      codon > 566 & codon <= 656 ~ row_number() + 13,
      TRUE ~ row_number() + 13
    )
  )

#Apending the row ID for plotting continuous positions
#Join back to DMSO_PTC_data_and_medians
DMSO_PTC_data_and_medians <- DMSO_PTC_data_and_medians %>%
  left_join(codon_map, by = "codon")

#Tick labels: every 5 unique codons
tick_labels_df <- codon_map %>%
  filter(row %% 5 == 0)

#Plot
gp_data <- ggplot(DMSO_PTC_data_and_medians, aes(x = row, y = norm_mean, fill = mutation)) +
  geom_point(size = 2.0, color="black", shape=21, stroke = 0.2, alpha=0.8) +
  scale_fill_manual(values = c("TAA" = "#66c2a5",
                              "TAG" = "#fc8d62",
                              "TGA" = "#8da0cb")) +
  geom_errorbar(aes(ymin = norm_mean - std_val, ymax = norm_mean + std_val),
                width = 0.2, linewidth = 0.1) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) +
  scale_y_continuous(limits = c(-0.2, 1),
                     breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(
    name = "Codon",
    breaks = tick_labels_df$row,
    labels = tick_labels_df$codon
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gp_data


 

