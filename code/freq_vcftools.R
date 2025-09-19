library(tidyverse)
library(cowplot)
library(vcfR)
library(ggplot2)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Git/rye_project")
# id_map <- read.csv("data/ID_data/species_id_map.csv", header = TRUE)
# 
 # write.table(
 #   subset(id_map, species == "Secale strictum")[["Old_ID"]],
 #   file = "data/ID_data/Secale_strictum_IDs.txt",
 #   quote = FALSE,
 #   row.names = FALSE,
 #   col.names = FALSE
 # )


###########
#Histogramm
###########
library(vcfR)
library(ggplot2)

lines <- readLines("data/vcftools/freq_vcf.frq")
lines <- lines[-1]

parse_line <- function(line, maf_smallest = TRUE) {
  fields <- unlist(strsplit(line, "\\s+"))
  chrom <- fields[1]
  pos <- as.integer(fields[2])
  n_alleles <- as.integer(fields[3])
  n_chr <- as.integer(fields[4])
  allele_freqs <- fields[5:length(fields)]
  freqs <- as.numeric(sub(".*:", "", allele_freqs))
  maf <- NA
  if (maf_smallest) {
    maf <- min(freqs)
  } else {
    freqs_sorted <- sort(freqs, decreasing = TRUE)
    maf <- if (length(freqs_sorted) >= 2) freqs_sorted[2] else freqs_sorted[1]
  } 
  return(list(
    CHROM = chrom,
    POS = pos,
    N_ALLELES = n_alleles,
    N_CHR = n_chr,
    MAF = maf
  ))
}

parsed <- lapply(lines, parse_line)
df <- data.frame(
  CHROM = sapply(parsed, `[[`, "CHROM"),
  POS = sapply(parsed, `[[`, "POS"),
  N_ALLELES = sapply(parsed, `[[`, "N_ALLELES"),
  N_CHR = sapply(parsed, `[[`, "N_CHR"),
  MAF = sapply(parsed, `[[`, "MAF")
)
df <- df[!is.na(df$MAF) & df$MAF >= 0.01, ]


vcf <- read.vcfR("data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.sylvestre.recode.vcf")
vcf_variants <- as.data.frame(getFIX(vcf))  # enthält CHROM und POS als Strings
vcf_variants$POS <- as.integer(vcf_variants$POS)
sylvestre_with_maf <- merge(vcf_variants, df[, c("CHROM", "POS", "MAF")],
                            by = c("CHROM", "POS"), all = FALSE)

df_plot <- bind_rows(
  data.frame(MAF = df$MAF,            Gruppe = "Global"),
  data.frame(MAF = sylvestre_with_maf$MAF, Gruppe = "Sylvestre")
)

p <- ggplot() +
  geom_histogram(
    data = df_plot[df_plot$Gruppe == "Global", ],
    aes(x = MAF, color = "Global", fill = "Global"),
    binwidth = 0.01, boundary = 0
  ) +
  geom_histogram(
    data = df_plot[df_plot$Gruppe == "Sylvestre", ],
    aes(x = MAF, color = "Sylvestre"),
    binwidth = 0.01, boundary = 0, linewidth = 1, fill = NA
  ) +
  scale_color_manual(
    name   = NULL,
    values = c("Global" = "#377EB8", "Sylvestre" = "#FFFF33"),
    labels = c("Global", expression(italic("Secale sylvestre")))
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Global" = "#377EB8", "Sylvestre" = NA)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        fill   = c("#377EB8", NA),          
        color  = c("#377EB8", "#FFFF33"),   
        alpha  = c(1, 1),                 
        linewidth = c(1, 1)                 
      )
    ),
    fill = "none" 
  ) +
  labs(
    x = "Minor Allele Frequency",
    y = "Number of SNPs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.ticks.length= unit(2, "pt"),
    panel.grid       = element_blank(),
    legend.position  = c(0.75, 0.75)
  )

pdf("results/vcftools/maf_histogram.pdf", width = 6, height = 4)
print(p)
dev.off()






#############################
# Vergleich mit Allem anderen
#############################

vcf <- read.vcfR("data/vcftools/minDP20_maxDP100_minQ40_missing0.1_non-sylvestre.recode.vcf")
vcf_variants <- as.data.frame(getFIX(vcf))  # enthält CHROM und POS als Strings
vcf_variants$POS <- as.integer(vcf_variants$POS)
non_sylvestre_with_maf <- merge(vcf_variants, df[, c("CHROM", "POS", "MAF")],
                                by = c("CHROM", "POS"), all = FALSE)

sylvestre_with_maf$key <- paste(sylvestre_with_maf$CHROM, sylvestre_with_maf$POS, sep = "_")
non_sylvestre_with_maf$key <- paste(non_sylvestre_with_maf$CHROM, non_sylvestre_with_maf$POS, sep = "_")

shared_keys <- intersect(sylvestre_with_maf$key, non_sylvestre_with_maf$key)
only_sylvestre <- setdiff(sylvestre_with_maf$key, shared_keys)
write.csv(data.frame(SNP = only_sylvestre), "only_sylvestre.csv", row.names = FALSE)
only_sylvestre

sylvestre_with_maf$Gruppe <- ifelse(sylvestre_with_maf$key %in% shared_keys, "beide", "sylvestre")
non_sylvestre_with_maf$Gruppe <- ifelse(non_sylvestre_with_maf$key %in% shared_keys, "beide", "non-sylvestre")

sylvestre_unique <- sylvestre_with_maf[sylvestre_with_maf$Gruppe == "sylvestre", ]
non_sylvestre_unique <- non_sylvestre_with_maf[non_sylvestre_with_maf$Gruppe == "non-sylvestre", ]
shared <- sylvestre_with_maf[sylvestre_with_maf$Gruppe == "beide", ]

combined_maf <- rbind(shared, sylvestre_unique, non_sylvestre_unique)

# Chromosomen-Plot
library(dplyr)


combined_maf_sub <- combined_maf[combined_maf$MAF >= 0.01 & combined_maf$MAF < 0.03, ]
combined_maf_sub$MAF_Category <- cut(
  combined_maf_sub$MAF,
  breaks = c(0.01, 0.02, 0.03),
  labels = c("0.01–0.02", "0.02–0.03"),
  include.lowest = TRUE,
  right = FALSE
)

bin_size <- 25000000
combined_maf_sub$BIN <- floor(combined_maf_sub$POS / bin_size) * bin_size

agg_df <- combined_maf_sub %>%
  group_by(CHROM, BIN, Gruppe, MAF_Category) %>%
  summarise(SNP_Count = n()) %>%
  ungroup()


agg_df$Gruppe <- factor(agg_df$Gruppe, levels = c("sylvestre", "beide", "non-sylvestre"))
agg_df$MAF_Category <- factor(agg_df$MAF_Category, levels = c("0.01–0.02", "0.02–0.03"))
agg_df$Group_MAF <- interaction(agg_df$Gruppe, agg_df$MAF_Category, lex.order = TRUE)
agg_df$BIN <- factor(agg_df$BIN)

ggplot(agg_df, aes(x = BIN, y = SNP_Count, fill = Group_MAF)) +
  geom_col(position = "stack") +
  facet_wrap(~CHROM, scales = "free_x") +
  scale_fill_manual(
    name = "Gruppe & MAF",
    values = c(
      "beide.0.01–0.02" = "violet",
      "beide.0.02–0.03" = "purple",
      "sylvestre.0.01–0.02" = "red",
      "sylvestre.0.02–0.03" = "darkred",
      "non-sylvestre.0.01–0.02" = "blue",
      "non-sylvestre.0.02–0.03" = "darkblue"
    )
  ) +
  labs(
    title = "Verteilung von SNPs mit MAF < 0.03 entlang der Chromosomen (Fenstergröße: 25Mbp)",
    x = "Genomposition (binned)",
    y = "Anzahl SNPs"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)  # Schräge Achsentexte
  )