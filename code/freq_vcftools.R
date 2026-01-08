library(tidyverse)
library(cowplot)
library(vcfR)
library(ggplot2)
library(dplyr)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Test_Git")
id_map <- read.csv("data/ID_data/species_id_map.csv", header = TRUE)
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

lines <- readLines("data/vcftools/freq_vcf.90.frq")
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


vcf <- read.vcfR("data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.sylvestre.recode.vcf")
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
    values = c("Global" = "#377EB8", "Sylvestre" = "#E41A1C"),
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
        color  = c("#377EB8", "#E41A1C"),   
        alpha  = c(1, 1),                 
        linewidth = c(1, 1)                 
      )
    ),
    fill = "none" 
  ) +
  labs(
    x = "Minor Allele Frequency",
    y = "Number of Variants"
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

vcf <- read.vcfR("data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.non-sylvestre.recode.vcf")
non_fix <- as.data.frame(getFIX(vcf)) %>%
  transmute(CHROM, POS = as.integer(POS))
non_sylvestre_with_maf <- non_fix %>%
  distinct(CHROM, POS) %>%
  inner_join(df %>% select(CHROM, POS, MAF), by = c("CHROM","POS"))

shared <- inner_join(sylvestre_with_maf, non_sylvestre_with_maf,
                     by = c("CHROM","POS"), suffix = c(".s",".n")) %>%
          transmute(CHROM, POS, MAF = MAF.s, Gruppe = "both")

sylv_only <- anti_join(sylvestre_with_maf, non_sylvestre_with_maf,
                       by = c("CHROM","POS")) %>%
             mutate(Gruppe = "sylvestre") %>% select(CHROM, POS, MAF, Gruppe)

nons_only <- anti_join(non_sylvestre_with_maf, sylvestre_with_maf,
                       by = c("CHROM","POS")) %>%
             mutate(Gruppe = "non-sylvestre") %>% select(CHROM, POS, MAF, Gruppe)

# 7 Positionen zwar global mal als SNP aufgetaucht, aber in beiden Subgruppen kein beobachtbares ALT-Allel -> Fallen weg
combined_maf <- bind_rows(shared, sylv_only, nons_only)


bin_size <- 25e6
combined_maf_sub <- combined_maf %>%
  filter(MAF >= 0.01, MAF < 0.03) %>%
  mutate(
    MAF_Category = cut(MAF,
                       breaks = c(0.01, 0.02, 0.03),
                       labels = c("0.01-0.02","0.02-0.03"),
                       include.lowest = TRUE, right = FALSE),
    BIN = floor(POS / bin_size) * bin_size
  )

agg_df <- combined_maf_sub %>%
  count(CHROM, BIN, Gruppe, MAF_Category, name = "SNP_Count") %>%
  mutate(
    Gruppe       = factor(Gruppe, levels = c("sylvestre","both","non-sylvestre")),
    MAF_Category = factor(MAF_Category, levels = c("0.01-0.02","0.02-0.03")),
    Group_MAF    = interaction(Gruppe, MAF_Category, lex.order = TRUE),
    BIN          = factor(BIN, levels = sort(unique(BIN)))
  )

fills <- c(
  "sylvestre.0.01-0.02" = "#E41A1C",
  "sylvestre.0.02-0.03" = "darkred",
  "both.0.01-0.02" = "#F781BF",
  "both.0.02-0.03" = "#984EA3",
  "non-sylvestre.0.01-0.02" = "#377EB8",
  "non-sylvestre.0.02-0.03" = "darkblue"
)
brks <- names(fills)

p <- ggplot(agg_df, aes(x = BIN, y = SNP_Count, fill = Group_MAF)) +
  geom_col(position = "stack") +
  facet_wrap(~CHROM, scales = "free_x") +
  scale_fill_manual(
    name   = "Group and MAF",
    values = fills,
    breaks = brks,  # <- wichtig
    labels = c(
      expression(italic("S. sylvestre")~"only, 0.01-0.02"),
      expression(italic("S. sylvestre")~"only, 0.02-0.03"),
      expression("Both, 0.01-0.02"),
      expression("Both, 0.02-0.03"),
      expression("Non-"*italic("S. sylvestre")*" only, 0.01-0.02"),
      expression("Non-"*italic("S. sylvestre")*" only, 0.02-0.03")
    )
  ) +
  scale_x_discrete(
    breaks = function(x) x[seq(1, length(x), by = 3)],
    labels = function(x) sprintf("%.0f", as.numeric(as.character(x)) / 1e6)  
  ) +
  labs(
    title = "Window size: 25Mbp",
    x = "Genom position in Mbp (binned)",
    y = "Number of Variants"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x       = element_text(angle = 60, hjust = 1, color = "black"),
    axis.text.y       = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(size = 14,color = "black"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.text       = element_text(size = 14),     
    legend.title      = element_text(size = 14),    
    legend.key.size   = unit(4, "mm"),              
    legend.spacing.y  = unit(2, "mm"),              
    legend.spacing.x  = unit(2, "mm") 
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)  
  ) +
  theme(legend.position = c(0.5, 0.13))

pdf("results/vcftools/maf_sylvestre_vs_non-sylvestre.pdf", width = 14, height = 8)
print(p)
dev.off()


#############################
# Vergleich mit Strictum
#############################

vcf <- read.vcfR("data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.strictum.recode.vcf")
strictum_fix <- as.data.frame(getFIX(vcf)) %>%
  transmute(CHROM, POS = as.integer(POS))
strictum_with_maf <- strictum_fix %>%
  distinct(CHROM, POS) %>%
  inner_join(df %>% select(CHROM, POS, MAF), by = c("CHROM","POS"))

shared <- inner_join(sylvestre_with_maf, strictum_with_maf,
                     by = c("CHROM","POS"), suffix = c(".s",".n")) %>%
  transmute(CHROM, POS, MAF = MAF.s, Gruppe = "both")

sylv_only <- anti_join(sylvestre_with_maf, strictum_with_maf,
                       by = c("CHROM","POS")) %>%
  mutate(Gruppe = "sylvestre") %>% select(CHROM, POS, MAF, Gruppe)

stri_only <- anti_join(strictum_with_maf, sylvestre_with_maf,
                       by = c("CHROM","POS")) %>%
  mutate(Gruppe = "strictum") %>% select(CHROM, POS, MAF, Gruppe)

combined_maf <- bind_rows(shared, sylv_only, stri_only)

fills <- c(
  "sylvestre" = "#E41A1C",
  "both" = "#984EA3",
  "strictum" = "#377EB8"
  )

p <- ggplot(combined_maf, aes(x = MAF, fill = Gruppe)) +
  geom_histogram(binwidth = 0.01, boundary = 0, color = "black", position = "stack") +
  scale_fill_manual(values = fills,
    breaks = names(fills), 
    labels = c(
      expression(italic("S. sylvestre")~"only"),
      expression("Both"),
      expression(italic("S. strictum")~" only")
    )) +
  labs(
    title = "",
    x = "Minor Allele Frequency (MAF)",
    y = "Number of Variants",
    fill = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x       = element_text(angle = 60, hjust = 1, color = "black"),
    axis.text.y       = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.text       = element_text(size = 14),     
    legend.title      = element_text(size = 14),    
    legend.key.size   = unit(4, "mm"),              
    legend.spacing.y  = unit(2, "mm"),              
    legend.spacing.x  = unit(2, "mm") 
  )+
  theme(legend.position = c(0.8, 0.8))
pdf("results/vcftools/histogramm_sylvestre_vs_strictum.pdf", width = 6, height = 4)
print(p)
dev.off()



bin_size <- 25e6
combined_maf_sub <- combined_maf %>%
  filter(MAF >= 0.01, MAF < 0.03) %>%
  mutate(
    MAF_Category = cut(MAF,
                       breaks = c(0.01, 0.02, 0.03),
                       labels = c("0.01-0.02","0.02-0.03"),
                       include.lowest = TRUE, right = FALSE),
    BIN = floor(POS / bin_size) * bin_size
  )

agg_df <- combined_maf_sub %>%
  count(CHROM, BIN, Gruppe, MAF_Category, name = "SNP_Count") %>%
  mutate(
    Gruppe       = factor(Gruppe, levels = c("sylvestre","both","strictum")),
    MAF_Category = factor(MAF_Category, levels = c("0.01-0.02","0.02-0.03")),
    Group_MAF    = interaction(Gruppe, MAF_Category, lex.order = TRUE),
    BIN          = factor(BIN, levels = sort(unique(BIN)))
  )

fills <- c(
  "sylvestre.0.01-0.02" = "#E41A1C",
  "sylvestre.0.02-0.03" = "darkred",
  "both.0.01-0.02" = "#F781BF",
  "both.0.02-0.03" = "#984EA3",
  "strictum.0.01-0.02" = "#377EB8",
  "strictum.0.02-0.03" = "darkblue"
)

p <- ggplot(agg_df, aes(x = BIN, y = SNP_Count, fill = Group_MAF)) +
  geom_col(position = "stack") +
  facet_wrap(~CHROM, scales = "free_x") +
  scale_fill_manual(
    name   = "Group and MAF",
    values = fills,
    breaks = names(fills),  
    labels = c(
      expression(italic("S. sylvestre")~"only, 0.01-0.02"),
      expression(italic("S. sylvestre")~"only, 0.02-0.03"),
      expression("Both, 0.01-0.02"),
      expression("Both, 0.02-0.03"),
      expression(italic("S. strictum")~" only, 0.01-0.02"),
      expression(italic("S. strictum")~" only, 0.02-0.03")
    )
  ) +
  scale_x_discrete(
    breaks = function(x) x[seq(1, length(x), by = 3)],
    labels = function(x) sprintf("%.0f", as.numeric(as.character(x)) / 1e6)  
  ) +
  labs(
    title = "Window size: 25Mbp",
    x = "Genom position in Mbp (binned)",
    y = "Number of Variants"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x       = element_text(angle = 60, hjust = 1, color = "black"),
    axis.text.y       = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    strip.background  = element_rect(fill = "white", color = NA),
    strip.text        = element_text(size = 14,color = "black"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.text       = element_text(size = 14),     
    legend.title      = element_text(size = 14),    
    legend.key.size   = unit(4, "mm"),              
    legend.spacing.y  = unit(2, "mm"),              
    legend.spacing.x  = unit(2, "mm") 
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)  
  ) +
  theme(legend.position = c(0.8, 0.13))

pdf("results/vcftools/maf_sylvestre_vs_strictum.pdf", width = 14, height = 8)
print(p)
dev.off()

