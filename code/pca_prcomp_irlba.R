library(vcfR)
library(irlba)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(cowplot)
library(grid)

setwd("~/rye_project")

vcf <- read.vcfR("data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted.vcf")
gt <- extract.gt(vcf,convertNA=TRUE,element='GT', as.numeric=FALSE, return.alleles=FALSE)
geno_num <- apply(gt, 2, function(x) {
  ifelse(x == "0/0", 0,
         ifelse(x %in% c("0/1", "1/0"), 1,
                ifelse(x == "1/1", 2, NA)))
})

geno_num[is.na(geno_num)] <- 0
geno_num <- t(geno_num)

pca_result <- prcomp_irlba(geno_num, n=20, center=TRUE, scale.=TRUE)
summary(pca_result)
rownames(pca_result$x) <- rownames(geno_num)

############ Screeplot ############

prop_var <- (pca_result$sdev^2) / pca_result$totalvar
df <- data.frame(
  PC = seq_along(prop_var),
  proportion = prop_var
)

p <- ggplot(df, aes(x = PC, y = proportion)) +
  geom_line(color = "#E41A1C") +
  geom_point(color = "#E41A1C") +
  labs(
    x = "Principal component",
    y = "Variance proportion",
    title = ""
  )+
  theme_classic(base_size = 13) +
  theme(
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.ticks.length= unit(2, "pt"),
    panel.grid       = element_blank()
  )

pdf("results/prcomp_irlba/pca_screeplot_prcomp_irlba.pdf", width = 6, height = 4)
print(p)
dev.off()

############ PCA-Plots ############

id_map <- read.csv("data/ID_data/species_id_map.csv", header = TRUE)

lab_pc1 <- paste0("PC1 (", round(prop_var[1],3), ")")
lab_pc2 <- paste0("PC2 (", round(prop_var[2],3), ")")
lab_pc3 <- paste0("PC3 (", round(prop_var[3],3), ")")

scores <- as_tibble(pca_result$x, rownames = "Old_ID") %>%
  select(Old_ID, PC1, PC2, PC3) %>%
  left_join(id_map %>% select(Old_ID, species), by = "Old_ID") %>%
  filter(!is.na(species), !is.na(PC1), !is.na(PC2), !is.na(PC3)) %>%
  mutate(species = factor(species))

coloures <- c(
  "#4DAF4A",  
  "darkblue",
  "#377EB8",  
  "#984EA3", 
  "#F781BF",  
  "#E41A1C",  
  "darkred",
  "#FF7F00"
)

base_opts <- list(
  scale_color_manual(values = coloures, drop = TRUE, na.translate = FALSE),
  theme_classic(base_size = 12),
  theme(
    legend.text  = element_text(face = "italic"),
    legend.title = element_text(face = "bold"),
    plot.title       = element_text(face = "bold"),
    plot.margin  = margin(4, 2, 2, 2),
    panel.spacing = unit(1, "mm"),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.ticks.length= unit(2, "pt"),
    panel.grid       = element_blank()
  )
)

p12 <- ggplot(scores, aes(PC1, PC2, color = species)) +
  geom_point(alpha = 0.9, size = 1.7) +
  labs(x = lab_pc1, y = lab_pc2, title = "(a)", color = "Species") +
  base_opts +
  theme(legend.position = "none")

p13 <- ggplot(scores, aes(PC1, PC3, color = species)) +
  geom_point(alpha = 0.9, size = 1.7) +
  labs(x = lab_pc1, y = lab_pc3, title = "(b)", color = "Species") +
  base_opts +
  theme(
    legend.position      = c(0.99, 0.87),
    legend.justification = c(1, 0.5),
    legend.background    = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
    legend.key.size      = unit(4, "mm")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

p <- plot_grid(p12, p13, ncol = 2, align = "hv")
pdf("results/prcomp_irlba/pca_prcomp_irlba.pdf", width = 9.27, height = 5)
print(p)
dev.off()



############ Loadings ############

loadings_df <- as.data.frame(pca_result$rotation[, 1:3])
loadings_df$SNP <- rownames(loadings_df)
loadings_long <- loadings_df %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "PC",
               values_to = "Loading")

p <- ggplot(loadings_long, aes(x = PC, y = Loading)) +
  geom_violin(fill = "#377EB8", alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_x_discrete(labels = c("PC1", "PC2", "PC3")) +
  labs(x = "Principal component", 
       y = "SNP loading",
       title = "") +
  theme_classic(base_size = 13) +
  theme(
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.ticks.length= unit(2, "pt"),
    panel.grid       = element_blank(),
    legend.position  = c(0.75, 0.75)
  )
pdf("results/prcomp_irlba/pc_loadings_violins_prcomp_irlba.pdf", width = 9.27, height = 5)
print(p)
dev.off()


target_snps <- read.csv("data/ID_data/only_sylvestre_SNPs_maf.01.csv")[[1]]
sylvestre_loadings_df <- as_tibble(pca_result$rotation[, 1:3]) %>%   # PC1–PC3
  mutate(
    SNP = sub("^(([^_]*_[^_]*)).*", "\\1", names(pca_result$center)),
    abs_PC1 = abs(PC1),
    abs_PC2 = abs(PC2),
    abs_PC3 = abs(PC3)
  ) %>%
  select(SNP, everything()) %>% 
  filter(SNP %in% target_snps)
write.csv(sylvestre_loadings_df, "results/prcomp_irlba/pca_loadings_only_sylvestre_prcomp_irlba.csv", row.names = FALSE)


n_snps_total <- nrow(loadings_df)
n_snps_pc2 <- sum(loadings_df$PC2 < -0.03, na.rm = TRUE)
cat("PC2 SNPs <", -0.03, ":", n_snps_pc2, "von", n_snps_total,
    sprintf("(%.2f%%)", (n_snps_pc2 / n_snps_total * 100)), "\n")
# Res: PC2 SNPs < -0.03 : 288 von 21548 (1.34%)

# stat. Nachweis:
#sylvestre SNPs > 0.03
n <- 9
k <- nrow(sylvestre_loadings_df)
p_one_sided <- phyper(k - 1, n_snps_pc2, n_snps_total - n_snps_pc2, n, lower.tail = FALSE)
p_one_sided
# Res: p=0
# Wahrscheinlichkeit, dass zufällig SNPs gewählt wurden, die größer als 0.03 haben beträgt ~0

#hist(loadings_df$PC2,breaks = 80, main="Distribution of loadings", xlab="Loading")

