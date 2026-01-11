library(tidyverse)
library(PCAviz)
library(cowplot)
library(openxlsx)
library(stringr)
library(readr)
library(tibble)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/rye_project")

#################### Combine ID with species ####################

DIV_table <- read.xlsx("data/ID_data/12870_2025_6416_MOESM2_ESM.xlsx", sheet = 2)
Other_table <- read.xlsx("data/ID_data/20220603_Supplementary_File_1.xlsx", sheet = 1)
DIV_clean <- DIV_table %>%
  rename(ID = 1) %>%
  select(ID, accession, species, subspecies, domestication.level, countrycode, panel)
Other_clean <- Other_table %>%
  rename(ID = 1) %>%
  rename(
    subspecies = subtaxa,
    domestication.level = domestication.status
  ) %>%
  select(ID, accession, species, subspecies, domestication.level, countrycode, panel)

combined_table <- bind_rows(DIV_clean, Other_clean)
rm(DIV_clean,DIV_table,Other_clean,Other_table)

cleaning_ids <- function(old_id_vec) {
  x <- as.character(old_id_vec)
  x <- x %>%
    sub("\\..*$", "", .) %>%        # alles nach erstem Punkt weg
    sub("^[^_]*_", "", .)           # alles vor erstem "_" weg
  
  x <- vapply(x, function(id) {
    if (grepl("ME_L001", id)) {
      parts <- strsplit(id, "_")[[1]]
      if (length(parts) > 5) return(paste(head(parts, length(parts) - 5), collapse = "_"))
    } else if (grepl("DIV_22", id)) {
      parts <- strsplit(id, "_")[[1]]
      if (length(parts) > 4) return(paste(head(parts, length(parts) - 4), collapse = "_"))
    }
    id
  }, FUN.VALUE = character(1))
  x
}


files <- list.files(path = "data/pca_PLINK", pattern = "recode\\.eigenvec", full.names = TRUE)

old_ids <- lapply(files, function(f) {
  df <- read.table(f, stringsAsFactors = FALSE, colClasses = "character")
  tibble(Old_ID = df[[1]])
}) %>% bind_rows() %>% distinct()

id_map <- old_ids %>%
  mutate(ID = cleaning_ids(Old_ID)) %>%
  left_join(combined_table, by = "ID") %>%
  mutate(
    species = case_when(
      grepl("NPK_22|NF_22", ID)        ~ "Conduct",
      grepl("NPK_EG_23|NF_EG_23", ID)  ~ "Diversity Panel 23",
      grepl("LO7|Lo7|LO_7", ID)        ~ "Secale cereale",
      TRUE ~ species
    )
  ) %>%
  select(Old_ID, ID, species) %>%
  distinct()

id_map$old_species <- id_map$species

id_map <- id_map %>%
  mutate(
    species = case_when(
      grepl("NPK_22|NF_22", ID)        ~ "Secale cereale",
      grepl("NPK_EG_23|NF_EG_23", ID)  ~ "Secale cereale",
      grepl("LO7|Lo7|LO_7", ID)        ~ "Secale cereale",
      TRUE ~ species
    )
  )

# Kazakhstan samples form other_tab might be sylvestre as well and not strictum
ids_to_sylvestre <- c(
  "Sample_R1108_1_156_ME_L001_R1_001.unique.mapped.sorted.bam",
  "Sample_R1108_2_156_ME_L001_R1_001.unique.mapped.sorted.bam",
  "Sample_R1108_3_156_ME_L001_R1_001.unique.mapped.sorted.bam",
  "Sample_R1108_4_156_ME_L001_R1_001.unique.mapped.sorted.bam",
  "Sample_R1108_5_156_ME_L001_R1_001.unique.mapped.sorted.bam",
  "Sample_R1108_6_156_ME_L001_R1_001.unique.mapped.sorted.bam"
)
# context replacements
ids_to_cereale <- c(
  "3019347_R_1084_DIV_22_4_S102_L001_R1_001.unique.mapped.sorted.bam",
  "3019348_R_1185_DIV_22_4_S113_L001_R1_001.unique.mapped.sorted.bam",
  "3019501_R_1112_DIV_22_5_S345_L002_R1_001.unique.mapped.sorted.bam",
  "3019514_R_1112_DIV_22_5_S317_L002_R1_001.unique.mapped.sorted.bam",
  "3019531_R_1084_DIV_22_4_S329_L002_R1_001.unique.mapped.sorted.bam",
  "3019539_R_1185_DIV_22_4_S330_L002_R1_001.unique.mapped.sorted.bam"
)

id_map <- id_map %>%
  mutate(species = ifelse(Old_ID %in% ids_to_sylvestre, "Secale sylvestre", species)) %>%
  mutate(species = ifelse(Old_ID %in% ids_to_cereale, "Secale cereale", species))

write.csv(id_map, "data/ID_data/species_id_map.csv", row.names = FALSE)



#################### Plots ####################

id_map <- read.csv("data/ID_data/species_id_map.csv", header = TRUE)

#shapes <- c(16, 17, 15, 16, 17, 15, 16, 17, 15, 16, 17, 15)
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

nPCs <- 20
files <- list.files(path = "data/pca_PLINK", pattern = "recode\\.eigenvec", full.names = TRUE)
prefixes <- sub("\\.eigenvec", "", files)


maf_df <- tibble(prefix = prefixes) %>%
  mutate(
    maf_tag   = str_extract(prefix, "maf\\.[0-9.]+"),
    maf       = readr::parse_number(maf_tag),
    maf_label = ifelse(!is.na(maf), sprintf("MAF = %.2f", maf), NA_character_)
  ) %>%
  arrange(maf)

panel_letters <- letters[seq_len(nrow(maf_df) * 2)]

plots <- list()

for (i in seq_len(nrow(maf_df))) {
  prefix      <- maf_df$prefix[i]
  maf_label   <- maf_df$maf_label[i]
  eigvec_file <- paste0(prefix, ".eigenvec")
  eigval_file <- paste0(prefix, ".eigenval")
  
  pca <- read.table(eigvec_file, stringsAsFactors = FALSE)
  names(pca) <- c("Old_ID", "ID", paste0("PC", 1:nPCs))
  pca <- pca %>%
    left_join(id_map %>% select(Old_ID, species), by = "Old_ID") %>%
    filter(!is.na(species), !is.na(PC1), !is.na(PC2), !is.na(PC3))
  
  eigenval     <- sqrt(scan(eigval_file, quiet = TRUE))
  sum.eigenval <- sum(scan(eigval_file, quiet = TRUE))
  hgdp <- pcaviz(dat = pca, sdev = eigenval, var = sum.eigenval)
  
  base_opts <- list(
    scale_color_manual(values = coloures, drop = TRUE, na.translate = FALSE),
    theme(
      legend.position = "none",
      legend.text  = element_text(face = "italic"),
      legend.title = element_text(face = "bold"),
      plot.margin  = margin(4, 2, 2, 2),
      panel.spacing = unit(1, "mm")
    )
  )

  letter12 <- paste0("(", panel_letters[(i - 1) * 2 + 1], ")")
  letter13 <- paste0("(", panel_letters[(i - 1) * 2 + 2], ")")
  title12 <- bquote(bold(.(letter12)) ~ " " * .(maf_label))
  title13 <- bquote(bold(.(letter13)) ~ " " * .(maf_label))
  
  p12 <- plot(hgdp, coords = c("PC1","PC2"), color = "species",
              draw.points = TRUE, labels = FALSE, scale.pc.axes = 0.6) +
    labs(title = title12) + base_opts
  p13 <- plot(hgdp, coords = c("PC1","PC3"), color = "species",
              draw.points = TRUE, labels = FALSE, scale.pc.axes = 0.6) +
    labs(title = title13) + base_opts
  
  if (i == 3) {
    p13 <- p13 +
      theme(
        legend.position      = c(0.98, 0.75),
        legend.justification = c(1, 0.5),
        legend.background    = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size      = unit(4, "mm"),
        plot.margin          = margin(2, 8, 2, 2)
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  }
  
  plots <- append(plots, list(p12, p13))
}

outfile <- "results/PLINK/pca_PLINK.pdf"

pdf(outfile, width = 8.27, height = 11.69, onefile = TRUE)
g <- plot_grid(plotlist = plots, ncol = 2, align = "hv", 
                        rel_widths = c(1,1), rel_heights = rep(1, length(plots)/2))
print(g)
dev.off()



# Scree-Plot
multi_maf_screeplot <- function(csv_file) {
  df <- read.csv(csv_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  maf_names <- unlist(df[1, -1])             
  colnames(df) <- c("PC", maf_names)
  
  df <- df[-1, ]
  df$PC_num <- as.integer(sub("PC", "", df$PC))
  df[,-1] <- lapply(df[,-1], function(x) as.numeric(gsub(",", ".", x)))
  df$PC <- factor(df$PC, levels = df$PC)
  
  df_long <- df %>%
    pivot_longer(
      cols = -c(PC, PC_num),
      names_to = "MAF",
      values_to = "Eigenvalue"
    )
  df_long$MAF <- gsub(",", ".", df_long$MAF)
  
  p <- ggplot(df_long, aes(x = PC_num, y = Eigenvalue, color = MAF, group = MAF)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    labs(
      x = "Principal component",
      y = "Eigenvalue",
      title = ""
    ) +
    scale_color_manual(
      values = c(
        "0.01" = "#E41A1C",
        "0.02" = "#377EB8",
        "0.03" = "#4DAF4A",
        "0.04" = "#984EA3",
        "0.05" = "#FF7F00"
      ),
      name = "MAF"
    ) +
    scale_x_continuous(breaks = 1:20) +
    theme_classic(base_size = 14) +
    theme(
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.grid        = element_blank(),
      legend.title      = element_text(size = 12),
      legend.text       = element_text(size = 11),
      legend.position = c(0.85, 0.75),   # Position im Plot (0–1)
      legend.justification = c(1, 1),     # Bezugspunkt (rechts/oben)
      legend.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.3)
    )
  
  return(p)
}

p <- multi_maf_screeplot(csv_file = "results/PLINK/pca_PLINK_eigenval_PC1-20.csv")
pdf("results/PLINK/screeplot_PLINK.pdf", width = 8, height = 4, onefile = TRUE)
print(p)
dev.off()

