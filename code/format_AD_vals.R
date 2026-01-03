library(tidyverse)

setwd("/home/mie/Schreibtisch/Unistuff/Master/Semester2/Forschungsgruppenpraktikum_Steven/Test_Git/data/vcftools")

df <- read_tsv("sylvestre_snps.annotation_and_carriers.tsv", show_col_types = FALSE)


long <- df %>%
  select(VARIANT, CARRIERS_AD) %>%
  filter(!is.na(CARRIERS_AD), CARRIERS_AD != "") %>%
  separate_rows(CARRIERS_AD, sep = ";") %>%
  separate(CARRIERS_AD, into = c("sample", "AD"), sep = ":", fill = "right", extra = "merge")

ad_matrix <- long %>%
  pivot_wider(names_from = sample, values_from = AD) %>%
  arrange(VARIANT)

write_csv2(ad_matrix, "sylvestre_carriers_AD_matrix.csv")
