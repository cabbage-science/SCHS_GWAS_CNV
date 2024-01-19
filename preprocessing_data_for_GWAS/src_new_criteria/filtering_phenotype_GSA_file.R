library(tidyverse)
library(data.table)

# Filtering of phenotype file to exclude samples LRR_SD > 0.3, BAF_drift > 0.01, WF > |0.02| and CNV > 100.
# After editing the sample IDs such that only the sample number (?) (integer) remains, 
# and also removing ALL rows that had duplicates in their integer IDs

pheno <- data.frame(fread("../data/Phenotype_GSA.txt", header = TRUE))

pheno_cleaned <- pheno %>%
  mutate(IID = str_replace(IID, "225_", "")) %>%
  mutate(IID = str_replace(IID, "SCHS_", "")) %>%
  mutate(IID = str_replace(IID, "_01_GlobalScreening_iA1_1_GlobalScreening", "")) %>%
  mutate(IID = str_replace(IID, "_01_GlobalScreening_iA1_2_GlobalScreening", "")) %>%
  mutate(IID = str_replace(IID, "_03_GlobalScreening_iA1_1_GlobalScreening", "")) %>%
  mutate(IID = str_replace(IID, "_04_GlobalScreening_iA1_1_GlobalScreening", "")) %>%
  mutate(IID = str_replace(IID, "_02_GlobalScreening_iA1_1_GlobalScreening", "")) %>%
  mutate(IID = str_replace(IID, "startsnp", "")) %>%
  mutate(IID = str_replace(IID, "startsn", "")) %>%
  mutate(IID = str_replace(IID, "starts", "")) %>%
  mutate(IID = str_replace(IID, "start", "")) %>%
  mutate(IID = str_replace(IID, "star", "")) %>%
  #mutate(sample_ID = str_replace(sample_ID, "C", "")) %>%
  mutate(IID = as.integer(str_trim(IID)))

pheno_cleaned_test <- pheno_cleaned[!is.na(pheno_cleaned$IID),] %>% 
  arrange(desc(IID))

pheno_duplicated <- which(duplicated(pheno_cleaned_test$IID))
pheno_duplicated_ids <- pheno_cleaned_test[pheno_duplicated, "IID"]

pheno_final <- pheno_cleaned_test %>%
  filter(!IID %in% pheno_duplicated_ids) %>%
  filter(LRR_SD < 0.3, BAF_drift < 0.01, NumCNV < 100) %>%
  filter(WF > -0.02, WF < 0.02)

write.table(pheno_final, "../data/Phenotype_GSA_filtered.txt", row.names = FALSE, sep = "\t", quote = FALSE)
