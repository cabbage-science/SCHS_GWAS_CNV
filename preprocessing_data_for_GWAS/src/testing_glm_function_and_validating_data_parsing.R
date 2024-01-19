###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

#############################################
#### CHECKING IF DATA PARSING IS CORRECT ####
#############################################

# WORKS AS INTENDED!!!!!!!!!

column_names <- c("chr", "start", "end")
# Read BED file into a data frame
chromosomes <- paste0("chr",(seq(1:22)))
bed_data <- read.table("../data/binned_hg19_genome_200kb_10kb.bed", header = FALSE, col.names = column_names, stringsAsFactors = FALSE) %>%
  # Remove rows containing contigs/misc windows, and X/Y chromosomes
  filter(chr %in% chromosomes) %>%
  # Mini section to arrange chromosome numbers in ascending order
  mutate(no = as.integer(str_replace(chr,"chr",""))) %>%
  arrange(no) %>% select(-no) ; bed_data <- bed_data %>%
  # Creating a column to give each window a name
  mutate(window_name = paste0("CNVW",seq(1:nrow(bed_data))))


CNV_colnames <- c("cnv_id", "num_snp", "length", "state_cn", "indiv_ID", "start", "end")
all_cnv <- data.frame(fread("../data/CNV_batch_1_5_data_all_indivs.txt", col.names = CNV_colnames))
all_cnv_chr22 <- all_cnv %>%
  filter(str_detect(cnv_id, "chr22")) ### FIX THIS

matrix_chr22 <- data.frame(fread("../data/CNV_matrix_indivs_windows_chr22.txt"))

all_cnv_chr22[13,]

matrix_chr22[matrix_chr22$sample_ID == "225_SCHS_11499_01_GlobalScreening_iA1_1_GlobalScreening","CNVW285558"]
matrix_chr22[matrix_chr22$sample_ID == "225_SCHS_02990_01_GlobalScreening_iA1_1_GlobalScreening",4290]

bed_data[bed_data$window_name == "CNVW287264",]

##############################
#### TESTING GLM FUNCTION ####
##############################

?glm

telomere <- data.frame(fread("../data/Phenotype_GSA.txt")) %>%
  select(IID,tl_z,age_blood)

new_df <- telomere %>% left_join(matrix_chr22, by = c("IID"="sample_ID")) %>% 
  na.omit() %>% select(2,4301)

glm_test <- glm(new_df$tl_z ~ new_df[,2], family = gaussian()) 

list_a <- summary(glm_test) ; list_a

list_a$coefficients[paste0("new_df[, ",4299,"]"),"Pr(>|t|)"]

