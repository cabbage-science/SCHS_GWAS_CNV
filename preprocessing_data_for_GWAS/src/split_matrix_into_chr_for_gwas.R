###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

#######################
#### LOADING FILES ####
#######################

CNV_matrix_1 <- data.frame(fread("../data/CNV_matrix_indivs_windows_1.txt", header = TRUE))
CNV_matrix_2 <- data.frame(fread("../data/CNV_matrix_indivs_windows_2.txt", header = TRUE))
CNV_matrix_3 <- data.frame(fread("../data/CNV_matrix_indivs_windows_3.txt", header = TRUE))
CNV_matrix_4 <- data.frame(fread("../data/CNV_matrix_indivs_windows_4.txt", header = TRUE))
CNV_matrix_5 <- data.frame(fread("../data/CNV_matrix_indivs_windows_5.txt", header = TRUE))
CNV_matrix_6 <- data.frame(fread("../data/CNV_matrix_indivs_windows_6.txt", header = TRUE))
CNV_matrix_7 <- data.frame(fread("../data/CNV_matrix_indivs_windows_7.txt", header = TRUE))
CNV_matrix_8 <- data.frame(fread("../data/CNV_matrix_indivs_windows_8.txt", header = TRUE))
CNV_matrix_9 <- data.frame(fread("../data/CNV_matrix_indivs_windows_9.txt", header = TRUE))
CNV_matrix_10 <- data.frame(fread("../data/CNV_matrix_indivs_windows_10.txt", header = TRUE))

# Defining colnames for bed file
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

##########################
#### TESTING FUNCTION ####
##########################

#indiv_ID <- paste0("indiv_",seq(1:10))
#CNVW <- paste0("CNVW",seq(1:22))
#temp_df <- data.frame(matrix(2,10,22))
#colnames(temp_df) <- CNVW 
#test_df <- cbind(indiv_ID,temp_df)

#fake_bed_data <- bed_data[1:22,]
#fake_bed_data[,1] <- chromosomes

#for(k in 1:22) {
#  temp_df <- fake_bed_data %>%
#    filter(chr %in% paste0("chr",k))
#  windows <- temp_df$window_name
#  CNV_matrix_chr <- test_df %>%
#    select(1,windows)
#  write.table(CNV_matrix_chr, paste0("CNV_matrix_indivs_windows_chr",k,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
#}

##############################################
#### SPLTTING CNV MATRIX INTO CHROMOSOMES ####
##############################################

CNV_matrix_combined <- rbind(CNV_matrix_1,CNV_matrix_2,CNV_matrix_3,CNV_matrix_4,CNV_matrix_5,
                             CNV_matrix_6,CNV_matrix_7,CNV_matrix_8,CNV_matrix_9,CNV_matrix_10)

rm(CNV_matrix_1,CNV_matrix_2,CNV_matrix_3,CNV_matrix_4,CNV_matrix_5,
   CNV_matrix_6,CNV_matrix_7,CNV_matrix_8,CNV_matrix_9,CNV_matrix_10)

for(k in 1:22) {
  temp_df <- bed_data %>%
    filter(chr %in% paste0("chr",k))
  windows <- temp_df$window_name
  CNV_matrix_chr <- CNV_matrix_combined %>%
    select(1,windows)
  write.table(CNV_matrix_chr, paste0("CNV_matrix_indivs_windows_chr",k,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}

#unique(bed_data$chr)

#### FINISH ####