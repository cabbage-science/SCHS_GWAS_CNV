###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

#######################
#### LOADING FILES ####
#######################

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

# Loading phenotype data and covariates with individual ID
telomere <- data.frame(fread("../data/Phenotype_GSA.txt")) %>%
  select(IID,tl_z,gender,age_blood,PC1,PC2,PC3)

##############################################
#### ASSOCIATION TESTING FOR DUPLICATIONS ####
##############################################

output_df <- data.frame()

for (j in 1:22) {
  # This loop repeats once per chromosome
  # Read in CNV data for sliding windows for a specific chromosome
  CNV_matrix_by_chr <- data.frame(fread(paste0("../data/CNV_matrix_indivs_windows_chr",j,".txt")))
  # Retrieve the number of sliding windows - which would be the number of independent tests
  ite <- ncol(CNV_matrix_by_chr) - 1
  # Create a temporary empty data frame of 2 columns - 1 for sliding window ID, 1 for p-value
  temp_output_df <- data.frame(matrix(NA, ncol(CNV_matrix_by_chr) - 1, 2))
  colnames(temp_output_df) <- c("window_ID","p_value")
  # Merging phenotype data with CNV data, only for individuals in common
  merged_df <- telomere %>% 
    inner_join(CNV_matrix_by_chr, by = c("IID"="sample_ID")) %>%
    # Some samples were ommited due to missing telomere length data
    na.omit()
  # Run independent GWAS tests
  for (k in 1:ite) {
    # Select telomere data, all covariates and CNV data for a single sliding window
    temp_df <- merged_df %>% select(2,3,4,5,6,7,7+k)
    # Provide current window name to an object for addition to the temporary empty df
    window_ID <- colnames(temp_df)[7]
    # Select rows only with deletions, since duplications are tested separately
    temp_data <- temp_df[temp_df[,7] %in% c(0,1,2),]
    # Condition to only run GWAS if number of unique copy numbers >= 2 for the sliding window
    if(length(unique(temp_data[,7])) > 1) {
      # Extract summary of results of the GLM - linear regression
      results <- summary(glm(temp_data[,1] ~ temp_data[,7]+ temp_data[,2] + temp_data[,3] + temp_data[,4] + temp_data[,5] + temp_data[,6],
                             data = temp_data, family = "gaussian"))
      # Extract p-value of the CNV from the summary of the GLM results
      p_value <- results$coefficients[2,"Pr(>|t|)"]
      # Building of temporary df for - meaning that there would be NA values
      temp_output_df[k,1] <- window_ID
      temp_output_df[k,2] <- p_value
    }
  }
  # Bind temporary df to output df for every chromosome
  output_df <- rbind(output_df, temp_output_df)
  print(paste0("Association analysis for chromosome ",j," completed!"))
}

# Need to omit NA values, which are sliding windows with copy number = 1
final_output_df <- output_df %>% na.omit() %>% left_join(bed_data, by = c("window_ID"="window_name"))

write.table(final_output_df,"../results/GLM_telomere_allcovar_DEL_remove_windows.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#### FINISH ####
