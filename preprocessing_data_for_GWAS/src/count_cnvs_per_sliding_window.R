###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

# THIS R SCRIPT CONTAINS THE WORKING FUNCTION USED FOR HPC ANALYSIS

# THE TEST FUNCTION IS IN THE OTHER R SCRIPT
#   Testing performed using the hg19 genome divided into bins of 10 Mb pairs, with step size of 5 Mb.
#   And a portion of the CNV data file for batches 1-3 (first 2500 rows)

#######################
#### LOADING FILES ####
#######################

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
  

# Read CNV data tsv file, with allocated column names 
# Purpose of colnames is for function to parse data effectively
CNV_colnames <- c("cnv_id", "num_snp", "length", "state_cn", "indiv_ID", "start", "end")
cnv_raw <- data.frame(fread("../data/CNV_batch_1_5_data_all_indivs.txt", 
                            header = FALSE, col.names = CNV_colnames))

##########################
#### WRITING FUNCTION ####
##########################

# Argument 1: df - CNV data in dataframe format - output from PennCNV 
# Argument 2: bed_df - df of binned genome created from bed file
# Argument 3: Minimum % of overlap of window size to be counted. Default: 10%
# Argument 4: Size of sliding window. Default: 200 kb
#   Output should be a data frame divided into genomic bins, with the copy number for each
#   individual. Default value would be 2, so deletion = 1/0 and duplication would be 3/4

count_cnvs_per_sliding_window <- function(cnv_df,bed_df,overlap=10,window_size=200000) {
  # Tidying up raw CNV data frame - for effective parsing
  cnv_df <- cnv_df %>% 
    separate(cnv_id, into = c("chr","bp"), sep = ":", convert = TRUE) %>%
    separate(bp, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    rename(cn = state_cn) %>%
    mutate(cn = as.integer(str_sub(cn,-1))) %>%
    select(-num_snp) %>%
    mutate(start = as.integer(str_replace(start, "startsnp=", "")),
           end = as.integer(str_replace(end, "endsnp=", "")),
           length = as.numeric(gsub("\\,","",str_replace(length, "length=", "")))) %>%
    filter(length > window_size*(overlap/100))
  #return(cnv_df)
  
  # Creating a vector of unique individual IDs
  sample_ids <- unique(cnv_df$indiv_ID)
  
  # Creating an "empty" combined df of sample IDs and copy numbers for each window (default = 2)
  output_df <- data.frame(unique(cnv_df$indiv_ID), 
                          matrix(2,length(sample_ids),nrow(bed_df),
                                 dimnames = list(NULL,bed_df$window_name))) %>%
    rename(sample_ID=unique.cnv_df.indiv_ID.)
  #return(output_df)
  
  # Count number of CNVs for sliding windows - repeating for each individual
  for(i in 1:length(sample_ids)) {
    # Creating a temporary df for individual of interest for more efficient computation
    temp_cnv_df <- cnv_df %>%
      filter(indiv_ID == sample_ids[i])
    # Running a loop for each of the 22 chromosomes - reduce computational time
    for(j in 1:length(unique(bed_df$chr))) {
      # Creating a temporary df for sliding windows in a single chromosome
      temp_bed_df <- bed_df %>% 
        filter(chr == paste0("chr",j))
      # Creating a temporary df for CNV data in a single chromosome
      temp_cnv_df_by_chr <- temp_cnv_df %>%
        filter(chr == paste0("chr",j))
      # Repeat the loop for every sliding window - effectively testing for overlap with 
      # every CNV in a particular chromosome
      for(k in 1:nrow(temp_bed_df)) {
        for(m in 1:nrow(temp_cnv_df_by_chr)) {
          # In case there are no CNVs in a chromosome for an individual. Prevent code from malfunctioning
          if(nrow(temp_cnv_df_by_chr) > 0) {
            # If the CNV is completely within a window
            if(temp_cnv_df_by_chr[m,"start"] >= temp_bed_df[k,"start"] & temp_cnv_df_by_chr[m,"end"] <= temp_bed_df[k,"end"]) {
              output_df[i,temp_bed_df[k,"window_name"]] <- temp_cnv_df_by_chr[m,"cn"]
              # If the bp position of end of the CNV is greater than the bp position of the end of the window
            } else if(temp_cnv_df_by_chr[m,"start"] >= temp_bed_df[k,"start"] & temp_cnv_df_by_chr[m,"end"] >= temp_bed_df[k,"end"]) {
              # Ensure that the overlap is greater than n% of the window, as defined by the argument
              if(temp_bed_df[k,"end"] - temp_cnv_df_by_chr[m,"start"] >= window_size*(overlap/100)) {
                output_df[i,temp_bed_df[k,"window_name"]] <- temp_cnv_df_by_chr[m,"cn"]
              }
              # If the bp position of start of the CNV is smaller than the bp position of the start of the window
            } else if(temp_cnv_df_by_chr[m,"start"] <= temp_bed_df[k,"start"] & temp_cnv_df_by_chr[m,"end"] <= temp_bed_df[k,"end"]) {
              # Ensure that the overlap is greater than n% of the window, as defined by the argument
              if(temp_cnv_df_by_chr[m,"end"] - temp_bed_df[k,"start"] >= window_size*(overlap/100)) {
                output_df[i,temp_bed_df[k,"window_name"]] <- temp_cnv_df_by_chr[m,"cn"]
              }
            }
          }
        }
      } # Then, loop repeats for next chromosome
    }
    # Output to track progress of parsing 
    print(paste0("Analysis for sample number ",i," out of ",length(sample_ids)," samples completed"))
  }
  return(output_df)
}

# Running the function on the CNV data and sliding window data
CNV_matrix <- count_cnvs_per_sliding_window(cnv_raw, bed_data)

# Writing a txt file for the matrix of CNV counts per individual per window
write.table(CNV_matrix, "CNV_matrix_indivs_windows.txt", quote = FALSE, sep = "\t", row.names = FALSE)

##################