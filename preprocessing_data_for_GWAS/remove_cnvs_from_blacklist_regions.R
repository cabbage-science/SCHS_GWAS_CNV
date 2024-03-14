###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

#######################
#### LOADING FILES ####
#######################

# Read CNV data CNV data file for first 5 batches
CNV_colnames <- c("cnv_id", "num_snp", "length", "state_cn", "indiv_ID", "start", "end")
cnv_raw <- data.frame(fread("../data/CNV_batch_1_5_data_all_indivs.txt", 
                            header = TRUE, col.names = CNV_colnames))
# Subsetting of raw CNV data for testing
#cnv_raw <- cnv_raw[1:5000,]

# Read in the "ENCODE Blacklist V2" data file - edited
# - Addition of telomere regions - first and last 5 million bp, according to this link
# - https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&chromInfoPage= 
# - Removal of lines pertaining to ChrX/Y, and contigs/supercontigs.
blacklist_raw <- read.table("../data/blacklist.txt", header = TRUE) %>%
  rename(end = chromEnd,
         start = chromStart)

# Read in the N-masked regions, in case the "blacklist" data file did not cover all N-masked regions as well
n_masked <- read.table("../data/nmasked.txt", header = TRUE) %>%
  rename(end = chromEnd,
         start = chromStart,
         chr = chrom)

################################################
#### REMOVING CNVS FROM BLACKLISTED REGIONS ####
################################################

# Combining blacklisted and n-masked data files
combined_blacklist <- n_masked %>%
  # Remove centromere and telomere rows as theyhave already been accounted for in blacklist_raw
  filter(!type %in% c("centromere","telomere")) %>%
  # Remove rows of chrX, Y and supercontigs/contigs
  filter(chr %in% paste0("chr",seq(1:22))) %>%
  # Select essential rows
  select(chr,start,end) %>%
  # Select essential rows
  rbind(blacklist_raw[,-4])

### Function 1: Remove rows that did not meet minimum requirements
# Argument 1: CNV raw data df
# Argument 2: Minimum size of CNV - default = 20000
# Argument 3: Minimum nSNPs covering the CNV - default = 10
# Argument 4: Minimum SNP density of CNV - default = 0.0001 = 1 in 10000
remove_CNVs_fail_basic_filters <- function(cnv_df, size=20000, nsnp=10, snpdensity=0.0001) {
  # Creating empty vector - to tabulate row numbers of CNVs for exclusion later on
  tabulate_row_numbers_1 <- NULL
  
  # Creating temp df for output
  cnv_df_for_output <- cnv_df
  
  # Tidying up raw CNV data frame for data analysis
  cnv_df <- cnv_df %>% 
    separate(cnv_id, into = c("chr","bp"), sep = ":", convert = TRUE) %>%
    separate(bp, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    select(-state_cn, -indiv_ID) %>%
    mutate(start = as.integer(str_replace(start, "startsnp=", "")),
           end = as.integer(str_replace(end, "endsnp=", "")),
           num_snp = as.integer(str_replace(num_snp, "numsnp=", "")),
           length = as.numeric(gsub("\\,","",str_replace(length, "length=", "")))) %>%
    mutate(snp_density = num_snp/length)
  print("Data cleaning completed!")
  
  # Extract row numbers that did not meet minimum SNP size OR minimum nSNPs OR minimum SNP density
  for(i in 1:nrow(cnv_df)) {
    if(cnv_df[i,"length"] < size | cnv_df[i,"num_snp"] < nsnp | cnv_df[i,"snp_density"] < snpdensity) {
      tabulate_row_numbers_1 <- c(tabulate_row_numbers_1, i)
    }
    
    # Track progress
    if(i %% 50000 == 0) {
      print(paste0("Screening for ", i, " out of ", nrow(cnv_df), " completed (Phase 1)!"))
    }
  }
  
  # Remove duplicated entries for row numbers - just in case
  tabulate_row_numbers <- unique(tabulate_row_numbers_1)
  
  # Remove "tagged" rows from the CNV df
  cnv_df_for_output <- cnv_df_for_output[-tabulate_row_numbers,]

  return(cnv_df_for_output)
  #return(tabulate_row_numbers)
  #return(cnv_df)
}

# Remove/add the comment hashtag where appropriate.
cnv_data_intermediate <- remove_CNVs_fail_basic_filters(cnv_raw, size=20000, nsnp=10)
cnv_data_intermediate <- remove_CNVs_fail_basic_filters(cnv_raw, size=20000, nsnp=5)
cnv_data_intermediate <- remove_CNVs_fail_basic_filters(cnv_raw, size=10000, nsnp=10)
cnv_data_intermediate <- remove_CNVs_fail_basic_filters(cnv_raw, size=10000, nsnp=5)


### Function 2: Remove rows of CNVs that had > 50% of itself overlap with blacklisted regions
# Argument 1: CNV raw data
# Argument 2: Combined blacklist data frame
# Argument 3: Overlap percentage - default = 50% / 0.5
remove_CNVs_inside_blacklist_regions <- function(cnv_df, blacklist_df, overlap=0.5) {
  # Creating empty vector - to tabulate row numbers of CNVs for exclusion later on
  tabulate_row_numbers_2 <- NULL
  
  # Creating temp df for output
  cnv_df_for_final_output <- cnv_df
  
  # Tidying up raw CNV data frame for data analysis
  cnv_df <- cnv_df %>% 
    separate(cnv_id, into = c("chr","bp"), sep = ":", convert = TRUE) %>%
    separate(bp, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    select(-state_cn, -indiv_ID, -num_snp) %>%
    mutate(start = as.integer(str_replace(start, "startsnp=", "")),
           end = as.integer(str_replace(end, "endsnp=", "")),
           length = as.numeric(gsub("\\,","",str_replace(length, "length=", ""))))
  print("Data cleaning completed!")
  
  for (i in 1:nrow(cnv_df)) {
    chr_no <- cnv_df[i,"chr"]
    temp_blacklist_df <- blacklist_df %>% 
      filter(chr == chr_no)
    for(j in 1:nrow(temp_blacklist_df)) {
      if(cnv_df[i,"start"] > temp_blacklist_df[j,"start"] & cnv_df[i,"end"] < temp_blacklist_df[j,"end"]) {
        tabulate_row_numbers_2 <- c(tabulate_row_numbers_2, i)
        next
      } else if(cnv_df[i,"start"] > temp_blacklist_df[j,"start"] & cnv_df[i,"end"] > temp_blacklist_df[j,"end"]) {
        if((temp_blacklist_df[j,"end"] - cnv_df[i,"start"])/cnv_df[i,"length"] > overlap) {
          tabulate_row_numbers_2 <- c(tabulate_row_numbers_2, i)
        }
        next
      } else if(cnv_df[i,"start"] < temp_blacklist_df[j,"start"] & cnv_df[i,"end"] < temp_blacklist_df[j,"end"]) {
        if((cnv_df[i,"end"] - temp_blacklist_df[j,"start"])/cnv_df[i,"length"] > overlap) {
          tabulate_row_numbers_2 <- c(tabulate_row_numbers_2, i)
        }
        next
      } else if(cnv_df[i,"start"] < temp_blacklist_df[j,"start"] & cnv_df[i,"end"] > temp_blacklist_df[j,"end"]) {
        if((temp_blacklist_df[j,"end"] - temp_blacklist_df[j,"start"])/cnv_df[i,"length"] > overlap) {
          tabulate_row_numbers_2 <- c(tabulate_row_numbers_2, i)
        }
        next
      }
    }
    
    # Track progress
    if(i %% 1000 == 0) {
      print(paste0("Screening for ", i, " out of ", nrow(cnv_df), " completed! (Phase 2)"))
    }
  }
  
  # Remove duplicated entries for row numbers - just in case
  tabulate_row_numbers <- unique(tabulate_row_numbers_2)
  
  # Remove "tagged" rows from the CNV df
  cnv_df_for_final_output <- cnv_df_for_final_output[-tabulate_row_numbers,]
  
  return(cnv_df_for_final_output)
  #return(tabulate_row_numbers)
  #return(cnv_df)
  #return(temp_blacklist_df)
}

cnv_raw_data_final <- remove_CNVs_inside_blacklist_regions(cnv_data_intermediate, blacklist_raw)

# Remove/add the comment hashtag where appropriate.
write.table(cnv_raw_data_final, "../data/CNV_batch_1_5_all_indivs_blacklist_filtered.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cnv_raw_data_final, "../data/CNV_batch_1_5_all_indivs_blacklist_filtered_20kb_5nsnp.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cnv_raw_data_final, "../data/CNV_batch_1_5_all_indivs_blacklist_filtered_10kb_10nsnp.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cnv_raw_data_final, "../data/CNV_batch_1_5_all_indivs_blacklist_filtered_10kb_5nsnp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# 19969 unique individual IDs, 78408 rows of CNV data!
# Is the filtering too stringent?