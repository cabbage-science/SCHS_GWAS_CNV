###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(tidyverse)

#######################
#### LOADING FILES ####
#######################

# Initializing column names for pennCNV output table
CNV_colnames <- c("cnv_id", "num_snp", "length_bp", "state_cn", "indiv_ID", "start_snp", "end_snp", "confidence")

# Loading penn_CNV output file
pennCNV_output <- fread("../data/sample_cnv_data.txt", header = FALSE)
# Allocating column names
colnames(pennCNV_output) <- CNV_colnames

######################################################
#### FUNCTION FOR TIDYING UP PENNCNV OUTPUT TABLE ####
######################################################

separate_cols <- function(df) {
  df <- df %>% 
    mutate(temp_chrpos = cnv_id) %>%
    separate(temp_chrpos, into = c("chr","bp"), sep = ":", convert = TRUE) %>%
    separate(bp, into = c("start_bp", "end_bp"), sep = "-", convert = TRUE) %>%
    separate(state_cn, into = c("state", "cn"), sep = "," , convert = TRUE) %>%
    select(-c(confidence,state)) %>%
    mutate(num_snp = as.integer(str_replace(num_snp, "numsnp=", "")),
           length_bp = as.integer(str_replace(length_bp, "length=", "")),
           cn = as.integer(str_replace(cn, "cn=", "")),
           start_snp = str_replace(start_snp, "startsnp=", ""),
           end_snp = str_replace(end_snp, "endsnp=", ""))
  return(df)
}

# Calling function to tidy up pennCNV output into new object
new_pennCNV_output <- separate_cols(pennCNV_output)


###################################################################
#### FUNCTION TO CREATE MATRIX OF INDIVIDUALS AND COPY NUMBERS ####
###################################################################

create_indiv_cnv_df <- function(df) {
  df <- df %>%
    select(c(indiv_ID, cnv_id, cn)) %>%
    # Data is now in long format (tidy) - need to convert into wide table
    pivot_wider(names_from = cnv_id, 
                values_from = cn)
}

# Calling function on tidied pennCNV output to create a matrix of 
# individuals (rows) and copy numbers (columns)
CNV_data <- create_indiv_cnv_df(new_pennCNV_output)


#########################################################
#### FUNCTION TO CREATE RELATIONAL TABLE OF CNV INFO ####
#########################################################

create_cnv_info <- function(df) {
  df <- df %>%
    select(c(cnv_id, chr, start_bp, end_bp, num_snp, start_snp, end_snp))
  # Remove duplicated rows of CNVs
  df <- unique(df) %>%
    arrange(chr, start_bp)
}

# Calling function on tidied pennCNV output to "CNV_data"
CNV_info <- create_cnv_info(new_pennCNV_output)


################################################################
#### FUNCTION TO CALCULATE COUNT AND FREQUENCY (POPULATION) #### WORK IN PROGRESS !!!!!
################################################################

# Calculate count Variation 1 - creation of intermediate lists
calculate_CNV_counts <- function(df) {
  # Remove the individual ID column
  df <- df %>% select(-indiv_ID)
  # Some chatGPT magic that counts the number of individuals with specific copy numbers per CNV
  CNV_counts_long <- do.call(rbind, lapply(names(df), function(allele) {
    number_counts <- as.data.frame(table(df[[allele]]))
    number_counts$Allele <- allele
    return(number_counts)
  }))
  # Renaming the columns for clarity
  colnames(CNV_counts_long) <- c("Number", "Count", "CNV")
  # Converting into wider dataframe
  CNV_df <- CNV_counts_long %>%
    pivot_wider(names_from = Number, 
                values_from = Count)
  # Replacing all NA values with 0
  CNV_df[is.na(CNV_df)] <- 0
  return(CNV_df)
}

# ---------------------------------------------------

# Calculate count Variation 2 - no intermediate lists created, just a df
calculate_CNV_counts <- function(df) {
  # Remove the individual ID column
  df <- df %>% select(-indiv_ID)
  # Some chatGPT magic that counts the number of individuals with specific copy numbers per CNV
  # Initialize an empty dataframe
  CNV_counts_long <- data.frame(CNV = character(0), Number = character(0), Count = numeric(0), stringsAsFactors = FALSE)
  for (allele in names(df)) {
    number_counts <- as.data.frame(table(df[[allele]]))
    number_counts$Allele <- allele
    names(number_counts) <- c("Number", "Count", "CNV")
    # Bind this allele's data to the main dataframe
    CNV_counts_long <- rbind(CNV_counts_long, number_counts)
  }
  # Converting into wider dataframe
  CNV_df <- CNV_counts_long %>%
    pivot_wider(names_from = Number, 
                values_from = Count)
  # Replacing all NA values with 0
  CNV_df[is.na(CNV_df)] <- 0
  CNV_df <- CNV_df %>% mutate(Total = rowSums(CNV_df[,-1]))
  return(CNV_df)
}

CNV_counts <- calculate_CNV_counts(CNV_data)

# ===================================================
# Calculate frequency - based on the count table from above

calculate_CNV_frequency <- function(df) {
  df <- df %>% mutate(across(-c(CNV, Total), ~ ./Total)) %>%
    select(-Total)
  # Renaming column names for clarity 
  colnames(df) <- paste0("cn=",colnames(df))
  df <- rename(df, "CNV" = "cn=CNV")
  return(df)
}

CNV_frequency <- calculate_CNV_frequency(CNV_counts)
