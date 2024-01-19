library(qqman)
library(tidyverse)
library(data.table)

del_data <- data.frame(fread("../results/GLM_telomere_nocovar_DEL_remove_windows.txt", header=TRUE)) %>%
  mutate(chr = as.integer(str_replace(chr, "chr",""))) %>% na.omit()
dup_data <- data.frame(fread("../results/GLM_telomere_nocovar_DUP_remove_windows.txt", header=TRUE)) %>%
  mutate(chr = as.integer(str_replace(chr, "chr",""))) %>% na.omit()
del_data_allcovar <- data.frame(fread("../results/GLM_telomere_allcovar_DEL_remove_windows.txt", header=TRUE)) %>%
  mutate(chr = as.integer(str_replace(chr, "chr",""))) %>% na.omit()
dup_data_allcovar <- data.frame(fread("../results/GLM_telomere_allcovar_DUP_remove_windows.txt", header=TRUE)) %>%
  mutate(chr = as.integer(str_replace(chr, "chr",""))) %>% na.omit()
new_filtered_del_allcovar_data <- data.frame(fread("../results/GLM_telomere_allcovar_DEL_remove_windows_new_filtered.txt",
                                                   header = TRUE)) %>%
  mutate(chr = as.integer(str_replace(chr, "chr",""))) %>% na.omit()

Suggested_threshold <- 5

# DELETION NO COVAR
Mann_plot <- manhattan(
  del_data,
  chr = "chr",
  bp = "start",
  p = "p_value",
  snp = "window_ID",
  col = c("blue3","goldenrod"),
  annotateTop = T,
  suggestiveline = Suggested_threshold,
  main ="Manhattan Plot for Telomere Length - Deletions (n = 5996)",
  cex.main = 1.3, cex.axis = 1.2, cex.lab = 1.3, las = 1)

QQ_plot <- qq(del_data$p_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)


# DUPLICATION NO COVAR
Mann_plot <- manhattan(
  dup_data,
  chr = "chr",
  bp = "start",
  p = "p_value",
  snp = "window_ID",
  col = c("blue3","goldenrod"),
  annotateTop = T,
  suggestiveline = Suggested_threshold,
  main ="Manhattan Plot for Telomere Length - Duplications (n = 5996)",
  cex.main = 1.3, cex.axis = 1.2, cex.lab = 1.3, las = 1,
  ylim = c(0,6))

QQ_plot <- qq(dup_data$p_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)


# DELETION ALL COVAR 
Mann_plot <- manhattan(
  del_data_allcovar,
  chr = "chr",
  bp = "start",
  p = "p_value",
  snp = "window_ID",
  col = c("blue3","goldenrod"),
  annotateTop = T,
  suggestiveline = Suggested_threshold,
  main ="Manhattan Plot for Telomere Length - Deletions with all covariates (n = 5996)",
  cex.main = 1.3, cex.axis = 1.2, cex.lab = 1.3, las = 1)

QQ_plot <- qq(del_data_allcovar$p_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)



# DUPLICATION ALL COVAR 
Mann_plot <- manhattan(
  dup_data_allcovar,
  chr = "chr",
  bp = "start",
  p = "p_value",
  snp = "window_ID",
  col = c("blue3","goldenrod"),
  annotateTop = T,
  suggestiveline = Suggested_threshold,
  main ="Manhattan Plot for Telomere Length - Duplications with all covariates (n = 5996)",
  cex.main = 1.3, cex.axis = 1.2, cex.lab = 1.3, las = 1,
  ylim = c(0,8))

QQ_plot <- qq(dup_data_allcovar$p_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)



# WITH 19K SAMPLES DELETION ALL COVAR
Mann_plot <- manhattan(
  new_filtered_del_allcovar_data,
  chr = "chr",
  bp = "start",
  p = "p_value",
  snp = "window_ID",
  col = c("blue3","goldenrod"),
  annotateTop = T,
  suggestiveline = Suggested_threshold,
  main ="Manhattan Plot for Telomere Length - Duplications with all covariates (n = 5996)",
  cex.main = 1.3, cex.axis = 1.2, cex.lab = 1.3, las = 1)

QQ_plot <- qq(new_filtered_del_allcovar_data$p_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
