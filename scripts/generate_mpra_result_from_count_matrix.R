library(data.table)
library(tidyverse)
library(mpra)
library(magrittr)
library(config)

#######################
# LOADING ENV & FUNCS #
#######################
# loading config and getting wd
config <- config::get(file = "where_the_config_file_is_stored/config_CD.yml")
current_wd <- str_c(config$cd_file_foldername)

# CREATE MATRICES AND VARIANT INFO
###################
dna_count_matrix <- fread(config$cd_dna_count_file)
rna_count_matrix <- fread(config$cd_rna_count_file)

variant_ids <- rownames(dna_count_matrix)

dna_count_matrix %<>% .[, -1]
rna_count_matrix %<>% .[, -1]

variant_stat %<>% .[, 1:2]
variant_stat$rsid <- unlist(lapply(strsplit(variant_stat$name, split = "_"), "[[", 2))
variant_seqs <- variant_stat[base::match(variant_ids, variant_stat$rsid), "variant"]

# CREATE DESIGN MATRIX
###################
design_matrix <-
  data.frame(
    intercept = rep(1, 10 * 2),
    alt = c(
      rep(TRUE, 10),
      rep(FALSE, 10)
    )
  )

samples <- rep(1:10, 2)

# RUN MPRA
###################
mpra_set <- MPRASet(
  DNA = dna_count_matrix,
  RNA = rna_count_matrix,
  eid = variant_ids,
  eseq = variant_seqs,
  barcode = NULL
)

mpra_lm <- mpralm(
  object = mpra_set,
  design = design_matrix,
  plot = T,
  aggregate = "none",
  normalize = T,
  block = samples,
  model_type = "corr_groups"
)

mpra_result <- topTable(mpra_lm, coef = 2, number = Inf)
mpra_result$variant <- rownames(mpra_result)

setwd(config$output_file_dir)

fwrite(mpra_result, config$cd_mpra_output_file, sep = "\t")
