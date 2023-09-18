library(data.table)
library(tidyverse)
library(mpra)
library(magrittr)
library(config)

#######################
# LOADING ENV & FUNCS #
#######################
# loading config and getting wd
config <- config::get(file = "where_the_config_file_is_stored/config_scz.yml")
current_wd <- str_c(config$scz_file_foldername)

# loading functions
source(str_c(config$script_dir, config$read_file_function_scriptname))
source(str_c(config$script_dir, config$create_count_matrix_functions_scriptname))

# loading files
variant_stat <- fread(config$scz_variant_statistics_filename)
dna_count_file <- fread(config$scz_dna_count_filepath)

# LOADING BARCODE DATA
####################
setwd(current_wd)
bcdat <- fread(config$scz_barcode_rda_filename)

# grepping file name and generate list with files
barcode_filenames <- list.files(pattern = config$bc_file_pattern)
barcode_file_list <-
  vector(mode = "list", length = length(barcode_filenames)) %>%
  read_bc_count_file_in_list(barcode_filenames)

# MERGE AND FILTER
###################
combined_data <-
  dna_count_file %>%
  merge_dna_and_rna_bc_counts(., barcode_file_list, 10) %>%
  set_colnames(c("barcode", "dna", str_c("rna", 1:10)))

na_count <- apply(combined_data[, 3:ncol(combined_data)], 1, function(x) sum(is.na(x)))

combined_data_na_filtered <- 
  combined_data[na_count < 9, ] %>%
  mutate_if(is.numeric, funs(replace(., is.na(.), 0)))

rm(barcode_file_list)

# ADD VARIANT DATA
###################
combined_data_na_filtered$variant <-
  bcdat[
    match(
      combined_data_na_filtered$barcode,
      bcdat$barcode
    ),
    "variant"
  ]

combined_data_na_filtered <-
  combined_data_na_filtered %>%
  as.data.frame() %>%
  na.omit() %>%
  filter(grepl("*.scz_", variant))

combined_data_final <-
  filter_by_individual_barcode(combined_data_na_filtered, 10) %>%
  set_colnames(c(
    "variant",
    str_c(
      rep(c("rna", "dna", "rna_to_dna", "log_rna_to_dna"), 10),
      rep(1:10, each = 4)
    )
  ))

combined_data_final$SNP <-
  unlist(lapply(
    strsplit(combined_data_final$variant, split = "_"),
    "[[", 2
  ))

combined_data_final$allele <-
  unlist(lapply(
    strsplit(combined_data_final$variant, split = "_"),
    "[[", 4
  ))

# SNPS FILTER (exclude data with only 1 SNP)
###################
snps <- table(combined_data_final$SNP)
snps <- names(snps[snps == 1])
combined_data_final <-
  combined_data_final[!(combined_data_final$SNP %in% snps), ] %>%
  as.data.frame()

# CREATE MATRICES AND VARIANT INFO
###################
dna_count_matrix <- create_count_matrix(combined_data_final, "dna", 10)
rna_count_matrix <- create_count_matrix(combined_data_final, "rna", 10)

variant_ids <- rownames(dna_count_matrix)

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

fwrite(mpra_result, config$scz_mpra_output_file, sep = "\t")
