#author: sool_lee@unc.edu 

read_bc_count_file_in_list <- function(list_input, filename_input) {
  for (i in 1:length(filename_input)) {
    list_input[[i]] <- fread(filename_input[i])
  }

  return(list_input)
}

merging_batch <- function(set1, set2) {
  merged_data <- base::merge(set1, set2, by = "barcode", all = T)
  merged_data[is.na(merged_data)] <- 0
  merged_data$num_reads <- rowSums(merged_data[, -1])
  merged_data %<>% dplyr::select(barcode, num_reads)
  return(merged_data)
}

merge_count_files_in_list <- function(pooled_data_input, batch_list_input, pool_num_input, rep_num_input) {
  for (i in 0:(rep_num_input - 1)) {
    merged_data <- merging_batch(
      batch_list_input[[i * pool_num_input + 1]],
      batch_list_input[[i * pool_num_input + 2]]
    )

    if (pool_num > 2) {
      for (j in 3:pool_num) {
        if (!is.null(batch_list_input[i * pool_num + j][[1]])) {
          merged_data %<>% merging_batch(batch_list_input[[i * pool_num + j]])
        } else {}
      }
    }

    pooled_data_input[[i + 1]] <- merged_data
  }

  return(pooled_data_input)
}

merge_dna_and_rna_bc_counts <- function(combined_counts_input, pooled_data_list_input, rep_num_input) {
  for (i in 1:rep_num_input) {
    combined_counts_input <- left_join(combined_counts_input, pooled_data_list_input[[i]], by = "barcode")
  }

  return(combined_counts_input)
}

filter_by_individual_barcode <- function(data_input, rep_num_input) {
  for (i in 1:rep_num_input) {
    rna_num <- i + 2
    data_subset <- data_input[, c(1:2, rna_num, ncol(data_input))] %>%
      set_colnames(c("barcode", "dna", "rna", "variant"))
    data_subset %<>% filter(rna > 1)
    data_subset_by_batch <- data_subset %>%
      group_by(variant) %>%
      summarise(
        rna_sum = sum(rna), dna_sum = sum(dna),
        rna_to_dna = sum(rna) / sum(dna),
        log_rna_to_dna = log(sum(rna) / sum(dna))
      )

    if (i > 1) {
      filtered_data %<>% base::merge(data_subset_by_batch, by = "variant", all = TRUE)
    } else {
      filtered_data <- data_subset_by_batch
    }
  }

  filtered_data %<>% na.omit()

  return(filtered_data)
}
