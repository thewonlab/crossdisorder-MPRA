#author: @sool_lee@unc.edu

create_count_matrix <- function(data_input, data_type_input, rep_num_input) {
  if (data_type_input == "dna") {
    grep_pattern <- "^dna\\d"
  } else if (data_type_input == "rna") {
    grep_pattern <- "^rna\\d"
  }

  header_indeces <- colnames(data_input) %>% grep(pattern = grep_pattern)

  count_matrix <- cbind(
    data_input[seq(1, nrow(data_input) - 1, 2), header_indeces],
    data_input[seq(2, nrow(data_input), 2), header_indeces]
  ) %>%
    set_rownames(unique(data_input$SNP)) %>%
    set_colnames(c(
      paste(rep("alt", rep_num_input), paste0("s", 1:rep_num_input), sep = "_"),
      paste(rep("ref", rep_num_input), paste0("s", 1:rep_num_input), sep = "_")
    ))

  return(count_matrix)
}
