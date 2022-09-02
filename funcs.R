plot_distribution <- function(cur_stat) {
  
  ggplot(df, aes(x = eval(parse(text = cur_stat)))) +
    geom_histogram(bins = 40) +
    labs(x = cur_stat,
         y = "Count") +
    theme_bw()
  
}

freq_mapper <- function(cur_pos) {
  
  freq_array <- cur_pos %>%
    select(contains(paste0("Freq", cur_pos$con_base))) %>%
    as.vector() %>%
    as.numeric()
  
  return(freq_array)
  
}

max_mapper <- function(cur_pos, pattern) {
  
  max_val <- cur_pos %>%
    select(contains(paste0(pattern, "_", cur_pos$AbundMaxIdx, "."))) %>%
    as.vector() %>%
    as.numeric()
  
  return(max_val)
  
}

sample_loader <- function(sample_name) {
  
  sample_df <- read_csv2(sample_name) %>%
    mutate(sample = gsub(".csv", "", gsub("data/", "", sample_name)),
           genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)))
  
  return(sample_df)
  
}

sample_loader2 <- function(sample_name) {
  
  sample_df <- read_tsv(sample_name) %>%
    mutate(sample = gsub("_1.paired.fastq.gz_profile.IS_scaffold_info.tsv", "", gsub("data/scaffold_info/reference.fa.", "", sample_name)))
  
  return(sample_df)
  
}

gene_loader <- function(sample_name, gen_list = "", method = read_csv) {
  
  sample_df <- method(sample_name) %>%
    mutate(sample = gsub("data/gene_info/", "", gsub("_gene_info.tsv", "", sample_name))) %>%
    left_join(mapping)
  
  return(sample_df)
  
}

#' Hclust cannot handle matrices in which for some pairs of rows and columns,
#' only 1 or fewer shared values are non-NA. This function recurrently
#' identifies the most aggravating column/row, excludes that column/row and checks
#' whether more columns/rows need to be excluded
#'
#' @param mat Matrix to investigate
#' @param min_shared_fields Minimum number of positions that are not NA in both
#' vectors in order not to flag the vector pair as problematic
#'

identify_problematic_combs <- function(mat, min_shared_fields = 1) {
  exclude_rows <- NULL
  exclude_cols <- NULL
  stopifnot(is.matrix(mat))
  
  ## Loop over candidate removals
  for (k in 1:nrow(mat)) {
    candidate_rows <- setdiff(1:nrow(mat), exclude_rows)
    problem_row_combs <- NULL
    for (i in candidate_rows) {
      i_idx <- which(candidate_rows == i)
      for (j in candidate_rows[i_idx:length(candidate_rows)]) {
        if (sum(!is.na(mat[i, ]) & !is.na(mat[j, ])) <= min_shared_fields) {
          problem_row_combs <- rbind(problem_row_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_row_combs)) break
    exclude_rows <- c(exclude_rows,
                      as.integer(names(which.max(table(problem_row_combs)))))
  }
  
  for (k in 1:ncol(mat)) {
    candidate_cols <- setdiff(1:ncol(mat), exclude_cols)
    problem_col_combs <- NULL
    for (i in candidate_cols) {
      i_idx <- which(candidate_cols == i)
      for (j in candidate_cols[i_idx:length(candidate_cols)]) {
        if (sum(!is.na(mat[, i]) & !is.na(mat[, j])) <= min_shared_fields) {
          problem_col_combs <- rbind(problem_col_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_col_combs)) break
    exclude_cols <- c(exclude_cols,
                      as.integer(names(which.max(table(problem_col_combs)))))
  }
  
  return(list('row' = exclude_rows, 'column' = exclude_cols))
}


remove_problematic_combs <- function() {
  problematic_combs <- identify_problematic_combs(
    mat = mat, min_shared_fields = min_shared_fields)
  if (!is.null(problematic_combs$row)) {
    mat <- mat[-problematic_combs$row, ]
  }
  if (!is.null(problematic_combs$column)) {
    mat <- mat[, -problematic_combs$column]
  }
  return(mat)
}
formals(remove_problematic_combs) <- formals(identify_problematic_combs)
