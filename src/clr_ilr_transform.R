
library(compositions)
library(zCompositions)

clr_transform <- function(df, outcome_cols, meta_var){
  
  # Extract compositional data
  comp <- as.matrix(df[, outcome_cols])  # ensure numeric matrix
  comp <- apply(comp, 2, as.numeric)
  
  # Remove columns that are all 0
  all_zero_cols <- apply(comp, 2, function(x) all(is.na(x) | x == 0))
  comp <- comp[, !all_zero_cols, drop = FALSE]
  
  # Replace remaining zeros with small pseudo-count
  pseudo_count <- 1e-6
  comp[comp == 0] <- pseudo_count
  
  # Renormalize rows to sum to 1
  comp <- comp / rowSums(comp)
  
  # CLR transform
  clr_mat <- clr(comp)
  clr_df <- as.data.frame(clr_mat)
  
  # Combine with metadata
  out <- dplyr::bind_cols(
    df %>% dplyr::select(all_of(meta_var)),
    clr_df
  )
  
  return(out)
  
}

ilr_part_vs_rest_transform <- function(df,
                                       outcome_cols,
                                       meta_var,
                                       use_czm = FALSE,   # use multiplicative zero replacement if available
                                       pseudocount = 1e-6) {
  
  # Extract compositional data
  comp <- as.matrix(df[, outcome_cols, drop = FALSE])
  storage.mode(comp) <- "numeric"
  
  # Drop columns that are all 0/NA
  keep_cols <- apply(comp, 2, function(x) !all(is.na(x) | x == 0))
  comp <- comp[, keep_cols, drop = FALSE]
  
  # Drop rows that are all 0/NA
  keep_rows <- apply(comp, 1, function(x) !all(is.na(x) | x == 0))
  if (any(!keep_rows)) warning("Dropping rows with all-zero/NA parts.")
  comp <- comp[keep_rows, , drop = FALSE]
  df_kept <- df[keep_rows, , drop = FALSE]
  
  # Zero replacement + closure
  if (use_czm && requireNamespace("zCompositions", quietly = TRUE)) {
    comp <- zCompositions::cmultRepl(comp, method = "CZM", output = "prop")
  } else {
    comp[comp == 0 | is.na(comp)] <- pseudocount
    comp <- comp / rowSums(comp)
  }
  
  # Compute part-vs-rest balances for all parts
  D <- ncol(comp)
  L <- log(comp)
  
  # b_i = sqrt((D-1)/D) * (log x_i - mean log others)
  scale_fac <- sqrt((D - 1) / D)
  balances <- sapply(seq_len(D), function(i) {
    scale_fac * (L[, i] - rowMeans(L[, -i, drop = FALSE]))
  })
  colnames(balances) <- colnames(comp) 
  balances_df <- as.data.frame(balances)
  
  # Combine with metadata
  out <- dplyr::bind_cols(
    df_kept %>% dplyr::select(dplyr::all_of(meta_var)),
    balances_df
  )
  
  return(out)
}
