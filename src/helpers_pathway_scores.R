

compute_pathway_score_mean <- function(beta_gene_mat,   # genes × mice (rows = genes, cols = mice)
                                       pathway_genes,   # character vector of genes in the pathway
                                       baseline_samples # vector of sample names to define baseline (controls + cancer)
) {
  # Subset to pathway genes
  common_genes <- intersect(rownames(beta_gene_mat), pathway_genes)
  if(length(common_genes) == 0) stop("No overlapping genes between beta matrix and pathway_genes")
  
  beta_mat <- beta_gene_mat[common_genes, , drop = FALSE]
  
  #Compute average across genes for each sample
  pathway_score <- colMeans(beta_mat, na.rm = TRUE)
  
  # Scale using baseline samples
  baseline_score <- pathway_score[baseline_samples]
  mean_base <- mean(baseline_score, na.rm = TRUE)
  sd_base   <- sd(baseline_score, na.rm = TRUE)
  
  pathway_score_scaled <- (pathway_score - mean_base) / sd_base
  
  # Return named vector
  return(pathway_score_scaled)
}