
library(igraph)
library(dplyr)

# Rank all genes by Personalized PageRank proximity to RANKL core
# Steady‑state probability that a random walker, which repeatedly diffuses through the PPI but “restarts” at the core seeds, is found at that node.
# g: igraph PPI (undirected; weighted or unweighted)
# array_genes: optional character vector of genes to report (e.g., array content)
# damping: 1 - restart probability; default 0.85 is standard
rank_genes_by_rankl_ppr <- function(g, array_genes = NULL, damping = 0.85,
                                    seeds = rankl_core) {
  # Ensure undirected for PPI
  if (is_directed(g)) g <- as.undirected(g, mode = "collapse", edge.attr.comb = "mean")
  
  # Keep only nodes present in PPI
  seeds_in <- intersect(seeds, V(g)$name)
  if (length(seeds_in) == 0) stop("No RANKL core seeds found in the PPI node set.")
  if (length(seeds_in) < length(seeds)) {
    warning("Some seeds not found in PPI and were dropped: ",
            paste(setdiff(seeds, seeds_in), collapse = ", "))
  }
  
  # Personalized restart vector: uniform over seeds
  pers <- as.numeric(V(g)$name %in% seeds_in)
  pers <- pers / sum(pers)  # sums to 1 over seeds
  
  # Use weights if present, else igraph treats them as 1
  pr <- page_rank(g,
                  directed   = FALSE,
                  damping    = damping,
                  weights    = if ("weight" %in% edge_attr_names(g)) E(g)$weight else NULL,
                  personalized = pers)$vector
  
  # Graph-wide baseline n and 1/n (use the vector length from PR)
  n_graph  <- length(pr)
  uniform  <- 1 / n_graph
  
  df <- data.frame(
    gene = V(g)$name,
    ppr  = pr, # Personalized PageRank score with respect to core seeds (probability sums to 1)
    fold_over_uniform = pr / uniform,    # enrichment over 1/n
    log2_fold = log2(pmax(pr / uniform, 1e-6)),  # log2 enrichment
    stringsAsFactors = FALSE
  )
  
  # Restrict to array genes
  if (!is.null(array_genes)) {
    df <- df[df$gene %in% array_genes, , drop = FALSE]
  }
  
  # Percentile rank within the returned set
  df$percentile <- rank(-df$ppr, ties.method = "average") / nrow(df) * 100
  
  df[order(-df$ppr), ]
}

empirical_pvals_degree_binned <- function(g,
                                          seeds,
                                          array_genes = NULL,
                                          B = 5000,
                                          damping = 0.85,
                                          nbins = 10,
                                          exclude_true_seeds = TRUE,
                                          rng_seed = 1,
                                          verbose = TRUE) {
  seeds <- as.character(seeds)
  
  # Observed scores
  df_obs <- rank_genes_by_rankl_ppr(
    g, array_genes = array_genes, damping = damping, seeds = seeds
  )
  pr_obs <- setNames(df_obs$ppr, df_obs$gene)
  genes_universe <- df_obs$gene
  
  # Degree on graph nodes
  deg <- degree(g)
  names(deg) <- V(g)$name
  
  # Seeds present in the graph
  seeds_in <- intersect(seeds, names(deg))
  k <- length(seeds_in)
  if (k == 0) stop("No seeds present in the graph.")
  
  # Sampling pool
  pool <- V(g)$name
  if (!is.null(array_genes)) pool <- intersect(pool, array_genes)
  if (exclude_true_seeds)    pool <- setdiff(pool, seeds_in)
  if (length(pool) < k) stop("Sampling pool smaller than seed set size; relax options.")
  
  # Degree bins on the pool
  unique_deg <- sort(unique(deg[pool]))
  nbins_eff <- min(nbins, max(2, length(unique_deg)))  # avoid too many bins
  qs <- stats::quantile(deg[pool], probs = seq(0, 1, length.out = nbins_eff + 1), na.rm = TRUE)
  breaks <- unique(as.numeric(qs))
  if (length(breaks) < 3) {
    # fallback to equal-width bins
    rng <- range(deg[pool], na.rm = TRUE)
    breaks <- unique(seq(rng[1], rng[2], length.out = 3))
  }
  
  bin_pool  <- cut(deg[pool],  breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_seeds <- cut(deg[seeds_in], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Count seeds per bin
  seed_bin_counts <- table(bin_seeds)
  
  if (verbose) {
    ks <- suppressWarnings(stats::ks.test(deg[seeds_in], deg[pool]))
    message(sprintf("Degree-binned permutations: k=%d, nbins=%d; mean deg seeds=%.2f, pool=%.2f; KS p=%.3g",
                    k, length(breaks) - 1, mean(deg[seeds_in]), mean(deg[pool]), ks$p.value))
  }
  
  # Pool by bin: split pool node names by their bin
  pool_by_bin <- split(pool, bin_pool)
  
  # Allocate null matrix
  null_mat <- matrix(0, nrow = length(genes_universe), ncol = B,
                     dimnames = list(genes_universe, NULL))
  
  set.seed(rng_seed)
  for (b in seq_len(B)) {
    seeds_b <- character(0)
    
    for (bb in names(seed_bin_counts)) {
      need <- as.integer(seed_bin_counts[[bb]])
      bb_int <- as.integer(bb)
      
      # candidates in the same bin
      candidates <- setdiff(pool_by_bin[[bb]], seeds_b)
      
      # If sparse, borrow from neighboring bins
      if (length(candidates) < need) {
        expand <- 1
        candidates_all <- candidates
        while (length(candidates_all) < need && expand <= (length(breaks) - 1)) {
          left  <- as.character(bb_int - expand)
          right <- as.character(bb_int + expand)
          add <- character(0)
          if (left  %in% names(pool_by_bin)) add <- c(add, setdiff(pool_by_bin[[left]],  seeds_b))
          if (right %in% names(pool_by_bin)) add <- c(add, setdiff(pool_by_bin[[right]], seeds_b))
          candidates_all <- unique(c(candidates_all, add))
          expand <- expand + 1
        }
        candidates <- candidates_all
      }
      
      if (length(candidates) < need) {
        stop(sprintf("Not enough candidates to match degree bin %s (need %d, have %d). Try reducing nbins.", bb, need, length(candidates)))
      }
      
      seeds_b <- c(seeds_b, sample(candidates, size = need, replace = FALSE))
    }
    
    # Run PPR with permuted seeds
    df_null <- rank_genes_by_rankl_ppr(
      g, array_genes = array_genes, damping = damping, seeds = seeds_b
    )
    null_mat[df_null$gene, b] <- df_null$ppr
    
    if (verbose && (b %% 500 == 0)) message(sprintf("Permutation %d/%d", b, B))
  }
  
  # Empirical p-values (upper tail)
  obs <- pr_obs[genes_universe]
  ge_counts <- rowSums(null_mat >= obs)
  p_emp <- (1 + ge_counts) / (B + 1)
  
  df_obs %>%
    mutate(p_emp = p_emp[gene],
           padj  = p.adjust(p_emp, method = "BH")) %>%
    arrange(padj, p_emp, desc(ppr))
}

# Some helpers for plotting
plot_ppr_base_epi <- function(ranked, dat_base, mean_met_epi, promoter_var, seed_tit, top_genes = 10){
  
  pdat <- ranked %>%
    dplyr::left_join(dat_base, by = "gene_primary") %>%
    dplyr::left_join(mean_met_epi, by = "gene_primary") %>%
    na.omit()
  
  # Rank Genes by PPR
  top_genes <- pdat %>%
    filter(promoter == promoter_var) %>%
    arrange(p_emp) %>%
    slice_head(n = top_genes)
  
  
  p <- pdat %>%
    dplyr::filter(promoter == promoter_var) %>%
    ggplot(aes(
      x = mean_deltaM_gene,
      y = -log10(p_emp),
      fill = diff
    )) +
    scale_fill_gradient2(
      name = "Mean Base Δ Epi",
      low = "#4575b4",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      breaks = pretty_breaks(n = 3)
    ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.5,
      size = 5
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene_primary),
      size = 3,
      color = "black",
      max.overlaps = Inf,
      segment.color = "grey50",
      segment.size = 0.3,
      segment.alpha = 0.8,
      nudge_y = 0.5,
      nudge_x = 0.05,
      box.padding = 0.2,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    theme_classic() +
    xlab("Mean Baseline Δ Exposure") +
    ylab(bquote(atop(
      -log[10] ~ "empirical p-value", "(PPR proximity to " ~ .(seed_tit) ~ ")"
    ))) +
    theme(legend.position = "bottom")
  
  return(p)
  
}

plot_ppr_base_pval_logfold <- function(ranked, dat_base, promoter_var, seed_tit, top_genes = 10){
  
  pdat <- ranked %>%
    dplyr::left_join(dat_base, by = "gene_primary") %>%
    na.omit()
  
  # Select top genes
  top_genes_df <- pdat %>%
    dplyr::filter(promoter == promoter_var) %>%
    dplyr::arrange(p_emp) %>%
    dplyr::slice_head(n = top_genes)
  
  p <- pdat %>%
    dplyr::filter(promoter == promoter_var) %>%
    ggplot(aes(
      x = mean_deltaM_gene,
      y = -log10(p_emp),
      fill = log2_fold
    )) +
    scale_fill_viridis_c(
      name = "Log2 Fold Enrichment",
      option = "C"
    ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.5,
      size = 5
    ) +
    geom_text_repel(
      data = top_genes_df,
      aes(label = gene_primary),
      size = 3,
      color = "black",
      max.overlaps = Inf,
      segment.color = "grey50",
      segment.size = 0.3,
      segment.alpha = 0.8,
      nudge_y = 0.5,
      nudge_x = 0.05,
      box.padding = 0.2,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    theme_classic() +
    xlab("Gene-level Mean ΔExposure (IPW)") +
    ylab(bquote(atop(
      -log[10] ~ "empirical P value",
      "(PPR proximity to " ~ .(seed_tit) ~ ")"
    ))) +
    theme(legend.position = "bottom")
  
  return(p)
}
