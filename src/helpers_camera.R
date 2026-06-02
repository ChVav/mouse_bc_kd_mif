library(dplyr)
library(metafor)
library(limma)

camera_probe_bias_adjusted <- function(dmr_result, anno_mouse, all_cpg, my_sets, span = 0.75) {
  # Map
  map_bg <- anno_mouse %>% filter(Name %in% all_cpg, !is.na(gene_primary))
  gene_counts <- as.integer(table(map_bg$gene_primary))
  names(gene_counts) <- names(table(map_bg$gene_primary))
  G <- names(gene_counts)
  
  # Annotate
  dmr_annot <- dmr_result %>%
    inner_join(map_bg, by = c("outcome" = "Name")) %>%
    filter(!is.na(estimate), !is.na(SE))
  
  # Per-gene random-effects meta-analysis
  fit_gene_meta <- function(beta, se) {
    if (length(beta) == 1L) return(c(beta_gene = beta, SE_gene = se))
    res <- try(rma.uni(yi = beta, sei = se, method = "REML"), silent = TRUE)
    if (inherits(res, "try-error")) {
      w <- 1 / (se^2); c(sum(w * beta) / sum(w), sqrt(1 / sum(w))) |> setNames(c("beta_gene","SE_gene"))
    } else {
      c(beta_gene = as.numeric(res$b), SE_gene = res$se)
    }
  }
  
  gene_meta <- dplyr::group_by(dmr_annot, gene_primary) %>%
    dplyr::reframe({
      beta <- estimate; se <- SE
      out <- fit_gene_meta(beta, se)
      tibble(beta_gene = out["beta_gene"], SE_gene = out["SE_gene"], n_cpg = length(beta))
    }) %>%
    mutate(z_gene = beta_gene / SE_gene) %>%
    filter(gene_primary %in% G) %>%
    ungroup()
  
  # Residualize z on CpG count
  m <- gene_counts[gene_meta$gene_primary]
  dfz <- data.frame(z = gene_meta$z_gene, m = m)
  fit <- stats::loess(z ~ log1p(m), data = dfz, span = span, degree = 1)
  z_res <- dfz$z - predict(fit, newdata = dfz)
  gene_meta$z_adj <- z_res
  
  # Build indices and run cameraPR
  genes_order <- gene_meta$gene_primary
  index <- lapply(my_sets, function(g) which(genes_order %in% g))
  cam_res <- cameraPR(stat = setNames(gene_meta$z_adj, genes_order),
                      index = index, use.ranks = FALSE)
  
  list(
    gene_meta = gene_meta,  # beta_gene, SE_gene, z_gene, z_adj, n_cpg
    camera    = cam_res
  )
}

# Starting from dmr results
prep_camera_probe_bias_adjusted <- function(dmr_result, anno, gene_sets, region_var, region_val){
  # Filter promoter dmr
  dmr_result_flt <- dplyr::left_join(dmr_result, anno, by = c("outcome" = "Name")) %>%
    dplyr::select(outcome, estimate, SE, p.value, p.adj_cell, gene_primary, .data[[region_var]]) %>%
    dplyr::filter(.data[[region_var]] == region_val)
  
  # Filter hallmark CpGs & region dmr
  anno_mouse <- anno %>%
    dplyr::select(Name,gene_primary) %>%
    dplyr::full_join(gene_sets, by = "gene_primary", relationship = "many-to-many") %>%
    dplyr::filter(!is.na(gs_name)) %>% # all in gene_sets
    dplyr::filter(Name %in% dmr_result_flt$outcome) # filter by region
  
  # length(unique(anno_mouse$Name)) 
  # length(unique(anno_mouse$gs_name)) 
  # length(unique(anno_mouse$gene_primary)) 
  
  # Filter rest
  dmr_result_flt <- dmr_result_flt %>% dplyr::filter(outcome %in% anno_mouse$Name) %>% dplyr::select(-gene_primary, -.data[[region_var]])
  
  all_cpg <- dmr_result_flt %>% pull(outcome)
  
  gene_sets_flt <- gene_sets %>% dplyr::filter(gene_primary %in% anno_mouse$gene_primary)
  # length(unique(gene_sets_flt$gene_primary))
  # length(unique(gene_sets_flt$gs_name))
  
  anno_mouse <- anno_mouse %>% dplyr::select(-gs_collection, -gs_name) %>% distinct()
  
  sets_df <- gene_sets_flt %>%
    filter(!is.na(gs_collection), !is.na(gs_name), !is.na(gene_primary)) %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_primary)), .groups = "drop")
  my_sets <- setNames(sets_df$genes, sets_df$gs_name)

  return(list(dmr_result = dmr_result_flt, 
              anno_mouse = anno_mouse,
              all_cpg = all_cpg,
              my_sets = my_sets))
}

# Starting from gene summaries (for example generated with gene_meth_ipw)
run_camera <- function(gene_result, gene_sets) {
  
  gene_sets <- gene_sets %>% filter(!is.na(gs_collection), !is.na(gs_name), !is.na(gene_primary))
  keep <- intersect(unique(gene_result$gene_primary),unique(gene_sets$gene_primary))
  gene_meta <- gene_result %>% dplyr::filter(gene_primary %in% keep)
  gene_sets <- gene_sets %>% dplyr::filter(gene_primary %in% keep)
  
  sets_df <- gene_sets %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_primary)), .groups = "drop")
  my_sets <- setNames(sets_df$genes, sets_df$gs_name)
  
  genes_order <- gene_meta$gene_primary
  index <- lapply(my_sets, function(g) which(genes_order %in% g))
  cam_res <- cameraPR(stat = setNames(gene_meta$z, genes_order),
                      index = index, use.ranks = FALSE)
  
  list(
    gene_meta = gene_meta,  # mean_deltaM_gene, SE_gene, z, p_gene, n_cpg
    camera    = cam_res)
  
}

