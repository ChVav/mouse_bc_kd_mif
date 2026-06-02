
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggtext)

# depends on helpers_bipartite_network.R
# Would need some modification for changing fill/size parameters

wrap_bipartite_networks <- function(gene_meta, gene_effects, prob_thres = 0.7, min_genes, min_genes2, region_var, region_value, label_thres,
                                    pathway_size, fill_var2, size_var2, pathway_size2, layout, fill_name2, size_name2){
  
  ## Baseline #----
  prep <- prepare_bipartite_data(
    gene_meta,
    region_var = region_var,
    region_value = region_value,
    label_thres = label_thres,
    min_genes = min_genes 
  )
  
  prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
  g <- build_bipartite_graph(
    prep$edges,
    fill_var = "effect",
    size_var = "neglog10_p",
    label_var = "effect",
    label_threshold = prep$threshold,
    pathway_size = pathway_size
  )
  
  # Keep gene labels to label modifier plots consistently
  prev_labels <- igraph::vertex_attr(g, "label")
  gene_vertices <- igraph::vertex_attr(g, "type") == "gene" # exclude pathways
  label_genes <- prev_labels[gene_vertices]
  label_genes <- label_genes[!is.na(label_genes)]  # keep only labeled genes
  
  ## Ketogenic diet #----
  
  gene_effects_mod <- gene_effects %>% 
    dplyr::filter(promoter == region_value & effect_type == "KD")
  
  # Start of with major hubs affected by pellets 
  edges2_kd <- prep$edges %>% dplyr::select(gene_primary, gs_name) 
  edges2_kd <- dplyr::left_join(edges2_kd, gene_effects_mod) 
  
  # only edges for which at least x genes have a |mean P| > prob_thres
  # Count number of strongly affected genes per pathway
  pathway_counts <- edges2_kd %>%
    filter(abs(mean_prob) > prob_thres) %>%   # only strong genes
    group_by(gs_name) %>%
    summarise(n_strong_genes = n_distinct(gene_primary))  # count unique genes
  n <- min_genes2 
  edges2_kd <- edges2_kd %>%
    filter(gs_name %in% pathway_counts$gs_name[pathway_counts$n_strong_genes > n])
  
  #length(unique(edges2_kd$gs_name))
  
  g2 <- build_bipartite_graph(
    edges2_kd,
    fill_var = fill_var2,
    size_var = size_var2,
    label_genes = label_genes,
    label_mode = "include",
    pathway_size = pathway_size2
  )
  
  ## Mifepristone #----
  
  gene_effects_mod <- gene_effects %>% 
    dplyr::filter(promoter == region_value & effect_type == "MIF")
  
  # Start of with major hubs effected by pellets (40 with n >=10 genes at baseline)
  edges2_mf <- prep$edges %>% dplyr::select(gene_primary, gs_name) 
  edges2_mf <- dplyr::left_join(edges2_mf, gene_effects_mod) 
  
  # only edges for which at least x genes have a |mean P| > prob_thres
  # Count number of strongly affected genes per pathway
  pathway_counts <- edges2_mf %>%
    filter(abs(mean_prob) > prob_thres) %>%   # only strong genes
    group_by(gs_name) %>%
    summarise(n_strong_genes = n_distinct(gene_primary))  # count unique genes
  n <- min_genes2 
  edges2_mf <- edges2_mf %>%
    filter(gs_name %in% pathway_counts$gs_name[pathway_counts$n_strong_genes > n])
  
  g3 <- build_bipartite_graph(
    edges2_mf,
    fill_var = fill_var2,
    size_var = size_var2,
    label_genes = label_genes,
    label_mode = "include",
    pathway_size = pathway_size2
  )
  
  ## Ketogenic diet + Mifepristone #----
  
  gene_effects_mod <- gene_effects %>% 
    dplyr::filter(promoter == region_value & effect_type == "KD_MIF")
  
  # Start of with major hubs affected by pellets 
  edges2_kdmf <- prep$edges %>% dplyr::select(gene_primary, gs_name) 
  edges2_kdmf <- dplyr::left_join(edges2_kdmf, gene_effects_mod) 
  
  # only edges for which at least x genes have a |mean P| > prob_thres
  # Count number of strongly affected genes per pathway
  pathway_counts <- edges2_kdmf %>%
    filter(abs(mean_prob) > prob_thres) %>%   # only strong genes
    group_by(gs_name) %>%
    summarise(n_strong_genes = n_distinct(gene_primary))  # count unique genes
  n <- min_genes2 
  edges2_kdmf <- edges2_kdmf %>%
    filter(gs_name %in% pathway_counts$gs_name[pathway_counts$n_strong_genes > n])
  
  length(unique(edges2_kdmf$gs_name)) #124
  
  g4 <- build_bipartite_graph(
    edges2_kdmf,
    fill_var = fill_var2,
    size_var = size_var2,
    label_genes = label_genes,
    label_mode = "include",
    pathway_size = pathway_size2
  )
  
  ## Make plot #----
  # Do this separately, so that can merge legend using global scales
  
  global_fill_limits <- range(
    c(edges2_kd$mean_delta,
      edges2_mf$mean_delta,
      edges2_kdmf$mean_delta),
    na.rm = TRUE
  )
  
  p1 <- plot_bipartite_graph(
    g,
    fill_var = "fill_value",
    fill_name = "Mean Δ (P/D+ − P/D−)",
    size_name = expression(-Log[10](adjusted~p)),
    size_labels = scales::number_format(accuracy = 1),
    layout = layout
  ) + ggtitle("No modifier (baseline)") +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    )
  
  # dh/stress > kk > fr
  
  p2 <- plot_bipartite_graph(
    g2,
    fill_var = "fill_value",
    fill_low = "#6a00a8",
    fill_high = "#ff7f0e",
    layout = layout
  ) + ggtitle("KD") +
    scale_fill_gradient2(
      low = "#6a00a8",
      mid = "white",
      high = "#ff7f0e",
      limits = global_fill_limits,
      name = fill_name2
    ) +
    scale_size_continuous(
      limits = c(0,1),
      breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1),
      trans = scales::trans_new("stretch", transform = function(x) x^7, inverse = function(x) x^(1/7)),
      name = size_name2
    ) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    )
    
  p3 <- plot_bipartite_graph(
    g3,
    fill_var = "fill_value",
    fill_low = "#6a00a8",
    fill_high = "#ff7f0e",
    layout = layout
  ) + ggtitle("MIF") +
    scale_fill_gradient2(
      low = "#6a00a8",
      mid = "white",
      high = "#ff7f0e",
      limits = global_fill_limits,
      name = fill_name2
    ) +
    scale_size_continuous(
      limits = c(0,1),
      breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1),
      trans = scales::trans_new("stretch", transform = function(x) x^7, inverse = function(x) x^(1/7)),
      name = size_name2
    ) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    )
  
  p4 <- plot_bipartite_graph(
    g4,
    fill_var = "fill_value",
    fill_low = "#6a00a8",
    fill_high = "#ff7f0e",
    layout = layout
  ) + ggtitle("KD + MIF") +
    scale_fill_gradient2(
      low = "#6a00a8",
      mid = "white",
      high = "#ff7f0e",
      limits = global_fill_limits,
      name = fill_name2
    ) +
    scale_size_continuous(
      limits = c(0,1),
      breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1),
      trans = scales::trans_new("stretch", transform = function(x) x^7, inverse = function(x) x^(1/7)),
      name = size_name2
    ) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    )
  
  return(list(
    edge_base = prep$eges,
    edge_kd = edges2_kd,
    edge_mf = edges2_mf,
    edge_kdmf = edges2_kdmf,
    g_base = g,
    g_kd = g2,
    g_mf = g3,
    g_kdmf = g4,
    p_base = p1,
    p_kd = p2,
    p_mf = p3,
    p_kdmf = p4,
    params = as.list(environment())
  ))
  
}


