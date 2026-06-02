library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggtext)

# # Example
# prep <- prepare_bipartite_data(
#   gene_meta,
#   region_var = "promoter",
#   region_value = TRUE,
#   min_genes = 8
# )
# 
# prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
# g <- build_bipartite_graph(
#   prep$edges,
#   fill_var = "effect",
#   size_var = "neglog10_p",
#   label_var = "effect",
#   label_threshold = prep$threshold,
#   pathway_size = 10
# )
# 
# plot_bipartite_graph(
#   g,
#   fill_var = "fill_value",
#   fill_name = "Mean Δ (P/D+ − P/D−)",
#   size_name = expression(-Log[10](adjusted~p)),
#   size_labels = scales::number_format(accuracy = 1),
#   layout = "dh"
# )

prepare_bipartite_data <- function(gene_meta,
                                   region_var = "promoter",
                                   region_value = TRUE,
                                   p_cutoff = 0.05,
                                   min_genes = 5,
                                   label_thres = 0.75,
                                   effect_var = "mean_deltaM_gene") {
  
  # Filter region + significance
  df <- gene_meta %>%
    dplyr::filter(.data[[region_var]] == region_value) %>%
    dplyr::filter(p_adj_gene < p_cutoff)
  
  # Label threshold (top 100-x% effect magnitude)
  threshold <- quantile(abs(dplyr::pull(df, .data[[effect_var]])), label_thres, na.rm = TRUE)
  
  edges <- df %>%
    dplyr::filter(!is.na(gs_name)) %>%
    dplyr::select(gene_primary, gs_name,
                  effect = .data[[effect_var]],
                  p_adj_gene) %>%
    dplyr::mutate(label_name =
                    ifelse(abs(effect) > threshold,
                           gene_primary, NA))
  
  # Count genes per pathway
  pathway_counts <- edges %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")
  
  # Filter pathways
  edges <- edges %>%
    dplyr::filter(gs_name %in%
                    pathway_counts$gs_name[pathway_counts$n_genes > min_genes])
  
  list(edges = edges,
       threshold = threshold)
}

build_bipartite_graph <- function(edges,
                                  gene_col = "gene_primary",
                                  pathway_col = "gs_name",
                                  fill_var = NULL,
                                  size_var = NULL,
                                  label_var = NULL,
                                  label_genes = NULL,        # character vector of genes
                                  label_threshold = NULL,
                                  label_mode = c("threshold", "include", "exclude"), # Default threshold
                                  pathway_size = 10,
                                  up_pathways = NULL, 
                                  down_pathways = NULL) {
  
  label_mode <- match.arg(label_mode)
  
  # Build graph
  g <- igraph::graph_from_data_frame(
    edges[, c(gene_col, pathway_col)],
    directed = FALSE
  )
  
  # Node type
  V(g)$type <- ifelse(
    V(g)$name %in% edges[[pathway_col]],
    "pathway", "gene"
  )
  
  # Unique gene-level data
  gene_df <- edges %>%
    dplyr::distinct(.data[[gene_col]], .keep_all = TRUE)
  
  # Attach fill variable
  if (!is.null(fill_var)) {
    gene_vertices <- which(V(g)$type == "gene")
    matched <- match(V(g)$name[gene_vertices], gene_df[[gene_col]])
    V(g)$fill_value <- NA
    V(g)$fill_value[gene_vertices] <- gene_df[[fill_var]][matched]
  }
  
  # Attach size variable
  if (!is.null(size_var)) {
    gene_vertices <- which(V(g)$type == "gene")
    matched <- match(V(g)$name[gene_vertices], gene_df[[gene_col]])
    V(g)$size <- NA
    V(g)$size[gene_vertices] <- gene_df[[size_var]][matched]
  } else {
    V(g)$size <- NA
  }
  
  # Pathway fixed size
  V(g)$size[V(g)$type == "pathway"] <- pathway_size
  
  # Labels
  V(g)$label <- NA
  
  gene_vertices <- which(V(g)$type == "gene")
  gene_names <- V(g)$name[gene_vertices]
  
  if (label_mode == "threshold" && 
      !is.null(label_var) && 
      !is.null(label_threshold)) {
    
    matched <- match(gene_names, gene_df[[gene_col]])
    
    V(g)$label[gene_vertices] <- ifelse(
      abs(gene_df[[label_var]][matched]) > label_threshold,
      gene_df[[gene_col]][matched],
      NA
    )
  }
  
  if (label_mode == "include" && !is.null(label_genes)) {
    V(g)$label[gene_vertices] <- ifelse(
      gene_names %in% label_genes,
      gene_names,
      NA
    )
  }
  
  if (label_mode == "exclude" && !is.null(label_genes)) {
    V(g)$label[gene_vertices] <- ifelse(
      !gene_names %in% label_genes,
      gene_names,
      NA
    )
  }
  
  # Pathways always labeled
  V(g)$label[V(g)$type == "pathway"] <- V(g)$name[V(g)$type == "pathway"]
  
  # Label formatting
  V(g)$label_expr <- ifelse(
    V(g)$type == "pathway",
    paste0("bold('", V(g)$label, "')"),
    ifelse(!is.na(V(g)$label),
           paste0("italic('", V(g)$label, "')"),
           NA)
  )
  
  # Optional: pathway direction tagging
  V(g)$pathway_dir <- NA_character_
  if (!is.null(up_pathways)) {
    V(g)$pathway_dir[V(g)$type == "pathway" & V(g)$name %in% up_pathways] <- "up"}
  if (!is.null(down_pathways)) {
    V(g)$pathway_dir[V(g)$type == "pathway" & V(g)$name %in% down_pathways] <- "down"}
  
  return(g)
}

plot_bipartite_graph <- function(g,
                                 fill_var = "effect",
                                 fill_low = "#4575b4",
                                 fill_mid = "white",
                                 fill_high = "#d73027",
                                 fill_midpoint = 0,
                                 fill_name = "Estimate",
                                 size_name = "Size",
                                 size_labels = waiver(),
                                 pathway_size = 10,
                                 layout = "fr",
                                 seed = 123,
                                 pathway_colors = c(HypoM = "#4575b4", HyperM = "#d73027"),
                                 pathway_legend = FALSE,
                                 up_pathways = NULL, 
                                 down_pathways = NULL) {
  
  set.seed(seed)
  
  # Extract gene-only size values for scaling
  gene_sizes <- igraph::vertex_attr(g, "size")[
    igraph::vertex_attr(g, "type") == "gene"
  ]
  size_limits <- range(gene_sizes, na.rm = TRUE)
  
  # Pathway directions (levels must match names(pathway_colors))
  v_names <- igraph::V(g)$name
  v_type  <- igraph::vertex_attr(g, "type")
  path_idx   <- which(v_type == "pathway")
  path_names <- v_names[path_idx]
  
  pathway_dir_chr <- rep(NA_character_, length(v_names))
  if (!is.null(up_pathways)) {
    pathway_dir_chr[path_idx[path_names %in% up_pathways]] <- "HyperM"
  }
  if (!is.null(down_pathways)) {
    pathway_dir_chr[path_idx[path_names %in% down_pathways]] <- "HypoM"
  }
  pathway_levels <- names(pathway_colors)
  pathway_dir_fac <- factor(pathway_dir_chr, levels = pathway_levels)
  g <- igraph::set_vertex_attr(g, "pathway_dir", value = pathway_dir_fac)
  
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(color = "lightgray") +
    
    # Gene nodes (mapped to size + continuous fill)
    ggraph::geom_node_point(
      data = function(x) dplyr::filter(x, type == "gene"),
      aes(size = size, fill = .data[[fill_var]]),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_fill_gradient2(
      name = fill_name,
      low = fill_low,
      mid = fill_mid,
      high = fill_high,
      midpoint = fill_midpoint,
      na.value = "white",
      guide = guide_colorbar(order = 1)
    ) +
    # Make the size legend use hollow circles and a controlled size range
    scale_size_continuous(
      name   = size_name,
      limits = size_limits,
      labels = size_labels,
      guide  = guide_legend(
        order = 2,
        override.aes = list(
          fill   = NA,
          colour = "black",
          shape  = 21
        )
      )
    )
  
  # Pathway nodes:
  any_dir <- any(!is.na(pathway_dir_fac[path_idx]))
  if (any_dir) {
    p <- p +
      ggnewscale::new_scale_fill() +
      ggraph::geom_node_point(
        data = function(x) dplyr::filter(x, type == "pathway"),
        aes(fill = pathway_dir),
        size = pathway_size,
        shape = 21,
        color = "black",
        alpha = 0.25,
        show.legend = c(fill = pathway_legend, size = FALSE)
      ) +
      scale_fill_manual(
        name   = "Pathway",
        values = pathway_colors,
        limits = pathway_levels,
        drop   = FALSE,
        na.value = "white",
        guide  = guide_legend(order = 3,
                              override.aes = list(shape = 21, colour = "black"))
      )
  } else {
    p <- p +
      ggraph::geom_node_point(
        data = function(x) dplyr::filter(x, type == "pathway"),
        size = pathway_size,
        fill = "white",
        shape = 21,
        color = "black",
        show.legend = FALSE
      )
  }
  
  # Draw text last so it sits on top of nodes
  p <- p +
    ggraph::geom_node_text(
      aes(label = label_expr),
      repel = TRUE,
      parse = TRUE,
      size = 3
    ) +
    theme_void()
  
  return(p)
}
