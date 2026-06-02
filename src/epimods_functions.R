### EpiMod functions revised to current igraph package

# statLabel(g, stat) : returns weighted graph with edge weights obtained by averaging entries of stat along edges of iGraph g
# computeMods(g, alph = 0.05, gam = 0.5) : Returns list of vectors of ENTREZ IDs of spinglass modules
#                     found in graph g, with top alpha proportion of nodes (ranked by statistics)
#                     used as seeds and spinglass parameter gamma taken to be gam.
# nodeValidate(lv, g, tstats, nrandomizations = 1000) : Validate spinglass modules in lv on weighted iGraph object
#                     using provided tstats for node weights in validation. FDRs are estimated across
#                     nrandomizations number of randomizations. Returns list of length 3, 
#                     whose first element [["fdr"]] is the estimated false discovery rate.
# Note ignore warning, because newer implementation compatible with modern iGraph gives fundamentally different results
# renderModule(eid, g, pval, stat) : Renders module with ENTREZ IDs eid, on top of network g with 
#                     p-values given by pval and statistics stat.

statLabel <- function(g, stat) {
  
  # Ensure stat matches graph nodes
  stat <- stat[V(g)$name]
  
  # Convert graph to adjacency matrix
  A <- as_adjacency_matrix(g, sparse = FALSE)
  
  astat <- abs(stat)
  TMAX <- max(astat, na.rm = TRUE)
  
  # Compute weighted adjacency
  temp1 <- A * astat
  W <- (temp1 + t(temp1)) / (2 * TMAX)
  
  print("Generating weighted graph")
  
  # Build graph from weighted matrix
  G <- graph_from_adjacency_matrix(
    W,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # Assign node weights
  V(G)$weight <- astat
  V(G)$name <- rownames(W)
  
  return(G)
}

computeMods = function(g, alph = 0.05, gam = 0.5 ) {
  # Returns list of vectors of modules found by spinglass algorithm applied to
  #   edge weighted iGraph object g, using top alph nodes as seeds with spinglass gamma parameter = gam
  temp = V(g)$weight; names(temp) = V(g)$name;
  temp = sort(temp, decreasing=TRUE);
  v = names(temp[1:as.integer(alph*length(temp))]); # select the top ALPHA genes

  out.lv = list();
  for (j in 1:length(v)) {
    print(paste("Computing modules: Vertex",j,"of",length(v)));
    out.lv[[j]] = V(g)$name[spinglass.community(g, weights = E(g)$weight,
                                                vertex = v[[j]], gamma = gam)$community];
    # WARNING: Older versions of iGraoh need were 0-indexed rather than 1-indexed
    #          As such, they require an additional +1 here to work properly!
    # Check with: sessionInfo()$otherPkgs$igraph$Version
  }; names(out.lv) = v; # label with gene names
  return(out.lv);
 }

computeMods_new <- function(g, alph = 0.05, gam = 0.5) {
  
  temp <- V(g)$weight
  names(temp) <- V(g)$name
  temp <- sort(temp, decreasing = TRUE)
  
  seeds <- names(temp[1:as.integer(alph * length(temp))])
  
  out.lv <- list()
  
  for (j in seq_along(seeds)) {
    print(paste("Computing modules: Vertex", j, "of", length(seeds)))
    
    seed_vertex <- V(g)[name == seeds[j]]
    
    sg <- cluster_spinglass(
      g,
      weights = E(g)$weight,
      vertex = seed_vertex,
      gamma = gam
    )
    
    comm_id <- sg$membership[seed_vertex]
    
    out.lv[[j]] <- V(g)$name[sg$membership == comm_id]
  }
  
  names(out.lv) <- seeds
  return(out.lv)
}

# given list lv of vectors, return symmetric matrix of jaccard coefficients between elements of lv
list.jaccard <- function(lv) {
  sapply(lv, function(v)
    sapply(lv, function(w)
      length(intersect(v, w)) / length(union(v, w))
    )
  )
}

Modularity = function(v, g)
  # Returns modularity score of subgraph induced by vertices v in iGraph g
{ 
  h = induced_subgraph(g, v);
  return(sum( E(h)$weight ))
}

nodeValidate = function(lv, g, tstats, nrandomizations = 1000)
  # Test elements of lv in g over node-randomizations of tstats 
{
  #A = get.adjacency(g);
  nrm.lv = list(); 
  nrp.v = vector();
  obs.v = vector();
  output = list();
  TMAX = max(tstats);
  
  # First restrict to g:
  mods = lapply(lv, function(v) {
    return( intersect(v,V(g)$name) )
  });
  print(sapply(mods, length));
  
  # Observed modularity values
  obs.v = as.vector(lapply(mods, function (j) { return(Modularity(j, g)) }));
  
  # Node-randomizations of modularity
  nmods = length(mods);
  nrm.lv = lapply(1:nmods, function (i) {
    j = mods[[i]];
    v = vector();
    h = induced.subgraph(g,j); # new iGraph uses induced.subgraph() not subgraph()
    B = get.adjacency(h, sparse=FALSE); # Do not want sparse object
    for (k in 1:nrandomizations) {
      print(paste("Testing significance: module",i,"of",nmods,"Randomization",k,"of",nrandomizations))
      atperm = sample(tstats, nrow(B) , replace=FALSE)
      temp1 = apply(B, 1, function(v) return(v*atperm))
      W = (temp1 + t(temp1))/(2*TMAX);
      v[k] = sum(W)/2; # W is weighted adj matrix (with diag=0) so every edge counted twice
    }
    return(v);
  }); names(nrm.lv) = names(mods);
  
  # Empirical p-values
  for (j in 1:length(mods)) {
    nrp.v[j] = length( which(nrm.lv[[j]] > obs.v[j]) )/nrandomizations;
  }; names(nrp.v) = names(mods);
  
  output[[1]] = nrp.v; output[[2]] = obs.v; output[[3]] = nrm.lv; 
  names(output) = c("fdr","Observed","Random");
  return(output);
}

# Recluster spinglass
computeMods_recluster <- function(g, alph = 0.05, gam = 0.5) {
  # g: pruned, connected backbone graph
  # alph: top fraction of nodes to use as "seed focus" (for reporting, not algorithm)
  # gam: spinglass gamma
  
  temp <- V(g)$weight
  names(temp) <- V(g)$name
  temp <- sort(temp, decreasing = TRUE)
  
  seeds <- names(temp[1:as.integer(alph * length(temp))])
  
  out.lv <- list()
  
  # Run spinglass globally on the backbone
  sg <- cluster_spinglass(g, weights = E(g)$weight, gamma = gam)
  
  for (j in seq_along(seeds)) {
    seed <- seeds[j]
    # extract the cluster containing this seed
    cluster_members <- V(g)$name[sg$membership == sg$membership[V(g)[name == seed]]]
    
    # fallback for very small clusters
    if (length(cluster_members) < 3) {
      cluster_members <- names(temp)[1:min(5, length(temp))]
    }
    
    out.lv[[j]] <- cluster_members
  }
  
  names(out.lv) <- seeds
  return(out.lv)
}

# Extract core modules from large hairballs
hairball_filter <- function(mods, g, stat.v, 
                            size_thresh = 50, 
                            edge_quantile = 0.75, # remove weakest 25% edges
                            gam = 0.5) {
  out <- list()
  
  base_names <- names(mods)
  if (is.null(base_names)) base_names <- rep("", length(mods))
  
  for (m in seq_along(mods)) {
    v <- mods[[m]]
    mod_name <- base_names[m]
    
    # Keep small modules unchanged (append immediately)
    if (length(v) <= size_thresh) {
      if (nzchar(mod_name)) {
        out <- c(out, setNames(list(v), mod_name))
      } else {
        out <- c(out, list(v))
      }
      next
    }
    
    message(paste("Refining large module:",
                  if (nzchar(mod_name)) mod_name else paste0("mod", m),
                  "size =", length(v)))
    
    # Induce subgraph
    subg <- igraph::induced_subgraph(g, v)
    
    # --- STEP 1: Edge backbone (keep strongest edges) ---
    ew <- igraph::E(subg)$weight
    cutoff <- stats::quantile(ew, edge_quantile, na.rm = TRUE)
    subg_backbone <- igraph::delete_edges(subg, igraph::E(subg)[weight < cutoff])
    
    # Remove isolates
    subg_backbone <- igraph::delete_vertices(subg_backbone, igraph::degree(subg_backbone) == 0)
    
    # If graph collapses, skip (no fallback)
    if (igraph::vcount(subg_backbone) < 10) next
    
    # Keep largest connected component
    comp <- igraph::components(subg_backbone)
    if (length(comp$csize) > 1) {
      largest <- which.max(comp$csize)
      subg_backbone <- igraph::induced_subgraph(
        subg_backbone,
        igraph::V(subg_backbone)[comp$membership == largest]
      )
    }
    
    # If still too small, skip (no fallback)
    if (igraph::vcount(subg_backbone) < 10) next
    
    # --- STEP 2: re-cluster ---
    spin2.lv <- computeMods_recluster(subg_backbone, alph = 0.05, gam = gam)
    
    # Filter sizes
    if (length(spin2.lv) > 0) {
      sizes <- vapply(spin2.lv, length, integer(1))
      spin2.lv <- spin2.lv[sizes >= 8 & sizes <= 500]
    }
    
    # If nothing valid -> skip (no fallback)
    if (length(spin2.lv) == 0) next
    
    # Remove perfectly identical modules
    n <- length(spin2.lv)
    if (n >= 2) {
    jac.m <- list.jaccard(spin2.lv)
    n <- length(spin2.lv)
    keep <- rep(TRUE, n)
    for (k in seq_len(n)) {
      if (!keep[k]) next
      dup.idx <- which(jac.m[k, ] == 1 & seq_len(n) != k)
      if (length(dup.idx) > 0) keep[dup.idx] <- FALSE
    }
    spin2.lv <- spin2.lv[keep]}
    
    # If all removed as duplicates -> skip (no fallback)
    if (length(spin2.lv) == 0) next
    
    # --- STEP 3: append all valid submodules as-is (no renaming) ---
    out <- c(out, spin2.lv)
  }
  
  return(out)
}

# Name fixing
ensure_unique_mod_names <- function(lv, prefix = "mod") {

  firsts <- vapply(
    lv,
    function(x) {
      if (length(x) > 0 && !is.na(x[1]) && nzchar(as.character(x[1]))) {
        as.character(x[1])
      } else {
        NA_character_
      }
    },
    character(1)
  )
  
  # Ensure uniqueness by appending suffixes only when needed
  names(lv) <- make.unique(firsts, sep = "_")
  lv
}

# Wrapper to detect and validate modules
run_epimod <- function(pin, stat.v, pval.v){
  
  # Match to pin
  common_genes <- intersect(V(pin)$name, names(stat.v))
  
  pin_sub <- induced_subgraph(pin, common_genes)
  
  stat.v <- stat.v[common_genes]
  pval.v <- pval.v[common_genes]
  
  # create weighted graph
  weight.g <- statLabel(pin_sub, stat.v)
  
  # Keep only largest connected component
  components <- components(weight.g)
  
  largest_comp <- which.max(components$csize)
  
  weight.g <- induced_subgraph(
    weight.g,
    which(components$membership == largest_comp)
  )
  
  # Find modules
  spin.lv <- computeMods(weight.g, alph = 0.05, gam = 0.5)
  
  # Filter sizes
  sizes <- sapply(spin.lv, length)
  spin.lv <- spin.lv[sizes >= 8 & sizes <= 500]
  
  # Remove perfectly identical modules
  n <- length(spin.lv)
  if (n >= 2) {
    jac.m <- list.jaccard(spin.lv)
    if (!is.matrix(jac.m)) jac.m <- as.matrix(jac.m)
    keep <- rep(TRUE, n)
    for (k in seq_len(n)) {
      if (!keep[k]) next
      dup.idx <- which(jac.m[k, ] == 1 & seq_len(n) != k)
      if (length(dup.idx) > 0) keep[dup.idx] <- FALSE
    }
    spin.lv <- spin.lv[keep]
  }
 
  # Deal with large clusters/hairballs
  spin.lv <- hairball_filter(spin.lv, weight.g, stat.v)
  spin.lv <- ensure_unique_mod_names(spin.lv)
  
  # Validate
  validate.l <- nodeValidate(spin.lv, weight.g, abs(stat.v), nrandomizations = 1000)
  fdr <- validate.l$fdr
  
  # Prune near-contained/near-duplicates using validation as tie-breaker
  spin.lv <- prune_by_overlapcoef(spin.lv, fdr = fdr, g = weight.g, coef_thresh = 0.80)
  spin.lv <- ensure_unique_mod_names(spin.lv)  # keep names consistent after pruning
  
  # Subset fdr to the kept modules and re-order
  fdr <- fdr[names(spin.lv)]
  
  # Significant modules?
  sig.idx <- which(fdr <= 0.05)
  sig.mods <- spin.lv[sig.idx]
  fdr.sig <- fdr[sig.idx]
  ord <- order(fdr.sig)
  sig.mods <- sig.mods[ord]
  fdr.sig <- fdr.sig[ord] #13 out of 19 :D
  
  return(list(
    all.mods = spin.lv,
    fdr = fdr,
    sig.mods = sig.mods,
    sig.fdr = fdr.sig,
    weight.g = weight.g,
    pval.v = pval.v, # filtered
    stat.v = stat.v #filtered
  ))
  
}

# Prune overlapping modules
prune_by_overlapcoef <- function(lv, fdr = NULL, g = NULL, coef_thresh = 0.80) {
  overlap_coef <- function(a, b) {
    inter <- length(intersect(a, b))
    denom <- min(length(a), length(b))
    if (denom == 0) return(0)
    inter / denom
  }
  mod_score <- function(v) {
    if (is.null(g)) return(NA_real_)
    h <- igraph::induced_subgraph(g, vids = which(igraph::V(g)$name %in% v))
    if (igraph::ecount(h) == 0) return(0)
    sum(igraph::E(h)$weight, na.rm = TRUE)
  }
  nm <- names(lv)
  if (is.null(nm)) nm <- paste0("mod", seq_along(lv))
  names(lv) <- nm
  if (!is.null(fdr)) fdr <- fdr[nm]
  
  keep <- rep(TRUE, length(lv))
  for (i in seq_along(lv)) {
    if (!keep[i]) next
    for (j in seq_along(lv)) {
      if (i == j || !keep[j]) next
      coef <- overlap_coef(lv[[i]], lv[[j]])
      if (coef >= coef_thresh) {
        cand <- c(nm[i], nm[j])
        fdr_cmp <- if (!is.null(fdr)) c(fdr[cand[1]], fdr[cand[2]]) else c(NA_real_, NA_real_)
        mod_cmp <- c(mod_score(lv[[cand[1]]]), mod_score(lv[[cand[2]]]))
        size_cmp <- c(length(lv[[cand[1]]]), length(lv[[cand[2]]]))
        # Prefer: lower fdr (NA last), then higher mod score, then larger size
        ord <- order(is.na(fdr_cmp), fdr_cmp, -mod_cmp, -size_cmp)
        winner <- cand[ord][1]
        loser  <- setdiff(cand, winner)
        keep[match(loser, nm)] <- FALSE
      }
    }
  }
  lv[keep]
}

renderModule_gg <- function(
    eid, g, pval, stat,
    vertex_cols = c("yellow", "blue", mid = "white"),
    k = 50,
    id2symbol = NULL,
    use_symbols = FALSE,
    module_id = NULL,
    sl_limits = NULL,
    w_limits  = NULL,
    w_breaks  = NULL,
    edge_width_range = c(0.3, 2),
    edge_use_abs = TRUE,
    show_legends = TRUE
) {
  library(igraph)
  library(ggraph)
  library(ggplot2)
  
  h <- induced_subgraph(g, vids = which(V(g)$name %in% eid))
  
  stat.v   <- stat[V(h)$name]
  pval.v   <- pval[V(h)$name]
  slpval.v <- sign(stat.v) * -log10(pval.v)
  
  if (!is.null(id2symbol) && use_symbols) {
    labels <- id2symbol[V(h)$name]
    labels[is.na(labels)] <- V(h)$name
  } else {
    labels <- V(h)$name
  }
  main_id   <- if (!is.null(module_id)) module_id else eid[1]
  label_col <- ifelse(V(h)$name == main_id, "orchid4", "darkslategrey")
  
  w_raw <- E(h)$weight
  w_map <- if (edge_use_abs) abs(w_raw) else w_raw
  
  if (is.null(sl_limits)) {
    M <- max(abs(slpval.v), na.rm = TRUE)
    sl_limits <- c(-M, M)
  }
  if (is.null(w_limits)) {
    w_limits <- range(w_map, na.rm = TRUE)
  }
  
  if (is.null(w_breaks)) {
    w_breaks <- pretty(w_limits, n = 3)
    w_breaks <- unique(w_breaks[w_breaks >= min(w_limits) & w_breaks <= max(w_limits)])
    if (length(w_breaks) < 2) {
      w_breaks <- unique(c(min(w_limits), mean(w_limits), max(w_limits)))
    }
  }
  
  n <- gorder(h)
  vertex_size <- if (n < 50) 5 else max(2, 5 * 50 / n)
  label_size  <- if (n < 50) 5 else max(1.5, 5 * 50 / n)
  
  p <- ggraph(h, layout = "fr") +
    geom_edge_link(
      aes(width = w_map, colour = w_map),
      alpha = 0.35,
      linetype = 1,
      show.legend = TRUE,
      key_glyph = ggplot2::draw_key_path
    ) +
    geom_node_point(
      aes(fill = slpval.v),
      size = vertex_size,
      shape = 21,          
      colour = "black",    
      stroke = 0.3,       
      show.legend = TRUE
    ) +
    geom_node_text(
      aes(label = labels, colour = I(label_col)),
      size = label_size,
      repel = TRUE,
      fontface = "italic",
      show.legend = FALSE
    ) +
    scale_edge_colour_gradient(
      name    = "Edge weight",
      low     = "grey80",
      high    = "black",
      limits  = w_limits,
      breaks  = w_breaks,
      labels = scales::label_number(accuracy = 0.01),
      guide   = guide_legend(
        order = 2,
        override.aes = list(
          shape = NA, # should remove point glyph
          linetype = 1,
          linewidth = scales::rescale(w_breaks, to = edge_width_range),
          colour = scales::col_numeric(c("grey80", "black"), w_limits)(w_breaks),
          alpha = 1
        )
      )
    ) +
    # disable width legend
    scale_edge_width(
      limits = w_limits,
      range  = edge_width_range,
      guide  = "none"
    ) +
    scale_fill_gradient2(
      name     = "Signed -log10(P value)",
      low      = vertex_cols[1],
      mid      = vertex_cols[3],
      high     = vertex_cols[2],
      midpoint = 0,
      limits   = sl_limits,
      guide    = guide_colorbar(order = 1)
    ) +
    theme_void() +
    theme(
      legend.position = if (show_legends) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
  
  return(p)
}

# Helper to plot gene count and FDR for the detected modules
plot_fdr_genecount <- function(out_mod, id2symbol, gtit) {
  
  library(tidyverse)
  library(patchwork)
  
  df <- data.frame(
    ENTREZID = names(out_mod[[2]]),
    FDR = unname(out_mod[[2]]))
  df$gene_primary <- id2symbol[df$ENTREZID]
  gene_counts <- sapply(out_mod[[1]], length) 
  df_counts <- data.frame(
    ENTREZID = names(gene_counts),
    gene_count = unname(gene_counts)
  )
  df <- df %>% dplyr::left_join(df_counts)
  df$sig <- as.character(ifelse(df$FDR < 0.05, 1, 0))
  
  pdat <- df %>%
    mutate(logFDR = -log10(FDR)) %>%
    arrange(logFDR)
  
  p_left <- ggplot(pdat, aes(x = gene_count, y = reorder(gene_primary, logFDR), fill = sig)) +
    geom_col() +
    scale_fill_manual(values = c("azure3", "darkslategrey")) +
    scale_x_reverse() +  # makes bars go to the left
    labs(x = "Gene count", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y.left = element_blank(),
          legend.position = "none") +
    ggtitle(gtit)
  p_mid <- ggplot(pdat, aes(y = reorder(gene_primary, logFDR), x = 1)) +
    geom_text(aes(label = gene_primary), fontface = "italic") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    ggtitle("")
  p_right <- ggplot(pdat, aes(x = logFDR, y = reorder(gene_primary, logFDR), fill = sig)) +
    geom_col() +
    scale_fill_manual(values = c("azure3", "darkslategrey")) +
    labs(x = "-log10(FDR)", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y.left = element_blank(),
          legend.position = "none") +
    ggtitle("")
  p <- p_left + p_mid + p_right + plot_layout(widths = c(2, 0.35, 2))
  
  return(p)
  
}
