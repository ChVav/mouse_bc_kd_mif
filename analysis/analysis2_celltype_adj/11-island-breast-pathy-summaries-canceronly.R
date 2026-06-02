## All CpGs limma estimates were summarized using inverse-Variance Weighted gene-level meta-analysis (see 11)
## Significant pathways of each modifier compared to baseline were detected using competitive gene testing (camera, see 12 and reactome results dag-aware trimmed based on overlap see 13)

## Per database, for all significant pathways plot bipartite network signatures - genes; promoter and gene body regions separately
# Use a gene threshold?
# Alternatively dot plots?
## Per database, plot pathway-level significance (perhaps start next script; use again analytical posterior framework)

library(here)
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggtext)
library(patchwork)
library(showtext) # so arrows and delta are rendered properly

showtext_auto()

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/11-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_pathway_summary.R"))
source(here("src/helpers_bipartite_network.R"))
source(here("src/wrap_bipartite_networks.R"))

gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds")) %>%
  dplyr::select(gene_primary, gs_name)

# MIF #-----

## Promoter #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/mif_reactome_promoter_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks not consistent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (MIF − Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks inconsitent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"mif_breast_reactome_camera_promoter.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)


## Gene body #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/mif_reactome_genebody_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.025) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"mif_breast_reactome_camera_genebody.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)

# KD #-----

## Promoter #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kd_reactome_promoter_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks not consistent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD − Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks inconsitent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kd_breast_reactome_camera_promoter.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)


## Gene body #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kd_reactome_genebody_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.025) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kd_breast_reactome_camera_genebody.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)

# KD + MIF #-----

## Promoter #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kdmif_reactome_promoter_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks not consistent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF − Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks inconsitent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kdmif_breast_reactome_camera_promoter.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)


## Gene body #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kdmif_reactome_genebody_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.025) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kdmif_breast_reactome_camera_genebody.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)

# KD + MIF vs MIF #-----

## Promoter #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kdmif_vs_mif_reactome_promoter_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks not consistent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF − Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks inconsitent
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = TRUE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Promoter, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kdmif_vs_mif_breast_reactome_camera_promoter.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)


## Gene body #-----
load(file.path(results_dir,"analysis2_celltype_adj/10-output/kdmif_vs_mif_reactome_genebody_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p1 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hyperM")

### Hypomethylated #-----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.025) %>% dplyr::filter(Direction == "Down") %>% rownames()

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
prep <- prepare_bipartite_data(
  gene_meta2,
  region_var = "promoter",
  region_value = FALSE,
  min_genes = 0,
  p_cutoff = 1,
  label_thres = 0.75
)

prep$edges$neglog10_p <- -log10(prep$edges$p_adj_gene)
g <- build_bipartite_graph(
  prep$edges,
  fill_var = "effect",
  size_var = "neglog10_p",
  label_var = "effect",
  label_threshold = prep$threshold,
  pathway_size = 10
)

p2 <- plot_bipartite_graph(
  g,
  fill_var = "fill_value",
  fill_name = "Mean Δ (KD + MIF - Baseline)",
  size_name = expression(-Log[10](adjusted~p)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("Gene body, hypoM") 

p <- p1 / p2 +
  plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"kdmif_vs_mif_breast_reactome_camera_genebody.pdf"),
       width = 183*2,
       height = 200*2,
       units = "mm"
)



