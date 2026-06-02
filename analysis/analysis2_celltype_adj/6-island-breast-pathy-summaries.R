
# Focus on Reactome
## All CpGs were summarized using inverse-Variance Weighted gene-level meta-analysis
## Significant pathways at baseline were detected using competitive gene testing and dag-aware trimmed

## for paths with p < 0.01 plot bipartite network signatures - genes; promoter and gene body regions separately
## Group significant reactome pathways into themes
## test both paths with p < 0.05 and 0.01 compute pathway scores (for limma/analytical posterior testing)

library(here)
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggtext)
library(patchwork)
library(stringr)
library(purrr)
library(showtext) # so arrows and delta are rendered properly
showtext_auto()
library(plotly)
library(singscore)
library(GSEABase)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/6-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}
dirOut2 <- file.path(dirOut,"p0.05")
if(!dir.exists(dirOut2)){dir.create(dirOut2, recursive = TRUE)}

source(here("src/helpers_pathway_summary.R"))
source(here("src/helpers_bipartite_network.R"))
source(here("src/wrap_bipartite_networks.R"))

# Reactome #-----
gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds")) %>%
  dplyr::select(gene_primary, gs_name)

## Promoter #-----
load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
#sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames() #51
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Up") %>% rownames() #26
#sig <- out$camera_pruned %>% dplyr::filter(FDR < 0.25) %>% dplyr::filter(Direction == "Up") %>% rownames() #11

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
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
  fill_name = "Mean ΔExposure",
  size_name = expression(-Log[10](adj.~P)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("pCGI, hyperM")

### Hypomethylated #-----
#sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames() #371
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Down") %>% rownames() #4
#sig <- out$camera_pruned %>% dplyr::filter(FDR < 0.25) %>% dplyr::filter(Direction == "Down") %>% rownames() # none

gene_meta2 <- gene_meta %>% dplyr::filter(gs_name %in% sig)
gene_meta2$gs_name <- normalize_name(gene_meta2$gs_name)

# Looks cool
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
  fill_name = "Mean ΔExposure",
  size_name = expression(-Log[10](adj.~P)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("pCGI, hypoM") 

des <- c("AA
          AA
          AA
          BC")

p <- p1 + p2 + plot_spacer() +
  plot_annotation(tag_levels = "a") +
  plot_layout(design = des)

ggsave(p,
       file = file.path(dirOut,"breast_reactome_camera_promoter.pdf"),
       width = 183*2,
       height = 150*2,
       units = "mm"
)

ggsave(p,
       file = file.path(dirOut,"breast_reactome_camera_promoter.png"),
       width = 183*2,
       height = 150*2,
       units = "mm"
)

## Gene body #-----
load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"))
gene_meta <- out$gene_meta
gene_meta$p_adj_gene <- p.adjust(gene_meta$p_gene, method = "BH")
gene_meta <- gene_meta %>% dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many")

### Hypermethylated #-----
# sig <- out$camera_pruned %>% dplyr::filter(FDR < 0.25) %>% dplyr::filter(Direction == "Up") %>% rownames() #1
# sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames() #56
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Up") %>% rownames() #17

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
  fill_name = "Mean ΔExposure",
  size_name = expression(-Log[10](adj.~P)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("oCGI, hyperM")

### Hypomethylated #-----
# sig <- out$camera_pruned %>% dplyr::filter(FDR < 0.25) %>% dplyr::filter(Direction == "Down") %>% rownames() # none
# sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames() #43
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Down") %>% rownames() #9

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
  fill_name = "Mean ΔExposure",
  size_name = expression(-Log[10](adj.~P)),
  size_labels = scales::number_format(accuracy = 1),
  layout = "dh"
) + ggtitle("oCGI, hypoM") 

des <- c("A
          B")

p <- p1 + p2 + 
  plot_annotation(tag_levels = "a") +
  plot_layout(design = des)

ggsave(p,
       file = file.path(dirOut,"breast_reactome_camera_genebody.pdf"),
       width = 183*2,
       height = 150*2,
       units = "mm"
)

ggsave(p,
       file = file.path(dirOut,"breast_reactome_camera_genebody.png"),
       width = 183*2,
       height = 150*2,
       units = "mm"
)

# Group and summarize reactome pathways #----
# Focus on p < 0.0.1

load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"))
prom_camera_pruned <- out$camera_pruned
load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"))
body_camera_pruned <- out$camera_pruned

sig_prom <- prom_camera_pruned %>% dplyr::filter(PValue < 0.01) %>% rownames()
sig_body <- body_camera_pruned %>% dplyr::filter(PValue < 0.01) %>% rownames()
sig <- unique(c(sig_prom,sig_body))

# sig <- prom_camera_pruned %>% dplyr::filter(PValue < 0.05) %>% rownames()
# sig2 <- body_camera_pruned %>% dplyr::filter(PValue < 0.05) %>% rownames()
# sig <- unique(c(sig,sig2))  #166 terms; let's narrow down later
# rm(sig2);gc()

# sig <- prom_camera_pruned %>% dplyr::filter(FDR < 0.25) %>% rownames()
# sig2 <- body_camera_pruned %>% dplyr::filter(FDR < 0.25) %>% rownames()
# sig <- unique(c(sig,sig2))  #17 terms; let's narrow down later
# rm(sig2);gc()

df <- data.frame(
  gs_name = sig,
  gs_name_norm = normalize_name(sig))

df_thematic <- add_pathway_theme(df, name_col = "gs_name_norm") 
saveRDS(df_thematic, file = file.path(dirOut, "reactome_themes.Rds"))

# # Check which genes fall ounder physiological factors
# df_thematic %>% dplyr::filter(is.na(theme_primary))
# body_camera_pruned %>% dplyr::filter(rownames(.) == "REACTOME_PHYSIOLOGICAL_FACTORS")
# load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"))
# gene_meta <- out$gene_meta
# gene_meta <- gene_meta %>% 
#   dplyr::left_join(gene_reactome_tbl, relationship = "many-to-many") %>%
#   dplyr::filter(gs_name == "REACTOME_PHYSIOLOGICAL_FACTORS") %>%
#   dplyr::filter(p_gene < 0.05)
# # In this case more or less hormone receptors

# Save also themes for supplementary info
df_thematic <- df_thematic %>% dplyr::select(theme_primary,gs_name_norm, gs_name) %>% dplyr::arrange(theme_primary)
colnames(df_thematic) <- c("Functional.category","Pathway","Reactome.gs.name")
df_thematic %>% write.csv(file.path(dirOut,"reactome_themes.csv"), row.names = FALSE)

# Theme summary, number of pathways and gene region
df_prom <- data.frame(gs_name_norm = normalize_name(sig_prom), gene_region = "Promoter") %>% add_pathway_theme(., name_col = "gs_name_norm")
df_body <- data.frame(gs_name_norm = normalize_name(sig_body), gene_region = "Gene body") %>% add_pathway_theme(., name_col = "gs_name_norm")
df <- bind_rows(df_prom, df_body)
saveRDS(df, file = file.path(dirOut, "themes_paths_0.01.Rds"))

# Compute pathway scores p <0.01 #----

anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds"))
cpg_island <- anno %>% dplyr::filter(Relation_to_Island == "Island") %>% pull(Name) # 29,022 islands

load('data/beta_final.Rdata')
load("data/pheno.Rdata")

pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast'),])

## Islands promoters #----

load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"))
promotor_cg <- anno %>% dplyr::filter(promoter == TRUE) %>% pull(Name)

# Create sample-level gene summary
beta_sub <- beta_final[rownames(beta_final) %in% intersect(cpg_island,promotor_cg),
                       pheno_sub$basename,
                       drop = FALSE]
cpg_annot <- anno %>%
  dplyr::filter(Name %in% rownames(beta_sub)) %>%
  dplyr::select(Name, gene_primary)

gene_meth <- as.data.frame(beta_sub) %>%
  mutate(
    across(
      where(is.numeric),
      ~ {x <- .
        log2(x / (1 - x))
      }
    )
  )
gene_meth$cpg_id <- rownames(beta_sub) #14645 islands in promoter region
gene_meth <- gene_meth %>%
  left_join(cpg_annot, by = c("cpg_id" = "Name")) %>%
  dplyr::filter(!is.na(gene_primary)) %>%
  group_by(gene_primary) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            n_cpg = n(),
            .groups = "drop") #7675 genes summarized
saveRDS(gene_meth, file = file.path(dirOut, "island_breast_genemeanmeth_promoters.Rds"))


mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

### Up #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Up") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut, "breast_reactome_promoters_sigup_singscore.Rdata"))

# Check dispersion, quite some variability here; clustering by pathway not phenotype
plot(scoredf$Scores, scoredf$Dispersion)

df_long <- as.data.frame(scoredf$Scores) %>%
  mutate(pathway = rownames(scoredf$Scores)) %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "score") %>%
  left_join(
    as.data.frame(scoredf$Dispersions) %>%
      mutate(pathway = rownames(scoredf$Dispersions)) %>%
      pivot_longer(-pathway, names_to = "sample", values_to = "dispersion"),
    by = c("pathway", "sample")
  )

plot_ly(
  df_long,
  x = ~score,
  y = ~dispersion,
  color = ~pathway,
  type = "scatter",
  mode = "markers",
  text = ~paste("Pathway:", pathway,
                "<br>Sample:", sample,
                "<br>Score:", round(score, 3),
                "<br>Dispersion:", round(dispersion, 1)),
  hoverinfo = "text"
)

### Down #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Down") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

# flip direction for ranking hypoM (negative M-values) first
mat <- mat * -1
ranked <- rankGenes(mat) 
scoredf <- multiScore(ranked, upSetColc = gsc)

# reverse again scores
scoredf$Scores <- scoredf$Scores * -1

save(scoredf, file = file.path(dirOut, "breast_reactome_promoters_sigdown_singscore.Rdata"))

## Islands gene body #----

load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"))
body_cg <- anno %>% dplyr::filter(promoter == FALSE) %>% pull(Name)

# Create sample-level gene summary
beta_sub <- beta_final[rownames(beta_final) %in% intersect(cpg_island,body_cg),
                       pheno_sub$basename,
                       drop = FALSE]
cpg_annot <- anno %>%
  dplyr::filter(Name %in% rownames(beta_sub)) %>%
  dplyr::select(Name, gene_primary)

gene_meth <- as.data.frame(beta_sub) %>%
  mutate(
    across(
      where(is.numeric),
      ~ {x <- .
      log2(x / (1 - x))
      }
    )
  )
gene_meth$cpg_id <- rownames(beta_sub) #14377 islands in gene body
gene_meth <- gene_meth %>%
  left_join(cpg_annot, by = c("cpg_id" = "Name")) %>%
  dplyr::filter(!is.na(gene_primary)) %>%
  group_by(gene_primary) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            n_cpg = n(),
            .groups = "drop") #6575 genes summarized
saveRDS(gene_meth, file = file.path(dirOut, "island_breast_genemeanmeth_genebody.Rds"))


mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

### Up #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Up") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut, "breast_reactome_genebody_sigup_singscore.Rdata"))

### Down #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == "Down") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

# flip direction
mat <- mat * -1
ranked <- rankGenes(mat) 
scoredf <- multiScore(ranked, upSetColc = gsc)

# reverse again scores
scoredf$Scores <- scoredf$Scores * -1

save(scoredf, file = file.path(dirOut, "breast_reactome_genebody_sigdown_singscore.Rdata"))


# Compute pathway scores p <0.05 #----

## Islands promoters #----

load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"))
gene_meth <- readRDS(file.path(dirOut, "island_breast_genemeanmeth_promoters.Rds"))

mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

### Up #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut2, "breast_reactome_promoters_sigup_singscore.Rdata"))

### Down #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

# flip direction for ranking hypoM (negative M-values) first
mat <- mat * -1
ranked <- rankGenes(mat) 
scoredf <- multiScore(ranked, upSetColc = gsc)

# reverse again scores
scoredf$Scores <- scoredf$Scores * -1

save(scoredf, file = file.path(dirOut2, "breast_reactome_promoters_sigdown_singscore.Rdata"))

## Islands gene body #----

load(file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"))
gene_meth <- readRDS( file.path(dirOut, "island_breast_genemeanmeth_genebody.Rds"))

mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

### Up #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Up") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut2, "breast_reactome_genebody_sigup_singscore.Rdata"))

### Down #----
sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == "Down") %>% rownames()

gs_list <- gene_reactome_tbl %>%
  filter(gs_name %in% sig) %>%
  group_by(gs_name) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
gs_list <- setNames(gs_list$genes, gs_list$gs_name)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(gs_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = gs_list[[nm]]
    )
  })
)

# flip direction
mat <- mat * -1
ranked <- rankGenes(mat) 
scoredf <- multiScore(ranked, upSetColc = gsc)

# reverse again scores
scoredf$Scores <- scoredf$Scores * -1

save(scoredf, file = file.path(dirOut2, "breast_reactome_genebody_sigdown_singscore.Rdata"))