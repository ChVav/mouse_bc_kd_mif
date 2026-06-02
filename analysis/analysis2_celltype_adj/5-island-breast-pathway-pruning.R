
# For Reactome(/GO) separately, for each gene region and CAMERA directional results, do DAG-aware results pruning
# Recompute FDR on remaining pathways

library(here)
library(tidyverse)
library(patchwork)
library(igraph)
library(tidygraph)
library(ggraph)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/5-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_pathway_summary.R"))

# Reactome #----

# Build a pseudo-DAG from all genes on array we have data on (in beta matrix)
gene_sets <- readRDS(file = file.path(results_dir,"0-output/reactome_gene_universe.Rds")) 
g_pseudo <- build_pseudo_dag(gene_sets)
saveRDS(g_pseudo, file =file.path(dirOut,"pseudo_dag_reactome.Rds"))

## Promoter #----

load(file.path(results_dir,"analysis2_celltype_adj/4-output/out_reactome_promoter.Rdata"))
results_camera <- as.data.frame(out$camera)
gene_meta <- out$gene_meta

### Up #-----

sig <- results_camera %>% rownames_to_column(var = "Pathways") %>%
  dplyr::filter(Direction=="Up" & PValue < 0.05)
pruned_up <- prune_paths_from_dag(g_pseudo, sig$Pathways, approach = "general")
# 2 not in provided DAG.
# REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS2 not in provided DAG.
# REACTOME_CARDIOGENESIS

# Save visual check pseudo-subgraph
pdag_prom_up <- plot_dag(g_pseudo, paths_sub = sig$Pathways, layout = "kk")

# Update pathway summary
non_sig <- results_camera %>% dplyr::filter(Direction=="Up" & PValue >= 0.05) %>% rownames()
results_camera_pruned_up <- results_camera %>% dplyr::filter(rownames(.) %in% c(non_sig,pruned_up))

### Down #-----
sig <- results_camera %>% rownames_to_column(var = "Pathways") %>%
  dplyr::filter(Direction=="Down" & PValue < 0.05)
pruned_down <- prune_paths_from_dag(g_pseudo, sig$Pathways, approach = "general")

# Save visual check pseudo-subgraph
pdag_prom_down <- plot_dag(g_pseudo, paths_sub = sig$Pathways, layout = "kk")

# Update pathway summary
non_sig <- results_camera %>% dplyr::filter(Direction=="Down" & PValue >= 0.05) %>% rownames()
results_camera_pruned_down <- results_camera %>% dplyr::filter(rownames(.) %in% c(non_sig,pruned_down))

### Recombine, FDR #---------
results_camera_pruned <- rbind(results_camera_pruned_up,results_camera_pruned_down) %>%
  dplyr::mutate(FDR = stats::p.adjust(PValue, "BH"))

out <- list(gene_meta = gene_meta,
            camera = results_camera,
            camera_pruned = results_camera_pruned)
save(out, file = file.path(dirOut, "out_reactome_promoter_pruned.Rdata"))

## Gene body #----

load(file.path(results_dir,"analysis2_celltype_adj/4-output/out_reactome_genebody.Rdata"))
results_camera <- as.data.frame(out$camera)
gene_meta <- out$gene_meta

### Up #-----

sig <- results_camera %>% rownames_to_column(var = "Pathways") %>%
  dplyr::filter(Direction=="Up" & PValue < 0.05)
pruned_up <- prune_paths_from_dag(g_pseudo, sig$Pathways, approach = "general")
# 2 not in provided DAG.
# REACTOME_TRANSCRIPTIONAL_REGULATION_OF_TESTIS_DIFFERENTIATION2 not in provided DAG.
# REACTOME_CARDIOGENESIS

# Save visual check pseudo-subgraph
pdag_body_up <- plot_dag(g_pseudo, paths_sub = sig$Pathways, layout = "kk")

# Update pathway summary
non_sig <- results_camera %>% dplyr::filter(Direction=="Up" & PValue >= 0.05) %>% rownames()
results_camera_pruned_up <- results_camera %>% dplyr::filter(rownames(.) %in% c(non_sig,pruned_up))

### Down #-----
sig <- results_camera %>% rownames_to_column(var = "Pathways") %>%
  dplyr::filter(Direction=="Down" & PValue < 0.05)
pruned_down <- prune_paths_from_dag(g_pseudo, sig$Pathways, approach = "general")
# 1 not in provided DAG.
# REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION

# Save visual check pseudo-subgraph
pdag_body_down <- plot_dag(g_pseudo, paths_sub = sig$Pathways, layout = "kk")

# Update pathway summary
non_sig <- results_camera %>% dplyr::filter(Direction=="Down" & PValue >= 0.05) %>% rownames()
results_camera_pruned_down <- results_camera %>% dplyr::filter(rownames(.) %in% c(non_sig,pruned_down))

### Recombine, FDR #---------
results_camera_pruned <- rbind(results_camera_pruned_up,results_camera_pruned_down) %>%
  dplyr::mutate(FDR = stats::p.adjust(PValue, "BH"))

out <- list(gene_meta = gene_meta,
            camera = results_camera,
            camera_pruned = results_camera_pruned)
save(out, file = file.path(dirOut, "out_reactome_genebody_pruned.Rdata"))

## Plot DAG all significant pathways #----

p1 <- pdag_prom_up + ggtitle("pCGI, hyperM")
p2 <- pdag_prom_down + ggtitle("pCGI, hypoM")
p3 <- pdag_body_up + ggtitle("oCGI, hyperM")
p4 <- pdag_body_down + ggtitle("oCGI, hypoM")

des <- c("AABB
          AABB
          CCDD")

p <- p1 + p3 + p2 + p4 + plot_layout(design = des) + plot_annotation(tag_levels = "a")

ggsave(p,
       file = file.path(dirOut,"breast_islands_reactome_sig_dag.pdf"),
       width = 183*4,
       height = 200*4,
       units = "mm"
)

ggsave(p,
       file = file.path(dirOut,"FigS1_DAGpruning.png"),
       width = 183*4,
       height = 200*4,
       units = "mm"
)


