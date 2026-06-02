# Calculate hub means in baseline cell, tca compartments for pathways discovered from celltype adjusted data
# TCA gives latent signals, these are interpreted per CpG and per cell type: higher is higher methylation, lower is lower
# These were then scaled for cross-cg and cross-sample comparison

library(here)
library(tidyverse)
library(singscore)
library(GSEABase)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis3_tensor/3-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

helper_hub_means <- function(gene_meth_path, out_path, promoter_var, gene_reactome_tbl, file_out) {
  
  gene_meth <- readRDS(file = gene_meth_path) %>%
    column_to_rownames(var = "gene_primary") %>%
    dplyr::select(-n_cpg) %>%
    as.matrix()
  
  hub_list <- readRDS(out_path) %>%
    filter(promoter == promoter_var) %>%
    group_by(seed) %>%
    summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
  hub_list <- setNames(hub_list$genes, hub_list$seed)
  
  hub_mean_scores <- lapply(names(hub_list), function(seed) {
    hub_genes <- intersect(hub_list[[seed]], rownames(gene_meth))
    scores <- colMeans(gene_meth[hub_genes, , drop = FALSE], na.rm = TRUE) 
    tibble::tibble(
      Hub = seed,
      Sample  = names(scores),
      Score   = as.numeric(scores),
      n_genes = length(hub_genes)
    )
  }) |> dplyr::bind_rows()
  
  # Save
  save(hub_mean_scores, file = file_out)

}

# Data load #----

gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds")) %>%
  dplyr::select(gene_primary, gs_name)

# Epithelial #----

## Islands promoters #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                     promoter_var = TRUE,
                     gene_reactome_tbl = gene_reactome_tbl,
                     file_out = file.path(dirOut,"base_epi_hub_promoters_pathmean.Rdata"))

## Gene body #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_promoters.Rds"),
                 out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                 promoter_var = FALSE,
                 gene_reactome_tbl = gene_reactome_tbl,
                 file_out = file.path(dirOut,"base_epi_hub_genebody_pathmean.Rdata"))

# Immune #----

## Islands promoters #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_promoters.Rds"),
                 out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                 promoter_var = TRUE,
                 gene_reactome_tbl = gene_reactome_tbl,
                 file_out = file.path(dirOut,"base_imm_hub_promoters_pathmean.Rdata"))

## Gene body #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_promoters.Rds"),
                 out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                 promoter_var = FALSE,
                 gene_reactome_tbl = gene_reactome_tbl,
                 file_out = file.path(dirOut,"base_imm_hub_genebody_pathmean.Rdata"))

# Stromal #----

## Islands promoters #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_promoters.Rds"),
                 out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                 promoter_var = TRUE,
                 gene_reactome_tbl = gene_reactome_tbl,
                 file_out = file.path(dirOut,"base_stro_hub_promoters_pathmean.Rdata"))

## Gene body #----

helper_hub_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_promoters.Rds"),
                 out_path = file.path(results_dir,"analysis2_celltype_adj/13-output/ppr_proximal_genes.Rds"),
                 promoter_var = FALSE,
                 gene_reactome_tbl = gene_reactome_tbl,
                 file_out = file.path(dirOut,"base_stro_hub_genebody_pathmean.Rdata"))
