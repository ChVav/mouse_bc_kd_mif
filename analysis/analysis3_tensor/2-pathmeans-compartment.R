# Calculate signed pathway means in baseline cell, tca compartments for pathways discovered from celltype adjusted data
# TCA gives latent signals, these are interpreted per CpG and per cell type: higher is higher methylation, lower is lower
# These were then scaled for cross-cg and cross-sample comparison
# Note that I tried singscore, but doesnt give good results (scaling issues up/down?)

library(here)
library(tidyverse)
library(singscore)
library(GSEABase)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis3_tensor/2-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}
dirOut2 <- file.path(dirOut,"p0.05")
if(!dir.exists(dirOut2)){dir.create(dirOut2, recursive = TRUE)}

helper_pathway_means <- function(gene_meth_path, out_path, gene_reactome_tbl, dir_var, filepaths = c(dirOut,dirOut2), filename_out) {
  
  gene_meth <- readRDS(file = gene_meth_path) %>%
    column_to_rownames(var = "gene_primary") %>%
    dplyr::select(-n_cpg) %>%
    as.matrix()
  load(out_path)
  
  # p_thres = 0.01
  sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.01) %>% dplyr::filter(Direction == dir_var) %>% rownames()
  
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
  
  path_mean_scores <- lapply(names(gs_list), function(gs) {
    gs_genes <- intersect(gs_list[[gs]], rownames(gene_meth))
    scores <- colMeans(gene_meth[gs_genes, , drop = FALSE], na.rm = TRUE) 
    tibble::tibble(
      Pathway = gs,
      Sample  = names(scores),
      Score   = as.numeric(scores),
      n_genes = length(gs_genes)
    )
  }) |> dplyr::bind_rows()
  
  # Save
  save(path_mean_scores, file = file.path(filepaths[1],filename_out))
  
  # p_thres = 0.05
  sig <- out$camera_pruned %>% dplyr::filter(PValue < 0.05) %>% dplyr::filter(Direction == dir_var) %>% rownames()
  
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
  
  path_mean_scores <- lapply(names(gs_list), function(gs) {
    gs_genes <- intersect(gs_list[[gs]], rownames(gene_meth))
    scores <- colMeans(gene_meth[gs_genes, , drop = FALSE], na.rm = TRUE) 
    tibble::tibble(
      Pathway = gs,
      Sample  = names(scores),
      Score   = as.numeric(scores),
      n_genes = length(gs_genes)
    )
  }) |> dplyr::bind_rows()
  
  # Save
  save(path_mean_scores, file = file.path(filepaths[2],filename_out))
  
}

# Data load #----

gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds")) %>%
  dplyr::select(gene_primary, gs_name)

# Epithelial #----

## Islands promoters #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_epi_reactome_promoters_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_epi_reactome_promoters_sigdown_pathmean.Rdata")

## Gene body #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_epi_reactome_genebody_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_epi_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_epi_reactome_genebody_sigdown_pathmean.Rdata")

# Immune #----

## Islands promoters #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_imm_reactome_promoters_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_imm_reactome_promoters_sigdown_pathmean.Rdata")

## Gene body #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_imm_reactome_genebody_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_imm_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_imm_reactome_genebody_sigdown_pathmean.Rdata")

# Stromal #----

## Islands promoters #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_stro_reactome_promoters_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_promoters.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_promoter_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_stro_reactome_promoters_sigdown_pathmean.Rdata")

## Gene body #----

### Up #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Up",
                     filename_out = "base_stro_reactome_genebody_sigup_pathmean.Rdata")

### Down #-----------
helper_pathway_means(gene_meth_path = file.path(results_dir, "analysis3_tensor/1-output/base_stro_genemeanmeth_genebody.Rds"),
                     out_path = file.path(results_dir,"analysis2_celltype_adj/5-output/out_reactome_genebody_pruned.Rdata"),
                     gene_reactome_tbl = gene_reactome_tbl,
                     dir_var = "Down",
                     filename_out = "base_stro_reactome_genebody_sigdown_pathmean.Rdata")


