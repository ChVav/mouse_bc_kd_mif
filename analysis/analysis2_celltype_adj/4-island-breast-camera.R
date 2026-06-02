
# CAMERA on the adjusted per‑gene statistics for a competitive, correlation‑aware, threshold‑free gene set analysis.
# Check only pathways altered with cancer in the mouse model (baseline, no modifier)
# Use gene summaries already made in previous script (IPW with CpG bias correction)

library(here)
library(tidyverse)
library(metafor)
library(limma)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/4-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_camera.R"))

# Data load #----
anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds"))
gene_meth <- readRDS(file.path(results_dir,"analysis2_celltype_adj/3-output/gene_base_meta_universe.Rds"))

# Hallmark
gene_H_tbl <- readRDS(file.path(results_dir,"0-output/hallmark_gene_universe.Rds")) #6961 genes
colnames(gene_H_tbl) <- tolower(colnames(gene_H_tbl))

# Oncogenic signatures
gene_C6_tbl <- readRDS(file.path(results_dir,"0-output/oncsign_gene_universe.Rds"))
colnames(gene_C6_tbl) <- tolower(colnames(gene_C6_tbl))

# Reactome
gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds"))
colnames(gene_reactome_tbl) <- tolower(colnames(gene_reactome_tbl))

# Test #----
## Hallmark promoter #----
gene_meta <- gene_meth %>% dplyr::filter(promoter == TRUE)
out <- run_camera(gene_meta, gene_H_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") 
save(out, file = file.path(dirOut, "out_hallmark_promoter.Rdata"))

## Hallmark gene_body #----
gene_meta <- gene_meth %>% dplyr::filter(promoter != TRUE)
out <- run_camera(gene_meta, gene_H_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") 

save(out, file = file.path(dirOut, "out_hallmark_genebody.Rdata"))

## Onc sign promoter #----
gene_meta <- gene_meth %>% dplyr::filter(promoter == TRUE)
out <- run_camera(gene_meta, gene_C6_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") 
save(out, file = file.path(dirOut, "out_oncsign_promoter.Rdata"))

## Onc sign gene_body #----
gene_meta <- gene_meth %>% dplyr::filter(promoter != TRUE)
out <- run_camera(gene_meta, gene_C6_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") 
save(out, file = file.path(dirOut, "out_oncsign_genebody.Rdata"))

## Reactome promoter #----
gene_meta <- gene_meth %>% dplyr::filter(promoter == TRUE)
out <- run_camera(gene_meta, gene_reactome_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman")
save(out, file = file.path(dirOut, "out_reactome_promoter.Rdata"))

## Reactome gene_body #----
gene_meta <- gene_meth %>% dplyr::filter(promoter != TRUE)
out <- run_camera(gene_meta, gene_reactome_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") 
save(out, file = file.path(dirOut, "out_reactome_genebody.Rdata"))
