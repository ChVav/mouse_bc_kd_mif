
# CAMERA on the adjusted per‑gene statistics for a competitive, correlation‑aware, threshold‑free gene set analysis.
# Pathways altered in KD, MIF and KD + MIF groups compared to baseline (no modifier), as well as KD + MIF versus KD cancer-only

library(here)
library(tidyverse)
library(metafor)
library(limma)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/9-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_camera.R"))

# Data load #----
anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds"))
dmr_result <- readRDS(file.path(results_dir,"analysis2_celltype_adj/7-output/breast_limma_simple.Rds"))
dmr_result$outcome <- gsub("_logit","",dmr_result$outcome)
gene_meta <- readRDS(file.path(results_dir,"analysis2_celltype_adj/8-output/gene_meta_universe.Rds"))

# Hallmark
gene_H_tbl <- readRDS(file.path(results_dir,"0-output/hallmark_gene_universe.Rds")) #6961 genes
colnames(gene_H_tbl) <- tolower(colnames(gene_H_tbl))

# Oncogenic signatures
gene_C6_tbl <- readRDS(file.path(results_dir,"0-output/oncsign_gene_universe.Rds"))
colnames(gene_C6_tbl) <- tolower(colnames(gene_C6_tbl))

# Reactome
gene_reactome_tbl <- readRDS(file.path(results_dir,"0-output/reactome_gene_universe.Rds"))
colnames(gene_reactome_tbl) <- tolower(colnames(gene_reactome_tbl))

# Test MIF #----
dmr_result2 <- dmr_result %>% dplyr::filter(group == "KD- MIF+")
gene_meta2 <- gene_meta %>% dplyr::filter(group == "KD- MIF+")

## Promoter #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == TRUE)

### Hallmark #----
# prep <- prep_camera_probe_bias_adjusted(dmr_result2, anno, gene_H_tbl, "promoter", TRUE)
# out <- do.call(camera_probe_bias_adjusted, prep)
# 
# # Check CpG count adjustment
# gene_meta <- out$gene_meta
# gene_meta$p_gene <- 2 * pnorm(-abs(gene_meta$z_gene))
# cor(gene_meta$n_cpg, -log10(gene_meta$p_gene), method = "spearman") # 0.25 bias

# Compare to filtered gene_meta generated with IPW
out <- run_camera(gene_meta3, gene_H_tbl)
cor(out$gene_meta$n_cpg, -log10(out$gene_meta$p_gene), method = "spearman") # -0.04 better

# lets stick with simpler approach of IPW
save(out, file = file.path(dirOut, "mif_hallmark_promoter.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "mif_oncsign_promoter.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "mif_reactome_promoter.Rdata"))

## Gene_body #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == FALSE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "mif_hallmark_genebody.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "mif_oncsign_genebody.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "mif_reactome_genebody.Rdata"))

# KD #----
gene_meta2 <- gene_meta %>% dplyr::filter(group == "KD+ MIF-")

## Promoter #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == TRUE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kd_hallmark_promoter.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kd_oncsign_promoter.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kd_reactome_promoter.Rdata"))

## Gene_body #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == FALSE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kd_hallmark_genebody.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kd_oncsign_genebody.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kd_reactome_genebody.Rdata"))

# KD + MIF #----
gene_meta2 <- gene_meta %>% dplyr::filter(group == "KD+ MIF+")

## Promoter #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == TRUE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kdmif_hallmark_promoter.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kdmif_oncsign_promoter.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kdmif_reactome_promoter.Rdata"))

## Gene_body #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == FALSE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kdmif_hallmark_genebody.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kdmif_oncsign_genebody.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kdmif_reactome_genebody.Rdata"))

# KD+MIF vs MIF #----
dmr_result <- readRDS(file.path(results_dir,"analysis2_celltype_adj/7-output/breast_limma_simple2.Rds"))
dmr_result$outcome <- gsub("_logit","",dmr_result$outcome)
gene_meta <- readRDS(file.path(results_dir,"analysis2_celltype_adj/8-output/gene_meta_universe2.Rds"))

gene_meta2 <- gene_meta %>% dplyr::filter(group == "KD+ MIF+")

## Promoter #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == TRUE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_hallmark_promoter.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_oncsign_promoter.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_reactome_promoter.Rdata"))

## Gene_body #----

gene_meta3 <- gene_meta2 %>% dplyr::filter(promoter == FALSE)

### Hallmark #----
out <- run_camera(gene_meta3, gene_H_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_hallmark_genebody.Rdata"))

### Oncogenic signatures #----
out <- run_camera(gene_meta3, gene_C6_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_oncsign_genebody.Rdata"))

### Reactome #----
out <- run_camera(gene_meta3, gene_reactome_tbl)
save(out, file = file.path(dirOut, "kdmif_vs_mif_reactome_genebody.Rdata"))

