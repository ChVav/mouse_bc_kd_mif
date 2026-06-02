
# Compute inverse-Variance Weighted (IVW) effect for all genes and gene regions. 
# Both gene promoter and gene body region, CpG count correction

library(here)
library(tidyverse)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/8-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_exposure-modifier-gene.R")) 

anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds")) # gene_meth_ipw function, can stay as is

# all vs baseline #----
dmr_result <- readRDS(file.path(results_dir,"analysis2_celltype_adj/7-output/breast_limma_simple.Rds"))
dmr_result$outcome <- gsub("_logit","",dmr_result$outcome)
dmr_result <- dplyr::left_join(dmr_result, anno, by = c("outcome" = "Name"))

# # Promoter
gene_meta <- gene_meth_ipw(dmr_result, gene_col = "gene_primary", group_cols = c("gene_primary", "group", "promoter"), adjust_n_cpg = TRUE)

check <- gene_meta %>% dplyr::filter(group == "KD- MIF+" & promoter == TRUE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.07 ok

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF-" & promoter == TRUE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.1 ok

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF+" & promoter == TRUE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.16 bias

check <- gene_meta %>% dplyr::filter(group == "KD- MIF+" & promoter == FALSE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.04 ok

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF-" & promoter == FALSE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.15 still light bias

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF+" & promoter == FALSE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.04 ok

saveRDS(gene_meta, file = file.path(dirOut, "gene_meta_universe.Rds"))

# KD+MIF vs KD #----
dmr_result <- readRDS(file.path(results_dir,"analysis2_celltype_adj/7-output/breast_limma_simple2.Rds"))
dmr_result$outcome <- gsub("_logit","",dmr_result$outcome)
dmr_result <- dplyr::left_join(dmr_result, anno, by = c("outcome" = "Name"))

# # Promoter
gene_meta <- gene_meth_ipw(dmr_result, gene_col = "gene_primary", group_cols = c("gene_primary", "group", "promoter"), adjust_n_cpg = TRUE)

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF+" & promoter == TRUE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.12 bias

check <- gene_meta %>% dplyr::filter(group == "KD+ MIF+" & promoter == FALSE)
cor(check$n_cpg, -log10(check$p_gene), method = "spearman") # -0.18 ok

saveRDS(gene_meta, file = file.path(dirOut, "gene_meta_universe2.Rds"))

