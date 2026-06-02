
# Baseline stats only for island CpGs
# Compute inverse-Variance Weighted (IVW) effect for all genes and gene regions. CpG count correction
# Compute also personalized page rank (ppr) score with respect to RANKL seeds

library(here)
library(tidyverse)
library(igraph)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/3-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_exposure-modifier-gene.R"))
source(here("src/epimods_functions.R"))
source(here("src/helpers_pageRank.R"))

anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds"))

# Inverse-Variance Weighted Gene-Level Meta-Analysis #----

## Baseline #-------
dmr_result <- readRDS(file.path(results_dir,"analysis2_celltype_adj/2-output/breast_limma_celltype_adj.Rds")) %>% #29,022 Island CpGs with annotated genes
  dplyr::filter(keto == "KD-" & mifepristone == "MIF-") 
dmr_result$outcome <- gsub("_logit","",dmr_result$outcome)
dmr_result <- dplyr::left_join(dmr_result, anno, by = c("outcome" = "Name")) %>%
  dplyr::select(outcome, estimate, SE, p.value, p.adj_cell, gene_primary, promoter)

gene_meta <- gene_meth_ipw(dmr_result, gene_col = "gene_primary", group_cols = c("gene_primary", "promoter"))

# Check for bias
cor(gene_meta$n_cpg, -log10(gene_meta$p_gene), method = "spearman") # 0.23  CpG bias
prom <- gene_meta %>% dplyr::filter(promoter == TRUE)
body <- gene_meta %>% dplyr::filter(promoter != TRUE)
plot(prom$n_cpg,-log10(prom$p_gene)) # all good
cor(prom$n_cpg, -log10(prom$p_gene), method = "spearman") # 0.28 CpG bias
plot(body$n_cpg,-log10(body$p_gene)) # 2 genes with > 30 cpGs with very small p-value
cor(body$n_cpg, -log10(body$p_gene), method = "spearman") # 0.2

gene_meta <- gene_meth_ipw(dmr_result, gene_col = "gene_primary", group_cols = c("gene_primary", "promoter"), adjust_n_cpg = TRUE)
cor(gene_meta$n_cpg, -log10(gene_meta$p_gene), method = "spearman") # bias gone
prom <- gene_meta %>% dplyr::filter(promoter == TRUE)
body <- gene_meta %>% dplyr::filter(promoter != TRUE)
plot(prom$n_cpg,-log10(prom$p_gene)) # all good
cor(prom$n_cpg, -log10(prom$p_gene), method = "spearman") 
plot(body$n_cpg,-log10(body$p_gene)) 
cor(body$n_cpg, -log10(body$p_gene), method = "spearman") 

saveRDS(gene_meta, file = file.path(dirOut, "gene_base_meta_universe.Rds"))

# PPR #-----

## PPI, methylation and annotation #----
biogrid <- readRDS(file.path(results_dir,"1-output/biogrid_mouse_ppi_entrez.Rds"))
intact <- readRDS(file.path(results_dir,"1-output/intact_mouse_ppi_entrez.Rds")) 

comb <- rbind(biogrid,intact) %>% #51126
  distinct() #49765

# Annotation file of genes
annot <- anno %>%
  dplyr::select(gene_primary, ENTREZID) %>%
  distinct() #27031

id2symbol <- setNames(
  annot$gene_primary,
  annot$ENTREZID
)

meth_genes <- unique(annot$ENTREZID)
meth_genes <- meth_genes[!is.na(meth_genes)]
meth_genes <- as.character(meth_genes)

# Make pin
edges_filtered <- comb[
  comb$A %in% meth_genes &
    comb$B %in% meth_genes,
]

pin <- graph_from_data_frame(edges_filtered, directed = FALSE)
pin <- simplify(pin)

gene_meta_annot <- left_join(gene_meta, annot, by = "gene_primary")

# Core RANKL biology
rankl_core <- c("Tnfsf11","Tnfrsf11a","Tnfrsf11b","Traf6","Nfkb1","Rela")
seeds <- annot %>% dplyr::filter(gene_primary %in% rankl_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()

## PPR #----

# Promoter
z_meth <- gene_meta_annot %>% dplyr::filter(promoter == TRUE) %>% na.omit()
z_meth <- setNames(z_meth$z, z_meth$ENTREZID)

ranked_prom <- rank_genes_by_rankl_ppr_methyl(pin, z_vec = z_meth,
                                              seeds = seeds$ENTREZID,
                                              mode = "abs", lambda = 0.2) #18383 not in ppi
ranked_prom <- left_join(ranked_prom,annot, by = c("gene" = "ENTREZID"))
saveRDS(ranked_prom, file = file.path(dirOut, "ppr_ranked_prom.Rds"))

# Gene body
z_meth <- gene_meta_annot %>% dplyr::filter(promoter == FALSE) %>% na.omit()
z_meth <- setNames(z_meth$z, z_meth$ENTREZID)
ranked_body  <- rank_genes_by_rankl_ppr_methyl(pin, z_vec = z_meth,
                                               seeds = seeds$ENTREZID,
                                               mode = "abs", lambda = 0.2)
ranked_body <- left_join(ranked_body,annot, by = c("gene" = "ENTREZID"))
saveRDS(ranked_body , file = file.path(dirOut, "ppr_ranked_body.Rds"))

