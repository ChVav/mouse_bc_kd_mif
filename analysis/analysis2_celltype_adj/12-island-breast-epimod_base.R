# Detect cancer epimods: delta-exposure at baseline
# Use Biogrid, intact combined ppi
# note the newer FEM package does not seem to be maintained
# Have revised the original epimod implementation, to be compatible with my setup and newer iGraph version
# Source code revised from original epimod implementation
#https://pmc.ncbi.nlm.nih.gov/articles/PMC3620664/#sec10


library(igraph)
library(marray)
library(tidyverse)
library(here)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis2_celltype_adj/12-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/epimods_functions.R"))

# Load filtered ppis #----
biogrid <- readRDS(file.path(results_dir,"1-output/biogrid_mouse_ppi_entrez.Rds"))
intact <- readRDS(file.path(results_dir,"1-output/intact_mouse_ppi_entrez.Rds")) 

comb <- rbind(biogrid,intact) %>% #51126
 distinct() #49765

# Promoters #----
# Annotation file of genes I have beta-values on, keep only promoters
annot <- readRDS("results/0-output/anno2_genes_filtered.Rds") %>%
  dplyr::filter(promoter == TRUE) #56342 Cpgs, though more than one CpG per gene.

id2symbol <- setNames(
  annot$gene_primary,
  annot$ENTREZID
)

# Genes in methylation data
meth_genes <- unique(annot$ENTREZID)
meth_genes <- meth_genes[!is.na(meth_genes)]
meth_genes <- as.character(meth_genes)

# Stats
gene_meth <- readRDS(file.path(results_dir,"analysis2_celltype_adj/3-output/gene_base_meta_universe.Rds")) %>%
  dplyr::filter(promoter == TRUE)
annot <- annot %>% dplyr::select(gene_primary,ENTREZID) %>% na.omit() %>% distinct()
gene_meth <- left_join(gene_meth,annot)

# p-values
pval.v <- gene_meth$p_gene
names(pval.v) <- gene_meth$ENTREZID

# statistics (signed), scale?
stat.v <- gene_meth$z
#stat.v <- scale(stat.v)
names(stat.v) <- gene_meth$ENTREZID

# Make pin
edges_filtered <- comb[
  comb$A %in% meth_genes &
    comb$B %in% meth_genes,
]

pin <- graph_from_data_frame(edges_filtered, directed = FALSE)
pin <- simplify(pin)

# Sanity check
vcount(pin)   # 8311 genes
ecount(pin)   # 33826 interactions

# Detect and validate epimods
set.seed(123)
out <- run_epimod(pin,stat.v,pval.v)
save(out, file = file.path(dirOut, "epimods_promoter.Rdata"))

# Check some modules
sig.mods <- out[[3]] 
weight.g <- out[[5]]
pval.v <- out[[6]]
stat.v <- out[[7]]

check_genes <- annot %>% dplyr::filter(ENTREZID %in% sig.mods[[1]]) %>% pull(gene_primary)

# Node colours give significance (P < 0.05) of hypermethylation and hypomethylation. 
# Edge (aka link) colours range from white to black, indicating increasing weight with black denoting the top 2% of edges ranked by weight.
renderModule_gg(sig.mods[[3]], weight.g, pval.v, stat.v,
             vertex_cols = c("#4575b4","#d73027",mid="white"),
             k = 50,
             id2symbol = id2symbol,
             use_symbols = TRUE)
renderModule_gg(sig.mods[[4]], weight.g, pval.v, stat.v,
                vertex_cols = c("#4575b4","#d73027",mid="white"),
                k = 50,
                id2symbol = id2symbol,
                use_symbols = TRUE)

# jaccard <- function(a, b) {
#   inter <- length(intersect(a, b))
#   union <- length(union(a, b))
#   if (union == 0) return(0)
#   inter / union
# }
# 
# overlap_coef <- function(a, b) {
#   inter <- length(intersect(a, b))
#   denom <- min(length(a), length(b))
#   if (denom == 0) return(0)
#   inter / denom
# }
# 
# # Hypergeometric p-value for the overlap
# # N = universe size (e.g., number of genes in weight.g)
# overlap_p <- function(a, b, N) {
#   k <- length(intersect(a, b))
#   m <- length(a)
#   n <- N - m
#   # P(X >= k), X ~ Hypergeom(N, m, |b|)
#   phyper(k - 1, m, N - m, length(b), lower.tail = FALSE)
# }
# 
# # Example usage
# A <- sig.mods[[3]]
# B <- sig.mods[[4]]
# N <- length(V(weight.g))  # or length(unique(unlist(sig.mods)))
# cat("Jaccard:", jaccard(A, B), "\n")
# cat("Overlap coef:", overlap_coef(A, B), "\n")
# cat("Overlap p:", overlap_p(A, B, N), "\n")
# lower overlap coeff to 0.8


# Gene body #----
# Annotation file of genes I have beta-values on, keep only promoters
annot <- readRDS("results/0-output/anno2_genes_filtered.Rds") %>%
  dplyr::filter(promoter == FALSE) #159732 Cpgs, though more than one CpG per gene.

id2symbol <- setNames(
  annot$gene_primary,
  annot$ENTREZID
)

# Genes in methylation data
meth_genes <- unique(annot$ENTREZID)
meth_genes <- meth_genes[!is.na(meth_genes)]
meth_genes <- as.character(meth_genes)

# Stats
gene_meth <- readRDS(file.path(results_dir,"analysis2_celltype_adj/3-output/gene_base_meta_universe.Rds")) %>%
  dplyr::filter(promoter == FALSE)
annot <- annot %>% dplyr::select(gene_primary,ENTREZID) %>% na.omit() %>% distinct()
gene_meth <- left_join(gene_meth,annot)

# p-values
pval.v <- gene_meth$p_gene
names(pval.v) <- gene_meth$ENTREZID

# statistics (signed), scale?
stat.v <- gene_meth$z
#stat.v <- scale(stat.v)
names(stat.v) <- gene_meth$ENTREZID

# Make pin
edges_filtered <- comb[
  comb$A %in% meth_genes &
    comb$B %in% meth_genes,
]

pin <- graph_from_data_frame(edges_filtered, directed = FALSE)
pin <- simplify(pin)

# Sanity check
vcount(pin)   # 10762 genes
ecount(pin)   # 43244 interactions

# Detect and validate epimods
set.seed(123)
out <- run_epimod(pin,stat.v,pval.v)
save(out, file = file.path(dirOut, "epimods_genebody.Rdata"))

# Check some modules
sig.mods <- out[[3]] 
weight.g <- out[[5]]
pval.v <- out[[6]]
stat.v <- out[[7]]

renderModule_gg(sig.mods[[1]], weight.g, pval.v, stat.v,
                vertex_cols = c("#4575b4","#d73027",mid="white"),
                k = 50,
                id2symbol = id2symbol,
                use_symbols = TRUE)

# Node colours give significance (P < 0.05) of hypermethylation and hypomethylation. 
# Edge (aka link) colours range from white to black, indicating increasing weight with black denoting the top 2% of edges ranked by weight.

check_genes <- annot %>% dplyr::filter(ENTREZID %in% sig.mods[[1]]) %>% pull(gene_primary)
