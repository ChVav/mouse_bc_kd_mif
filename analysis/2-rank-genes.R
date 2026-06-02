# Rank all genes by Personalized PageRank proximity to different seeds
# damping: 1 - restart probability; default 0.85 is standard
# This is done without injecting any methylation evidence
# Cite
#Köhler et al., 2008, American Journal of Human Genetics (“Walking the Interactome”). Introduced random walk with restart on the interactome for gene prioritization from seed sets. Widely cited foundational paper for RWR/PPR on PPIs.
#Vanunu et al., 2010, Molecular Systems Biology. “Associating genes and protein complexes with disease via network propagation.” General network propagation framework closely related to PPR; shows how seed information diffuses over PPIs.
# Example application: https://link.springer.com/article/10.1186/1471-2105-6-233
# Progestin-driven tumorigenesis: pr, RANKL, wnt and ccnd1
# Modulation or sustaining: prl_stat5, pi3k_akt_mtor, survival_nfkb
# prl_stat5 <- c("Prlr","Jak2","Stat5a","Stat5b","Cish","Socs2")
# pi3k_akt_mtor <- c("Pik3ca","Pik3cb","Pik3r1","Akt1","Pten","Mtor","Rptor","Rheb")
# survival_nfkb <- c("Bcl2","Bcl2l1","Mcl1","Birc5")

library(here)
library(tidyverse)
library(patchwork)
library(igraph)

results_dir <- here("results")
dirOut <- file.path(results_dir,"2-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/epimods_functions.R"))
source(here("src/helpers_pageRank.R"))

# Load filtered ppis #----
biogrid <- readRDS(file.path(results_dir,"1-output/biogrid_mouse_ppi_entrez.Rds"))
intact <- readRDS(file.path(results_dir,"1-output/intact_mouse_ppi_entrez.Rds")) 

comb <- rbind(biogrid,intact) %>% #51126
  distinct() #49765

# Annotation file of genes I have beta-values on
annot <- readRDS("results/0-output/anno2_genes_filtered.Rds") %>%
  dplyr::select(gene_primary, ENTREZID) %>%
  distinct() #27031

id2symbol <- setNames(
  annot$gene_primary,
  annot$ENTREZID
)

# Genes in methylation data
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

# Core RANKL biology #----
rankl_core <- c("Tnfsf11","Tnfrsf11a","Tnfrsf11b","Traf6","Nfkb1","Rela")
seeds <- annot %>% dplyr::filter(gene_primary %in% rankl_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()

# Rank and save
ranked <- rank_genes_by_rankl_ppr(pin, array_genes = meth_genes, seeds = seeds$ENTREZID) # Tnfrsf11b (18383 not in PPI)
ranked <- left_join(ranked,annot, by = c("gene" = "ENTREZID"))

#saveRDS(ranked, file = file.path(dirOut, "ppr_ranked.Rds"))

# Run permutation test for pvalues
# Test whether each gene’s PPR score (with RANKL seeds) is higher than expected by chance, 
# after accounting for the fact that RANKL-seeds are high-degree hubs.
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)

# Attach symbols and save
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))

# Inspect top hits
head(res_perm, 20)
saveRDS(res_perm, file = file.path(dirOut, "ppr_rankl.Rds"))

# Extended RANKL core #----
rankl_nfbk_core <- c("Tnfsf11","Tnfrsf11a","Tnfrsf11b","Traf6","Nfkb1",
                     "Rela","Chuk","Ikbkb","Ikbkg","Nfkb2","Relb",
                     "Map3k14","Map3k7","Tab2","Traf2","Traf3",
                     "Nfkbia","Birc2","Birc3")

seeds <- annot %>% dplyr::filter(gene_primary %in% rankl_nfbk_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()

# Rank and save
ranked <- rank_genes_by_rankl_ppr(pin, array_genes = meth_genes, seeds = seeds$ENTREZID) # Tnfrsf11b (18383 not in PPI)
ranked <- left_join(ranked,annot, by = c("gene" = "ENTREZID"))

#saveRDS(ranked, file = file.path(dirOut, "ppr_ranked.Rds"))

# Run permutation test for pvalues
# Test whether each gene’s PPR score (with RANKL seeds) is higher than expected by chance, 
# after accounting for the fact that RANKL-seeds are high-degree hubs.
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)

# Attach symbols and save
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))

# Inspect top hits
head(res_perm, 20)
saveRDS(res_perm, file = file.path(dirOut, "ppr_rankl_nfbk.Rds"))

# Core pr biology #----
pr_core <-  c("Pgr","Ncoa1","Ncoa3","Ncor1","Ncor2","Pelp1","Kdm4b")
seeds <- annot %>% dplyr::filter(gene_primary %in% pr_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_pr.Rds"))

# Extended pr biology #----
pr_cor_ext <- c("Pgr","Ncoa1","Ncoa2","Ncoa3","Ncor1","Ncor2","Pelp1","Kdm4b",
                "Ep300","Crebbp","Med1","Kdm6b","Hsp90aa1","Hsp90ab1","Fkbp4","Fkbp5")
seeds <- annot %>% dplyr::filter(gene_primary %in% pr_cor_ext) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_pr_ext.Rds"))

# core WNT signaling #----
wnt_core <- c("Wnt4","Rspo1","Lgr5","Fzd2","Fzd7","Lrp5","Lrp6","Dvl2","Ctnnb1","Tcf7","Lef1","Axin2","Porcn")
seeds <- annot %>% dplyr::filter(gene_primary %in% wnt_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_wnt.Rds"))

#WNT extended #----
wnt_ext <- c("Wnt4","Rspo1","Lgr4","Lgr5","Fzd2","Fzd7","Lrp5","Lrp6","Dvl2","Ctnnb1","Tcf7","Tcf7l2","Lef1",
              "Axin2","Porcn","Gsk3b","Apc","Csnk1a1","Rnf43","Znrf3")
seeds <- annot %>% dplyr::filter(gene_primary %in% wnt_ext) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_wnt_ext.Rds"))

# CCND1 #----
ccnd1_core <- c("Ccnd1","Cdk4","Cdk6","Rb1","E2f1","E2f3")
seeds <- annot %>% dplyr::filter(gene_primary %in% ccnd1_core) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_ccnd1.Rds"))

# CCND1 extended #----
ccnd1_ext <- c("Ccnd1","Cdk4","Cdk6","Rb1","E2f1","E2f2","E2f3","Cdkn1b","Cdkn1a")
seeds <- annot %>% dplyr::filter(gene_primary %in% ccnd1_ext) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_ccnd1_ext.Rds"))

# Modulation or sustaining #----
prl_stat5 <- c("Prlr","Jak2","Stat5a","Stat5b","Cish","Socs2","Socs3","Elf5")
seeds <- annot %>% dplyr::filter(gene_primary %in% prl_stat5) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_prl_stat5.Rds"))

pi3k_akt_mtor <- c("Pik3ca","Pik3cb","Pik3r1","Pik3r2","Akt1","Pdpk1","Pten","Tsc1","Tsc2",
                   "Mtor","Rheb","Rptor","Rictor","Rps6kb1","Eif4ebp1")
seeds <- annot %>% dplyr::filter(gene_primary %in% pi3k_akt_mtor) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_pi3k_akt_mtor.Rds"))

survival_nfkb <- c("Bcl2","Bcl2l1","Mcl1","Birc5","Bcl2a1","Xiap","Traf1")
seeds <- annot %>% dplyr::filter(gene_primary %in% survival_nfkb) %>% dplyr::select(gene_primary, ENTREZID) %>% distinct()
set.seed(42)
res_perm <- empirical_pvals_degree_binned(
  g = pin,
  seeds = seeds$ENTREZID,
  array_genes = meth_genes,
  B = 50000,          
  damping = 0.85,
  nbins = 10,        
  exclude_true_seeds = TRUE,
  rng_seed = 42,
  verbose = TRUE
)
res_perm <- res_perm %>%
  left_join(annot, by = c("gene" = "ENTREZID")) %>%
  arrange(padj, p_emp, desc(ppr))
saveRDS(res_perm, file = file.path(dirOut, "ppr_survival_nfkb.Rds"))


