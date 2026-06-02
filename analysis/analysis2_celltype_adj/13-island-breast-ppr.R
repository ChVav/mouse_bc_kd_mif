
# Note on PPRs
# when using the mass injected scores (paths below), nicely see increase of ppr with abs(mean delta) :)
# ranked_prom <- readRDS("results/analysis2_celltype_adj/3-output/ppr_ranked_prom.Rds")
# ranked_body <- readRDS("results/analysis2_celltype_adj/3-output/ppr_ranked_body.Rds")
# but the other ppr score can be reused for bayes posterior analysis

# Stick with extended seeds (see diagnostic plots)
# Save gene sets proximal to different seeds (supplementary plot)
# Singscore for hyperM and hypoM gene sets seperately

library(tidyverse)
library(grid)
library(gridExtra)
library(ggpubr)
library(here)
library(ggrepel)
library(showtext) # so arrows and delta are rendered properly
showtext_auto()
library(scales)
library(patchwork)
library(singscore)
library(GSEABase)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis2_celltype_adj/13-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_pageRank.R"))
source(here("src/helpers_pathway_summary.R"))

# Pheno #----
load(here("data/pheno.Rdata"))

pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast'),])

# Reformat for 2x2x2 design, and diet groups for raw plots
dat <- pheno_sub %>%
  mutate(
    exposure      = factor(ifelse(grepl("\\bMPA/DMBA\\+", Treatment), 1, 0),
                           levels = c(0, 1), labels = c("P/D-", "P/D+")),
    keto          = factor(ifelse(grepl("\\bKD\\+",       Treatment), 1, 0),
                           levels = c(0, 1), labels = c("KD-", "KD+")),
    mifepristone  = factor(ifelse(grepl("\\bMIF\\+",      Treatment), 1, 0),
                           levels = c(0, 1), labels = c("MIF-", "MIF+"))
  ) %>%
  mutate(
    diet4 = factor(
      paste0(keto, " ", mifepristone),
      levels = c("KD- MIF-",  "KD+ MIF-", "KD- MIF+", "KD+ MIF+")
    )
  ) %>%
  dplyr::select(basename,exposure,keto,mifepristone,diet4)

# DMP results baseline
dmp_base <- readRDS(here("results/analysis2_celltype_adj/2-output/breast_limma_celltype_adj.Rds")) %>%
  dplyr::filter(keto == "KD-" & mifepristone == "MIF+")
dmp_base$outcome <- gsub("_logit","", dmp_base$outcome)
anno <- readRDS(here("results/0-output/anno2_genes_filtered.Rds"))
dmp_base <- left_join(dmp_base,anno, by = c("outcome" = "Name"))

dmp_prom <- dmp_base %>% dplyr::filter(promoter == TRUE & p.adj_cell < 0.05) %>% pull(ENTREZID) %>% unique() %>% na.omit() #310 gene prim, 294 entrez
dmp_body <- dmp_base %>% dplyr::filter(promoter == FALSE & p.adj_cell < 0.05) %>% pull(ENTREZID) %>% unique() %>% na.omit() #498 gene prim, 468 entrez

# gene_primary, entrez
annot <- anno %>%
  dplyr::select(gene_primary, ENTREZID) %>%
  distinct() #27031

id2symbol <- setNames(
  annot$gene_primary,
  annot$ENTREZID
)

# Gene-level summary methy base
dat_base <- readRDS(here("results/analysis2_celltype_adj/3-output/gene_base_meta_universe.Rds")) %>%
  dplyr::left_join(annot, by = "gene_primary") %>%
  dplyr::filter(ENTREZID %in% c(dmp_prom, dmp_body)) # focus on genes with DMP

# gene summaries epithelial cell type compartments #-----

gene_met_epi_prom <- readRDS("results/analysis3_tensor/1-output/base_epi_genemeanmeth_promoters.Rds")
gene_met_epi_body <- readRDS("results/analysis3_tensor/1-output/base_epi_genemeanmeth_genebody.Rds")

mean_diff_genes <- function(expr_df, meta_df, group_var = "exposure",
                            ref = "P/D-", contrast = "P/D+") {
  
  dat <- expr_df %>%
    pivot_longer(
      cols = -gene_primary,
      names_to = "basename",
      values_to = "value"
    ) %>% 
    inner_join(meta_df, by = "basename")
  
  # Mean difference per gene
  res <- dat %>%
    filter(.data[[group_var]] %in% c(ref, contrast)) %>%
    group_by(gene_primary) %>%
    summarise(
      mean_ref = mean(value[.data[[group_var]] == ref], na.rm = TRUE),
      mean_contrast = mean(value[.data[[group_var]] == contrast], na.rm = TRUE),
      diff = mean_contrast - mean_ref,
      .groups = "drop"
    )
  
  return(res)
}

mean_met_epi_prom <- mean_diff_genes(gene_met_epi_prom, dat %>% filter(basename %in% colnames(gene_met_epi_prom)))
mean_met_epi_body <- mean_diff_genes(gene_met_epi_body, dat %>% filter(basename %in% colnames(gene_met_epi_body)))

configs <- list(
  list(file = "ppr_pr.Rds",              title = "PR core co-regulators"),
  list(file = "ppr_pr_ext.Rds",          title = "Extended PR regulatory network"),
  list(file = "ppr_rankl.Rds",           title = "RANKL"),
  list(file = "ppr_rankl_nfbk.Rds",      title = "RANKL–NF-κB"),
  list(file = "ppr_wnt.Rds",             title = "Wnt/β-catenin signaling"),
  list(file = "ppr_wnt_ext.Rds",         title = "Extended Wnt/β-catenin signaling"),
  list(file = "ppr_ccnd1.Rds",           title = "CCND1/cell-cycle"),
  list(file = "ppr_ccnd1_ext.Rds",       title = "Extended CCND1/cell-cycle"),
  list(file = "ppr_prl_stat5.Rds",       title = "PRL–JAK2–STAT5 signaling"),
  list(file = "ppr_pi3k_akt_mtor.Rds",   title = "PI3K–AKT–mTOR signaling"),
  list(file = "ppr_survival_nfkb.Rds",   title = "NF-κB survival")
)

base_path <- "results/2-output/"

# Pval, log2 fold enrichment #----

plot_list <- list()

for (cfg in configs) {
  ranked <- readRDS(file.path(base_path, cfg$file))
  
  plots <- list(
    plot_ppr_base_pval_logfold(
      ranked, dat_base,
      promoter_var = TRUE,
      seed_tit = cfg$title,
      top_genes = 10
    ) + ggtitle("Islands in promoter region"),
    
    plot_ppr_base_pval_logfold(
      ranked, dat_base,
      promoter_var = FALSE,
      seed_tit = cfg$title,
      top_genes = 10
    ) + ggtitle("Islands in gene body region")
  )
  
  plot_list <- c(plot_list, plots)
}

pdf(file.path(dirOut,"pval_logfold.pdf"), width = 16, height = 24)

wrap_plots(plot_list, ncol = 4) +
  plot_annotation(tag_levels = "a")

dev.off()

# Delta bulk + epi #-----

plot_list <- list()

for (cfg in configs) {
  ranked <- readRDS(file.path(base_path, cfg$file))
  
  plots <- list(
    plot_ppr_base_epi(
      ranked, dat_base, mean_met_epi_prom,
      promoter_var = TRUE,
      seed_tit = cfg$title,
      top_genes = 10
    ) + ggtitle("Islands in promoter region"),
    
    plot_ppr_base_epi(
      ranked, dat_base, mean_met_epi_body,
      promoter_var = FALSE,
      seed_tit = cfg$title,
      top_genes = 10
    ) + ggtitle("Islands in gene body region")
  )
  
  plot_list <- c(plot_list, plots)
}

pdf(file.path(dirOut,"pval_epi.pdf"), width = 16, height = 24)

wrap_plots(plot_list, ncol = 4) +
  plot_annotation(tag_levels = "a")

dev.off()

# Gather significantly proximal genes to seeds #----

configs <- list(
  list(file = "ppr_pr_ext.Rds",          title = "PR regulatory network"),
  list(file = "ppr_rankl_nfbk.Rds",      title = "RANKL–NF-κB"),
  list(file = "ppr_wnt_ext.Rds",         title = "Wnt/β-catenin signaling"),
  list(file = "ppr_ccnd1_ext.Rds",       title = "CCND1/cell-cycle"),
  list(file = "ppr_prl_stat5.Rds",       title = "PRL–JAK2–STAT5 signaling"),
  list(file = "ppr_pi3k_akt_mtor.Rds",   title = "PI3K–AKT–mTOR signaling"),
  list(file = "ppr_survival_nfkb.Rds",   title = "NF-κB survival")
)
base_path <- "results/2-output/"

df <- list()

for (cfg in configs) {
  ranked <- readRDS(file.path(base_path, cfg$file))
  
  dat <- ranked %>%
    dplyr::left_join(dat_base, by = "gene_primary") %>%
    na.omit() %>%
    dplyr::filter(p_emp < 0.05) %>%
    dplyr::mutate(seed = cfg$title) %>%
    dplyr::select(seed, ENTREZID, gene_primary, promoter)
  
  df[[length(df) + 1]] <- dat
  
}

df <- bind_rows(df)

saveRDS(df, file = file.path(dirOut, "ppr_proximal_genes.Rds"))

# Summarize gene (region) hubs using singscore #----

## Island promoters #----
# Reuse sample-level gene summary (made for pathway singscore summary)
gene_meth <- readRDS(file.path(results_dir,"analysis2_celltype_adj/6-output/island_breast_genemeanmeth_promoters.Rds"))
mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

hub_list <- df %>%
  filter(promoter == TRUE) %>%
  group_by(seed) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
hub_list <- setNames(hub_list$genes, hub_list$seed)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(hub_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = hub_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut, "breast_hubs_promoters_singscore.Rdata"))

## Island Gene body #----
# Reuse sample-level gene summary (made for pathway singscore summary)
gene_meth <- readRDS(file.path(results_dir,"analysis2_celltype_adj/6-output/island_breast_genemeanmeth_genebody.Rds"))
mat <- gene_meth %>% column_to_rownames(var = "gene_primary") %>% dplyr::select(-n_cpg)
ranked <- rankGenes(mat) 

hub_list <- df %>%
  filter(promoter == FALSE) %>%
  group_by(seed) %>%
  summarize(genes = list(sort(unique(gene_primary))), .groups = "drop")
hub_list <- setNames(hub_list$genes, hub_list$seed)

gsc <- GSEABase::GeneSetCollection(
  lapply(names(hub_list), function(nm) {
    GeneSet(
      setName = nm,
      geneIds = hub_list[[nm]]
    )
  })
)

scoredf <- multiScore(ranked, upSetColc = gsc)
save(scoredf, file = file.path(dirOut, "breast_hubs_genebody_singscore.Rdata"))



