
# Create sample level gene-summaries of scaled latent signals
# Only baselevel (no modifier)

library(here)
library(tidyverse)
library(GSEABase)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis3_tensor/1-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

helper_gene_summary <- function(Z_std, celltype, anno, promoter_var, filepath_out) {
  
  beta <- Z_std[celltype] %>% as.data.frame() 
  colnames(beta) <- gsub(paste0(celltype,"\\."),"",colnames(beta))
  beta <- beta %>% dplyr::select(all_of(dat$basename))
  
  keep <- anno %>% dplyr::filter(promoter == promoter_var) %>% pull(Name)
  
  # Create sample-level gene summary
  beta_sub <- beta[rownames(beta) %in% intersect(cpg_island,keep),
                   dat$basename,
                   drop = FALSE]
  # sum(is.na(beta_sub)) # none
  # sum(beta_sub == 0, na.rm = TRUE) # none
  cpg_annot <- anno %>%
    dplyr::filter(Name %in% rownames(beta_sub)) %>%
    dplyr::select(Name, gene_primary)
  
  gene_meth <- beta_sub %>%
    rownames_to_column(var = "Name") %>%
    left_join(cpg_annot) %>%
    dplyr::filter(!is.na(gene_primary)) %>%
    group_by(gene_primary) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE),
              n_cpg = n(),
              .groups = "drop")
  
  saveRDS(gene_meth, file = filepath_out)
  
}

# Data load #----
anno <- readRDS(file.path(results_dir,"0-output/anno2_genes_filtered.Rds"))
cpg_island <- anno %>% dplyr::filter(Relation_to_Island == "Island") %>% pull(Name) # 29,022 islands

load(here("data/pheno.Rdata"))

# Mammary gland
pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast'),])

# Reformat for 2x2x2 design, and diet groups for raw plots
dat <- as.data.frame(pheno_sub) %>%
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
      levels = c("KD- MIF-", "KD- MIF+", "KD+ MIF-", "KD+ MIF+")
    )
  ) %>%
  dplyr::filter(diet4 == "KD- MIF-")

# subset scaled latent signals

load(file.path(results_dir,'analysis3_tensor/0-output/tca.tensor.z.Rdata'))
cell_types <- names(Z_std)
Z_std <- lapply(Z_std, function(Zh) Zh[, dat$basename, drop = FALSE])

# Epithelial #----

## Islands promoters #----
helper_gene_summary(Z_std, "Epithelial", anno, TRUE, file.path(dirOut, "base_epi_genemeanmeth_promoters.Rds"))

## Islands gene body #----
helper_gene_summary(Z_std, "Epithelial", anno, FALSE, file.path(dirOut, "base_epi_genemeanmeth_genebody.Rds"))

# Immune #----

## Islands promoters #----
helper_gene_summary(Z_std, "Immune", anno, TRUE, file.path(dirOut, "base_imm_genemeanmeth_promoters.Rds"))

## Islands gene body #----
helper_gene_summary(Z_std, "Immune", anno, FALSE, file.path(dirOut, "base_imm_genemeanmeth_genebody.Rds"))

# Stromal #----

## Islands promoters #----
helper_gene_summary(Z_std, "Stromal", anno, TRUE, file.path(dirOut, "base_stro_genemeanmeth_promoters.Rds"))

## Islands gene body #----
helper_gene_summary(Z_std, "Stromal", anno, FALSE, file.path(dirOut, "base_stro_genemeanmeth_genebody.Rds"))

