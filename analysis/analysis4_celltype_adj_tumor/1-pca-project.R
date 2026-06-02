
# Project celltype compositions in tumor tissue on PCA space used for normal breast tissue
# Only baseline

library(here)
library(tidyverse)
library(patchwork)
library(showtext) # so arrows and delta are rendered properly
showtext_auto()
library(reshape2)

source(here("src/clr_ilr_transform.R"))

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis4_celltype_adj_tumor/1-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}


# Data load #----

load(here('data/beta_final.Rdata'))
load(here("data/pheno.Rdata"))


# Mammary gland
pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast', 'Breast tumor') & pheno$Treatment == "MPA/DMBA+ KD- MIF-",])
beta_sub <- beta_final[rownames(beta_final) %in% cpg_island,
                       pheno_sub$basename,
                       drop = FALSE]
# Remove rows with any NA
beta_sub <- beta_sub[!apply(beta_sub, 1, anyNA), ] 

# Combine fat and fibroblasts
dat <- pheno_sub
dat$EpiFibFatIc_Stromal <- dat$EpiFibFatIc_Fat + dat$EpiFibFatIc_Fib
dat <- dat[, !colnames(dat) %in% c("EpiFibFatIc_Fat", "EpiFibFatIc_Fib")]

outcomes <- intersect(cpg_island,rownames(beta_sub))

# PCA #---
# Note project both breast tissue at risk as well as tumor DNA to PCA space
celltype_cols <- colnames(dat)[grepl("EpiFibFatIc",colnames(dat))]
covars <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Stromal")
covar2 <- setdiff(celltype_cols, covars[1:2])

clr_mat <- clr_transform(dat,covar2, c("basename")) %>%
  column_to_rownames(var = "basename")

load(file.path(results_dir,"analysis2_celltype_adj/2-output/pca.Rdata"))
PC_new <- predict(pca, newdata = clr_mat)
pc_scores <- PC_new %>% as.data.frame() %>% rownames_to_column(var = "basename")

# Adjusted M-values #----

M <- t(beta_sub[rownames(beta_sub) %in% outcomes, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
M <- M %>%
  tibble::column_to_rownames(var = "basename") %>%
  t() %>%
  as.data.frame()

PCs <- pc_scores %>%
  filter(basename %in% colnames(M)) %>%
  tibble::column_to_rownames(var = "basename") %>%
  dplyr::select(PC1,PC2,PC3,PC4)

identical(rownames(PCs),colnames(M))

# Scale and regress PCs out
Z <- scale(PCs, center = TRUE, scale = FALSE)
Y <- t(M)                    # samples x CpGs
qrZ <- qr(Z)
coefZ <- qr.coef(qrZ, Y)     # effects of PCs on each CpG  (nPC x CpGs)
Y_adj <- Y - Z %*% coefZ     # remove PC-related variation
M_adj <- t(Y_adj)            # CpGs x samples


# Save M, M_adj, PCs and Mouse.IDs

saveRDS(M, file = file.path(dirOut, "M_12pairs.Rds"))
saveRDS(M_adj, file = file.path(dirOut, "Madj_12pairs.Rds"))

pheno_sub <- pheno_sub %>%
  dplyr::select(basename, Mouse.ID, Tissue.type) %>%
  dplyr::left_join(pc_scores, by = "basename") 
saveRDS(pheno_sub, file = file.path(dirOut,"pheno_12pairs.Rds"))
