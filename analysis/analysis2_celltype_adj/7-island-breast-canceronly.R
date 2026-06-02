# Run lm to compare diet/mife vs no mod; as well as KD+MIF vs MIF
# extract significant features

library(here)
library(tidyverse)
library(patchwork)
library(limma)
library(showtext) # so arrows and delta are rendered properly

showtext_auto()

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/7-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_modifier.R")) 

# Annotation #----
anno <- readRDS(file.path(results_dir,"0-output/anno_12.v1.mm10.Rds"))
cpg_island <- anno %>% dplyr::filter(Relation_to_Island == "Island") %>% rownames()
cpg_nonisland <- anno %>% dplyr::filter(Relation_to_Island != "Island") %>% rownames()
promotor_cg <- anno %>% dplyr::filter(promoter == TRUE) %>% pull(Name)
body_cg <- anno %>% dplyr::filter(primary_feature == "tss_body") %>% pull(Name)

# Data load #----

load('data/beta_final.Rdata')
load("data/pheno.Rdata")

# Mammary gland
pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast'),])
beta_sub <- beta_final[rownames(beta_final) %in% cpg_island,
                       pheno_sub$basename,
                       drop = FALSE]
# Remove rows with any NA
beta_sub <- beta_sub[!apply(beta_sub, 1, anyNA), ] 

# Reformat for 1x1x1 design, and diet groups for raw plots
dat <- as.data.frame(pheno_sub) %>%
  mutate(
    exposure      = factor(ifelse(grepl("\\bMPA/DMBA\\+", Treatment), 1, 0),
                           levels = c(0, 1), labels = c("P/D-", "P/D+")),
    keto          = factor(ifelse(grepl("\\bKD\\+",       Treatment), 1, 0),
                           levels = c(0, 1), labels = c("KD-", "KD+")),
    mifepristone  = factor(ifelse(grepl("\\bMIF\\+",      Treatment), 1, 0),
                           levels = c(0, 1), labels = c("MIF-", "MIF+"))
  ) %>%
  dplyr::filter(exposure == "P/D+") %>%
  mutate(
    diet4 = factor(
      paste0(keto, " ", mifepristone),
      levels = c("KD- MIF-", "KD- MIF+", "KD+ MIF-", "KD+ MIF+")
    )
  ) %>%
  droplevels() %>%
  dplyr::select(-exposure, -Treatment)

# Combine fat and fibroblasts
dat$EpiFibFatIc_Stromal <- dat$EpiFibFatIc_Fat + dat$EpiFibFatIc_Fib
dat <- dat[, !colnames(dat) %in% c("EpiFibFatIc_Fat", "EpiFibFatIc_Fib")]

outcomes <- intersect(cpg_island,rownames(beta_sub))

# Beta logit #----
beta_sub <- beta_sub[,dat$basename, drop = FALSE]
beta_selected <- t(beta_sub[rownames(beta_sub) %in% outcomes, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
outcomes <- paste0(outcomes,"_logit")
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))

# Add PCs for correction #----
pc_scores <- readRDS(file.path(results_dir,"analysis2_celltype_adj/2-output/celltype_pc_scores.Rds"))
dat_covar1 <- left_join(dat, pc_scores[,1:5]) %>%
  full_join(beta_selected)
covars <- colnames(pc_scores[,1:5])[-1]

# run Limma with KD- MIF- as baseline #----
res1 <- summarize_outcomes_limma_simple(
  dat = dat_covar1,
  outcome_cols = outcomes,
  moderate = T,
  covars = covars
)

res1 <- res1 %>%
  group_by(group) %>%
  mutate(p.adj_cell = p.adjust(p.value, "BH")) %>%
  ungroup()

res1$group <- gsub("diet4","",res1$group)

saveRDS(res1, file = file.path(dirOut, "breast_limma_simple.Rds"))

# Quick check

res1 %>% dplyr::filter(group == "KD+ MIF-") %>% dplyr::filter(p.adj_cell < 0.05) %>% nrow() #0
res1 %>% dplyr::filter(group == "KD- MIF+") %>% dplyr::filter(p.adj_cell < 0.05) %>% nrow() #147
res1 %>% dplyr::filter(group == "KD+ MIF+") %>% dplyr::filter(p.adj_cell < 0.05) %>% nrow() #67

res1 %>% dplyr::filter(group == "KD+ MIF-") %>% dplyr::filter(p.value < 0.05) %>% nrow() #2423
res1 %>% dplyr::filter(group == "KD- MIF+") %>% dplyr::filter(p.value < 0.05) %>% nrow() #5064
res1 %>% dplyr::filter(group == "KD+ MIF+") %>% dplyr::filter(p.value < 0.05) %>% nrow() #3474

# run Limma with KD+ MIF- as baseline #---
res1 <- summarize_outcomes_limma_simple2(
  dat = dat_covar1,
  outcome_cols = outcomes,
  moderate = T,
  covars = covars
)

res1 <- res1 %>%
  group_by(group) %>%
  mutate(p.adj_cell = p.adjust(p.value, "BH")) %>%
  ungroup()

res1$group <- gsub("diet4","",res1$group)

saveRDS(res1, file = file.path(dirOut, "breast_limma_simple2.Rds"))

# Quick check

res1 %>% dplyr::filter(group == "KD+ MIF+") %>% dplyr::filter(p.adj_cell < 0.05) %>% nrow() #0
res1 %>% dplyr::filter(group == "KD+ MIF+") %>% dplyr::filter(p.value < 0.05) %>% nrow() #1346


