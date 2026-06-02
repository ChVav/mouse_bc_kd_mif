
# Check cell type changes
# Significant in healthy mice vs breast tissue vs tumor tissue?
# 4 groups: baseline, KD, MIF, KD + MIF

library(here)
library(tidyverse)
library(limma)
library(compositions)

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/1-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/helpers_exposure-modifier.R"))
source(here("src/clr_ilr_transform.R"))

# Check transformations with raw percentage scale #----
source(here("src/plot_box_raw.R"))

load("data/pheno.Rdata")
dat <- droplevels(pheno[pheno$Experiment=='exp5'& pheno$Tissue.type %in% c('Breast'),])

dat <- dat %>%
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
  mutate(
    diet5 = case_when(
      diet4 == "KD- MIF-" ~ "No modifier",
      diet4 == "KD- MIF+" ~ "MIF",
      diet4 == "KD+ MIF-" ~ "KD",
      diet4 == "KD+ MIF+" ~ "KD + MIF"
    )
  )

dat$diet5 <- factor(dat$diet5, levels = c("No modifier", "KD", "MIF", "KD + MIF"))

# Combine fat and fibroblasts, set outcomes
dat$EpiFibFatIc_Stromal <- dat$EpiFibFatIc_Fat + dat$EpiFibFatIc_Fib
dat <- dat[, !colnames(dat) %in% c("EpiFibFatIc_Fat", "EpiFibFatIc_Fib")]

outcome_cols <- colnames(dat)[grepl("EpiFibFatIc",colnames(dat))]
outcome_cols1 <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Stromal")
outcome_cols2 <- setdiff(outcome_cols, outcome_cols1[1:2])

## Broad  #-----------------

pdat <- dat %>% dplyr::filter(exposure == "P/D-") %>% dplyr::filter(Tissue.type == "Breast")

plot_list <- c()
count <- 1
# Raw values

for (outcome in outcome_cols1){
  plot_list[[count]] <- plot_box_raw(droplevels(dat[dat$Tissue.type %in% c('Breast'),]), outcome) + ggtitle("Raw")
  count <- count + 1
}

p <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect")

pdf(file.path(dirOut, "check_composition_broad.pdf"), height = 4, width = 16)
print(p)
dev.off()

## Hepi #-----------------

plot_list <- c()
count <- 1
# Raw values

for (outcome in outcome_cols2){
  plot_list[[count]] <- plot_box_raw(droplevels(dat[dat$Tissue.type %in% c('Breast'),]), outcome) + ggtitle("Raw")
  count <- count + 1
}

p <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect")

pdf(file.path(dirOut, "check_composition_hepi.pdf"), height = 12, width = 16)
print(p)
dev.off()


# Control mice (normal tissue), does the modifier affect cell type composition #-----
# Yes, but only global F-test

## Broad #----
pdat <- dat %>% dplyr::filter(exposure == "P/D-") %>% dplyr::filter(Tissue.type == "Breast")
clr_df <- clr_transform(pdat,outcome_cols1, c("diet5"))
design <- model.matrix(~ diet5, data = clr_df)
colnames(design) <- levels(clr_df$diet5)
clr_mat <- t(clr_df %>% dplyr::select(-diet5))
fit <- lmFit(clr_mat, design)
fit <- eBayes(fit)

topTable(fit, coef="KD")
topTable(fit, coef="MIF")
topTable(fit, coef="KD + MIF")

# Global F-test any changes?
topTable(fit, coef=NULL, number=Inf, sort.by="F")

## Hepi #----
pdat <- dat %>% dplyr::filter(exposure == "P/D-") %>% dplyr::filter(Tissue.type == "Breast")
clr_df <- clr_transform(pdat,outcome_cols2, c("diet5"))
design <- model.matrix(~ diet5, data = clr_df)
colnames(design) <- levels(clr_df$diet5)
clr_mat <- t(clr_df %>% dplyr::select(-diet5))
fit <- lmFit(clr_mat, design)
fit <- eBayes(fit)

topTable(fit, coef="KD")
topTable(fit, coef="MIF")
topTable(fit, coef="KD + MIF")

# Global F-test any changes?
topTable(fit, coef=NULL, number=Inf, sort.by="F")

# Cancer delta #-----

# See PDF Sfig_celltype_cancerdelta.Rmd
