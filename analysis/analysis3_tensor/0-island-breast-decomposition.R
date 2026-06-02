
# Decompose methylation-signals samples by cell-type
# Islands only
# Note tested this on subpopulations, but can't run all due to collinearity issues, broad cell types gave least inflation in baseline mice
# Fat, may be fat cells and smooth muscle cells (contamination due to GSE184410 fat tissue)

library(here)
library(tidyverse)
library(TCA)

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis3_tensor/0-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

source(here("src/vignette_analysis_tca.R"))

# Annotation #----
anno <- readRDS(file.path(results_dir,"0-output/anno_12.v1.mm10.Rds"))
cpg_island <- anno %>% dplyr::filter(Relation_to_Island == "Island") %>% rownames()

# Data load #----

load(here('data/beta_final.Rdata'))
load(here("data/pheno.Rdata"))

# Issue with PC doubled columns
tab   <- table(names(pheno))
keep  <- names(tab)[tab == 1]        # names that occur once
pheno  <- pheno[, keep, drop = FALSE]

# Mammary gland
pheno_sub <- droplevels(pheno[pheno$Experiment=='exp5' & pheno$Tissue.type %in% c('Breast'),])
beta_sub <- beta_final[rownames(beta_final) %in% cpg_island, 
                       pheno_sub$basename,
                       drop = FALSE]

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
  )

# Remove rows with any NA
beta_sub <- beta_sub[!apply(beta_sub, 1, anyNA), ] # 30298 probes left

rm(pheno, pheno_sub,beta_final); gc()

# Test only baseline/exposure, what types to include #-----

dat_base <- dat %>% dplyr::filter(diet4 == "KD- MIF-") %>% droplevels()
beta_base <- beta_sub[,dat_base$basename, drop = FALSE]

covars <- colnames(dat_base)[grepl("EpiFibFatIc",colnames(dat_base))]
covar1 <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Fib","EpiFibFatIc_Fat")
covar2 <- setdiff(covars, covar1[1:2])

## Most fine-grained possible #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
identical(colnames(beta_sub),rownames(cell_types))
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
cor(cell_types)
# CD4T/CD8T 0.975
# B/CD4T 0.995
# MDSC/Nk 0.986

# Combine CD4T, CD8T, B cells into lymphocytes
cell_types$Lymphocyte <- cell_types$CD4T + cell_types$CD8T + cell_types$B
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B")]

# Combine MDSC and NK
cell_types$Innate <- cell_types$MDSC + cell_types$NK
cell_types <- cell_types[, !colnames(cell_types) %in% c("MDSC", "NK")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial.basal", "Fat.and.muscle", "Fibroblast", "Luminal.progenitor", "Mature.luminal", "Neutrophile", "Lymphocyte", "Innate.immune")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                           W = cell_types,
                           C1 = type,
                           debug = T)
# Does not run

## Epithelial split, Immune bigger compartements, stromal #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
identical(colnames(beta_sub),rownames(cell_types))
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
# Combine CD4T, CD8T, B cells into lymphocytes
cell_types$Lymphocyte <- cell_types$CD4T + cell_types$CD8T + cell_types$B
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B")]

# Combine MDSC and NK
cell_types$Innate <- cell_types$MDSC + cell_types$NK
cell_types <- cell_types[, !colnames(cell_types) %in% c("MDSC", "NK")]

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial.basal", "Luminal.progenitor", "Mature.luminal", "Neutrophile", "Lymphocyte", "Innate.immune", "Stromal")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)
# Still does not run

## Epithelial split only luminal basal, Immune bigger compartements, stromal #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
identical(colnames(beta_sub),rownames(cell_types))
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
# Combine CD4T, CD8T, B cells into lymphocytes
cell_types$Lymphocyte <- cell_types$CD4T + cell_types$CD8T + cell_types$B
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B")]

# Combine MDSC and NK
cell_types$Innate <- cell_types$MDSC + cell_types$NK
cell_types <- cell_types[, !colnames(cell_types) %in% c("MDSC", "NK")]

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

# Combine luminal
cell_types$Epithelial.luminal <- cell_types$luminal_progenitor + cell_types$mature_luminal
cell_types <- cell_types[, !colnames(cell_types) %in% c("luminal_progenitor", "mature_luminal")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial.basal", "Neutrophile", "Lymphocyte", "Innate.immune", "Stromal", "Epithelial.luminal")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)
# Still does not run

## Epithelial split only luminal basal, Immune bigger compartments, fibroblasts and fat separate #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
identical(colnames(beta_sub),rownames(cell_types))
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
# Combine CD4T, CD8T, B cells into lymphocytes
cell_types$Lymphocyte <- cell_types$CD4T + cell_types$CD8T + cell_types$B
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B")]

# Combine MDSC and NK
cell_types$Innate <- cell_types$MDSC + cell_types$NK
cell_types <- cell_types[, !colnames(cell_types) %in% c("MDSC", "NK")]


# Combine luminal
cell_types$Epithelial.luminal <- cell_types$luminal_progenitor + cell_types$mature_luminal
cell_types <- cell_types[, !colnames(cell_types) %in% c("luminal_progenitor", "mature_luminal")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial.basal", "Fat.and.muscle", "Fibroblast", "Neutrophile", "Lymphocyte", "Innate.immune", "Epithelial.luminal")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)
# Still does not run

## Epithelial combined, Immune bigger compartements, stromal #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
# Combine CD4T, CD8T, B cells into lymphocytes
cell_types$Lymphocyte <- cell_types$CD4T + cell_types$CD8T + cell_types$B
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B")]

# Combine MDSC and NK
cell_types$Innate <- cell_types$MDSC + cell_types$NK
cell_types <- cell_types[, !colnames(cell_types) %in% c("MDSC", "NK")]

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

# Combine Epithelial
cell_types$Epithelial <- cell_types$luminal_progenitor + cell_types$mature_luminal + cell_types$basal
cell_types <- cell_types[, !colnames(cell_types) %in% c("luminal_progenitor", "mature_luminal", "basal")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Neutrophile", "Lymphocyte", "Innate.immune", "Stromal", "Epithelial")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)

# Does not run

## Broad panel #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar1)) %>% 
  column_to_rownames(var = "basename") 
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial", "Immune", "Fibroblast", "Fat")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)

# Joint testing p-values
joint <- tca.joint.baseline$gammas_hat_pvals.joint[,"exposure"]

# Marginally testing p-values
marg <- tca.joint.baseline$gammas_hat_pvals

#qq-plots
pdf(file.path(dirOut,"qq_test_celltypes_baseline1.pdf"), width = 15, height = 10)
plot_qq(list(joint, marg[,1], marg[,2], marg[,3], marg[,4]),
        labels = c("Joint test", colnames(marg)),
        ggarrange.nrow = 2,
        ggarrange.ncol = 3,
        experiment_wide_line = FALSE)
dev.off()

## Epi, immune, stromal #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar1)) %>% 
  column_to_rownames(var = "basename") 
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial", "Immune", "Stromal")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)

# Joint testing p-values
joint <- tca.joint.baseline$gammas_hat_pvals.joint[,"exposure"]

# Marginally testing p-values
marg <- tca.joint.baseline$gammas_hat_pvals

#qq-plots
pdf(file.path(dirOut,"qq_test_celltypes_baseline2.pdf"), width = 10, height = 10)
plot_qq(list(joint, marg[,1], marg[,2], marg[,3]),
        labels = c("Joint test", colnames(marg)),
        ggarrange.nrow = 2,
        ggarrange.ncol = 2,
        experiment_wide_line = FALSE)
dev.off()

## Epithelial subsets, Immune, stromal #----
cell_types <- dat_base %>% 
  dplyr::select(basename,all_of(covar2)) %>% 
  column_to_rownames(var = "basename") %>% 
  dplyr::select(-EpiFibFatIc_Mono) # Remove monocytes
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Fix collinearity in cell types
# Combine immune cells
cell_types$Immune <- cell_types$CD4T + cell_types$CD8T + cell_types$B + cell_types$MDSC + cell_types$NK + cell_types$Neutro
cell_types <- cell_types[, !colnames(cell_types) %in% c("CD4T", "CD8T", "B", "MDSC", "NK", "Neutro")]

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial.basal", "Luminal.progenitor", "Mature.luminal", "Immune", "Stromal")

# Prepare covariates
type <- dat_base %>% dplyr::select(basename,exposure) %>% column_to_rownames(var = "basename")
identical(colnames(beta_base),rownames(type))
identical(dat_base$basename, rownames(type))
type <- as.matrix(type)
type <- ifelse(type == "P/D-", 0, 1)

# Test
tca.joint.baseline <- tca(X = beta_base,
                          W = cell_types,
                          C1 = type,
                          debug = T)

# Joint testing p-values
joint <- tca.joint.baseline$gammas_hat_pvals.joint[,"exposure"]

# Marginally testing p-values
marg <- tca.joint.baseline$gammas_hat_pvals

#qq-plots
pdf(file.path(dirOut,"qq_test_celltypes_baseline3.pdf"), width = 15, height = 10)
plot_qq(list(joint, marg[,1], marg[,2], marg[,3], marg[,4], marg[,5]),
        labels = c("Joint test", colnames(marg)),
        ggarrange.nrow = 2,
        ggarrange.ncol = 3,
        experiment_wide_line = FALSE)
dev.off()

# Decompose all samples Epi, immune, stromal #----

# Prepare cell types
covar1 <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Fib","EpiFibFatIc_Fat")
cell_types <- dat %>% dplyr::select(basename,all_of(covar1))  
rownames(cell_types) <- cell_types$basename 
cell_types <- cell_types %>% dplyr::select(-basename)
colnames(cell_types) <- gsub("EpiFibFatIc_","",colnames(cell_types))

# Combine fat and fibroblasts
cell_types$Stromal <- cell_types$Fat + cell_types$Fib
cell_types <- cell_types[, !colnames(cell_types) %in% c("Fat", "Fib")]

cor(cell_types)

cell_types <- as.matrix(cell_types)
colnames(cell_types) <- c("Epithelial", "Immune", "Stromal")

identical(colnames(beta_sub),rownames(cell_types))

cell_types <- as.matrix(cell_types)

tca.joint <- tca(X = beta_sub, W = cell_types)
tca.tensor <- tensor(beta_sub, tca.joint, scale = FALSE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE )
names(tca.tensor) <- colnames(cell_types)
save(tca.joint, file = file.path(dirOut,"tca.joint.Rdata"))
save(tca.tensor, file = file.path(dirOut,"tca.tensor.Rdata"))

# Standardize the tca latent signals for cross-cpg and cross-cell type analyses #----
# Compute and freeze per‑CpG bulk anchors μ and σ using all available samples once.
# Apply those same μ/σ to standardize each cell‑type latent matrix.

## Prep aligned matrices #-----------

# Decomposed latent signals
load(file.path(results_dir,'analysis3_tensor/0-output/tca.tensor.Rdata'))
cell_types <- names(tca.tensor)

# Extract, clean column names, and align to selected samples
Z_raw <- setNames(vector("list", length(cell_types)), cell_types)
for (ct in cell_types) {
  Z <- tca.tensor[[ct]]                         
  colnames(Z) <- sub(paste0("^", ct, "\\."), "", colnames(Z))  
  Z_raw[[ct]] <- Z[, dat$basename, drop = FALSE]               
}

# Keep only CpGs and samples present in all cell types
common_cpgs    <- Reduce(intersect, lapply(Z_raw, rownames))
common_samples <- Reduce(intersect, lapply(Z_raw, colnames))
Z_raw <- lapply(Z_raw, function(Z) Z[common_cpgs, common_samples, drop = FALSE])

# Align CpGs/samples across all objects bulk beta_sub and decomposed
common_cpgs    <- Reduce(intersect, c(list(rownames(beta_sub)), lapply(Z_raw, rownames)))
common_samples <- Reduce(intersect, c(list(colnames(beta_sub)), lapply(Z_raw, colnames)))
Xb   <- beta_sub[common_cpgs, common_samples, drop = FALSE]
Z_al <- lapply(Z_raw, function(Z) Z[common_cpgs, common_samples, drop = FALSE])

## Compute bulk anchors #-----------
mu_bulk <- rowMeans(Xb, na.rm = TRUE)
sd_bulk <- apply(Xb, 1, sd, na.rm = TRUE)
sd_bulk[!is.finite(sd_bulk) | sd_bulk < 1e-8] <- 1

# Save for possible reuse
saveRDS(list(mu = mu_bulk, sd = sd_bulk, cpgs = common_cpgs), "bulk_anchors.rds")

## Standardize latent signals with frozen anchors #--------
Z_std <- lapply(Z_al, function(Zh) sweep(sweep(Zh, 1, mu_bulk, `-`), 1, sd_bulk, `/`))

save(Z_std, file = file.path(dirOut, "tca.tensor.z.Rdata"))







