
# Run lm for 2x2x2 design across all island CpG M-values for the breast/cervix
# Correct for cell types, by adding 2 PCAs on cell subtypes as covariates
# extract significant features
# posterior probability calculation for the significant features
# About half of the DMPs retained
# These highlight cell-type-specific regulation (or may still be driven by small compositonal changes not captured by PCA)

library(here)
library(tidyverse)
library(patchwork)
library(limma)
library(showtext) # so arrows and delta are rendered properly
showtext_auto()
library(reshape2)

source(here("src/helpers_exposure-modifier.R"))
source(here("src/clr_ilr_transform.R"))
source(here("src/plot_box_raw.R"))

diet4_colors <- c('#009392','#9ccb86','#eeb479','#cf597e') 

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis2_celltype_adj/2-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

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

# Combine fat and fibroblasts
dat$EpiFibFatIc_Stromal <- dat$EpiFibFatIc_Fat + dat$EpiFibFatIc_Fib
dat <- dat[, !colnames(dat) %in% c("EpiFibFatIc_Fat", "EpiFibFatIc_Fib")]

outcomes <- intersect(cpg_island,rownames(beta_sub))

# PCA analysis cell types #----
# covars <- colnames(dat)[grepl("EpiFibFatIc",colnames(dat))]
# covar1 <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Fib","EpiFibFatIc_Fat")
# covar2 <- setdiff(covars, covar1[1:2])
# clr_df <- clr_transform(dat,covar2, c("basename","exposure","keto", "mifepristone", "diet4"))

celltype_cols <- colnames(dat)[grepl("EpiFibFatIc",colnames(dat))]
covars <- c("EpiFibFatIc_Epi","EpiFibFatIc_IC","EpiFibFatIc_Stromal")
covar2 <- setdiff(celltype_cols, covars[1:2])

clr_df <- clr_transform(dat,covar2, c("basename","exposure","keto", "mifepristone", "diet4"))

clr_mat <- clr_transform(dat,covar2, c("basename")) %>%
  column_to_rownames(var = "basename")

pca <- prcomp(clr_mat, center = TRUE, scale. = FALSE) # already scaled cause of clr
save(pca, file = file.path(dirOut, "pca.Rdata"))
s <- summary(pca) # PCA1 and PC2 capture 67 and 13% of variation
pc_scores <- pca$x %>% as.data.frame() %>% rownames_to_column(var = "basename")

pdat <- clr_df %>% 
  left_join(pc_scores) 

# Plot
pc_plot <- pdat %>%
  ggplot(aes(x = PC1, y = PC2, shape = exposure, color = diet4)) +
  geom_point() +
  scale_color_manual(values = diet4_colors) +
  theme_classic()

pdat <- clr_df %>%
  dplyr::select(-keto, -mifepristone) %>%
  dplyr::select(-any_of(covar2)) %>%
  left_join(pc_scores) %>%
  melt(id.vars = c("basename","exposure","diet4"), variable.name = "PC", value.name = "score")

relabel <- paste0(names(s$importance[2,])," (", round(unname(s$importance[2,])*100, digits = 0), "%)") 
names(relabel) <- names(s$importance[2,])

pc_plot <- pdat %>%
  dplyr::filter(PC %in% c("PC1","PC2","PC3","PC4")) %>%
  ggplot(aes(x = exposure, y = score, group = exposure, shape = exposure, color = diet4)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  facet_wrap(~PC, scales = "free_y", labeller = labeller(PC = relabel), ncol = 4) +
  scale_color_manual(values = diet4_colors) +
  labs(x = "Exposure", y = "PC score", color = "", shape = "") +
  theme_classic()

# Check loadings to find which cell type changes are dominant
loadings <- pca$rotation  

# Sort contributing cell types per PC
PC1_top <- sort(abs(loadings[, "PC1"]), decreasing = TRUE)
PC2_top <- sort(abs(loadings[, "PC2"]), decreasing = TRUE)
PC3_top <- sort(abs(loadings[, "PC3"]), decreasing = TRUE)
PC4_top <- sort(abs(loadings[, "PC4"]), decreasing = TRUE)

#Create a data frame for plotting
top_df <- data.frame(
  CellType = rep(names(c(PC1_top, PC2_top, PC3_top, PC4_top)), 2),
  Loading = c(as.numeric(PC1_top), as.numeric(PC2_top), as.numeric(PC3_top), as.numeric(PC4_top)),
  PC = rep(c("PC1","PC2","PC3","PC4"), each = 10)
)

# Keep sign for coloring
top_df$Sign <- ifelse(c(loadings[names(PC1_top), "PC1"], loadings[names(PC2_top), "PC2"],
                        loadings[names(PC1_top), "PC3"], loadings[names(PC2_top), "PC4"]) > 0, "Positive", "Negative")


top_df$CellType <- gsub("EpiFibFatIc_","",top_df$CellType)

# Barplot
bar_plot <- ggplot(top_df, aes(x = reorder(CellType, Loading), y = Loading, fill = Sign)) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC, scales = "free_x", ncol = 4 ) +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#1b9e77", "Negative" = "#d95f02")) +
  labs(x = "Cell type", y = "Absolute Loading") +
  theme_classic()

pdf(file.path(dirOut, "pca_celltypes.pdf"), height = 8, width = 16)
print(pc_plot/bar_plot)
dev.off()

# Linear contrasts for different diets #----
beta_selected <- t(beta_sub[rownames(beta_sub) %in% outcomes, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
outcomes <- paste0(outcomes,"_logit")
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))

# include top 4 PCs to correct for celltype differences
dat_covar1 <- full_join(dat, pc_scores[,1:5]) %>%
  full_join(beta_selected)
covars <- colnames(pc_scores[,1:5])[-1]

# save PCs for reuse
saveRDS(pc_scores, file = file.path(dirOut,"celltype_pc_scores.Rds"))

## Run linear model ###----

res1 <- summarize_outcomes_limma(
  dat = dat_covar1,
  outcome_cols = outcomes,
  moderate = T,
  covars = covars
)

# Do multiple testing rather only for baseline and diet groups separately
res1 <- res1 %>%
  group_by(keto, mifepristone) %>%
  mutate(p.adj_cell = p.adjust(p.value, "BH")) %>%
  ungroup()

res1$keto <- gsub("minus","-",res1$keto)
res1$keto <- gsub("plus","+",res1$keto)
res1$mifepristone <- gsub("minus","-",res1$mifepristone)
res1$mifepristone <- gsub("plus","+",res1$mifepristone)

outcome_sig <- res1 %>%
  filter(keto == "KD-", mifepristone == "MIF-", p.adj_cell < 0.05) %>%
  tidyr::drop_na() %>%
  pull(outcome) %>% unique()
#1439 out of 30,369 features are significant after FDR correction (< 0.05) for exposure alone (KD-, mife-) # if I correct for 3 out 4 major cell types directly got only 4
#4022, with FDR-adjusted P < 0.25

saveRDS(res1, file = file.path(dirOut, "breast_limma_celltype_adj.Rds"))

res1_sig <- res1 %>% dplyr::filter(outcome %in% outcome_sig)
saveRDS(res1_sig, file = file.path(dirOut, "breast_limma_sig_celltype_adj.Rds"))

# Modulation #----

outcomes <- res1_sig %>% pull(outcome) %>% unique()
keep_cpg <- gsub("_logit","",outcomes)

beta_selected <- t(beta_sub[rownames(beta_sub) %in% keep_cpg, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))
dat <- full_join(dat, pc_scores[,1:5]) %>%
  full_join(beta_selected)

modifiers <- c("KD","MIF","KD_MIF")

# Subset datasets per modifier
dat_base <- dat %>% filter(keto=="KD-" & mifepristone=="MIF-")
dat_KD   <- dat %>% filter(keto=="KD+" & mifepristone=="MIF-")
dat_MIF  <- dat %>% filter(keto=="KD-" & mifepristone=="MIF+")
dat_KD_MIF <- dat %>% filter(keto=="KD+" & mifepristone=="MIF+")

dat_mod_list <- list(KD = dat_KD, MIF = dat_MIF, KD_MIF = dat_KD_MIF)

# Define prior SD
deltaM_baseline <- res1 %>% dplyr::filter(keto == "KD-" & mifepristone == "MIF-") %>% pull(estimate) %>% na.omit()
prior_sd <- mad(deltaM_baseline, constant = 1)  # robust SD

# Loop over modifiers and outcomes to calculate directional posterior probability
res2 <- map_dfr(modifiers, function(mod) {
  dat_mod <- dat_mod_list[[mod]]
  map_dfr(outcomes, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd, covars = covars))
}) %>%
  mutate(modifier = rep(modifiers, each = length(outcomes))) %>%
  mutate(
    color_value = prob_strengthen - prob_attenuate  # -1 → attenuate, +1 → strengthen
  )
saveRDS(res2, file = file.path(dirOut, "breast_posterior_sig_celltype_adj.Rds"))

# Combine results and annotate #----

res1 <- readRDS(file.path(dirOut, "breast_limma_sig_celltype_adj.Rds")) %>%
  as.data.frame() %>%
  mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
  mutate(cell = case_when(
    cell == "KD- | MIF-" ~ "No modifier",
    cell == "KD- | MIF+" ~ "MIF",
    cell == "KD+ | MIF-" ~ "KD",
    cell == "KD+ | MIF+" ~ "KD + MIF"
  )) %>%
  dplyr::mutate(
    direction = ifelse(estimate > 0, "hyper", "hypo")
  )
res2 <- readRDS(file = file.path(dirOut, "breast_posterior_sig_celltype_adj.Rds")) %>%
  as.data.frame()
res2$modifier <- factor(gsub("KD_MIF","KD + MIF",res2$modifier), levels = c("MIF", "KD", "KD + MIF"))

res <- left_join(res2,res1, by = c("outcome","modifier" = "cell"))
res$Name <- gsub("_logit","",res$outcome)

saveRDS(res, file = file.path(dirOut, "breast_results.Rds"))

anno <- anno %>% dplyr::filter(Name %in% res$Name)
saveRDS(anno, file = file.path(dirOut, "breast_sig_anno.Rds"))
