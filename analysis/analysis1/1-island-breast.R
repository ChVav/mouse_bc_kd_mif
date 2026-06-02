  
# Run lm for 2x2x2 design across all island CpG M-values for the breast/cervix
# extract significant features
# posterior probability calculation for the significant features

library(here)
library(tidyverse)
library(patchwork)
library(limma)
library(showtext) # so arrows and delta are rendered properly

showtext_auto()

source(here("src/helpers_exposure-modifier.R"))

results_dir <- "results"
dirOut <- file.path(results_dir,"analysis1/1-output")
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

# # Mostly unmethylated islands
# beta_long <- as.data.frame(beta_sub) %>%
#   mutate(CpG = rownames(beta_sub)) %>%
#   pivot_longer(-CpG, names_to="Sample", values_to="Beta")
# 
# ggplot(beta_long, aes(x=Beta, color=Sample)) +
#   geom_density() +
#   theme_minimal() +
#   xlim(0,1) +
#   labs(title="Beta value densities per sample", x="Beta", y="Density") +
#   theme(legend.position="bottom", legend.text=element_text(size=6))

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

outcomes <- intersect(cpg_island,rownames(beta_sub))

# Linear contrasts for different diets #----
beta_selected <- t(beta_sub[rownames(beta_sub) %in% outcomes, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
outcomes <- paste0(outcomes,"_logit")
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))
dat <- full_join(dat,beta_selected)

## Benchmarking ###----

# test_outcomes <- outcomes[1:50]
# 
# beta_selected <- t(beta_sub[rownames(beta_sub) %in% test_outcomes, ]) %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "basename") %>%
#   logit_transform()
# test_outcomes <- paste0(test_outcomes,"_logit")
# beta_selected <- beta_selected %>% dplyr::select(basename,all_of(test_outcomes))
# dat <- full_join(dat,beta_selected)
# 
# library(microbenchmark)
# 
# microbenchmark(
#   lm_emmeans      = map_dfr(test_outcomes, ~summarize_outcome(.x, dat = dat, outcomes = test_outcomes)),
#   lm_fast         = map_dfr(test_outcomes, ~summarize_outcome_fast(.x, dat = dat)),
#   limma_vectorized = summarize_outcomes_limma(dat, outcome_cols = test_outcomes, moderate = TRUE),
#   times = 5
# )

# Unit: milliseconds
# expr       min        lq       mean    median        uq       max neval cld
# lm_emmeans 7071.4955 8007.4628 8062.16352 8231.5950 8366.3561 8633.9082     5 a  
# lm_fast 1384.3674 1432.2067 1452.78142 1459.0176 1492.8033 1495.5121     5  b 
# limma_vectorized   28.9509   37.0471   41.96366   37.3412   38.6799   67.7992     5   c

## Run linear model ###----

res1 <- summarize_outcomes_limma(
  dat = dat,
  outcome_cols = outcomes,
  moderate = T
)
# Warning message:
#   Partial NA coefficients for 19 probe(s) 

# Note this is very conservative. 
# res1 <- res1 %>%
#   dplyr::mutate(
#     p.adj = p.adjust(p.value, method = "BH")
#   ) 
#2718 out of 30,369 features are significant after FDR correction
#outcome_sig <- res1 %>% dplyr::filter(p.adj < 0.05 & keto == "KDminus" & mifepristone == "MIFminus") %>% pull(outcome)

# This looks ok
ggplot(res1, aes(x=beta_ref)) +
     geom_density() +
    theme_minimal() +
    xlim(0,1) +
     labs(title="Beta value densities", x="Beta", y="Density") 

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
#3176 out of 30,369 features are significant after FDR correction (< 0.05) for exposure alone (KD-, mife-)
#6005 , with FDR-adjusted P < 0.25

saveRDS(res1, file = file.path(dirOut, "breast_limma.Rds"))

res1_sig <- res1 %>% dplyr::filter(outcome %in% outcome_sig)
saveRDS(res1_sig, file = file.path(dirOut, "breast_limma_sig.Rds"))

## Top 50 hyper/hypo-methylated #----

res1_sig <- readRDS(file = file.path(dirOut, "breast_limma_sig.Rds"))

res1_sig <- res1_sig %>%
  dplyr::mutate(
    direction = ifelse(estimate > 0, "hyper", "hypo")
  )

# Top 50 hypermethylated
top_hyper <- res1_sig %>%
  filter(
    direction == "hyper",
    keto == "KD-",
    mifepristone == "MIF-"
  ) %>%
  arrange(desc(estimate)) %>%
  slice_head(n = 50) %>%
  pull(outcome)
res1_sig_top_hyper <- res1_sig %>% dplyr::filter(outcome %in% top_hyper)
top_hypo <- res1_sig %>%
  filter(
    direction == "hypo",
    keto == "KD-",
    mifepristone == "MIF-"
  ) %>%
  arrange(desc(abs(estimate))) %>%
  slice_head(n = 50) %>%
  pull(outcome)
res1_sig_top_hypo <- res1_sig %>% dplyr::filter(outcome %in% top_hypo)

res1_sig_top <- bind_rows(res1_sig_top_hypo,res1_sig_top_hyper)
res1_sig_top$outcome <- factor(res1_sig_top$outcome, levels = c(top_hypo, top_hyper))

pdf(file.path(dirOut,"breast_top100_lm-nomod.pdf"), width = 15, height = 10)
plot_heat_delta(res1_sig_top)
dev.off()

# Modulation #----

outcomes <- res1_sig %>% pull(outcome) %>% unique()
keep_cpg <- gsub("_logit","",outcomes)

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

beta_selected <- t(beta_sub[rownames(beta_sub) %in% keep_cpg, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))
dat <- full_join(dat,beta_selected)

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
  map_dfr(outcomes, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd))
}) %>%
  mutate(modifier = rep(modifiers, each = length(outcomes))) %>%
  mutate(
    color_value = prob_strengthen - prob_attenuate  # -1 → attenuate, +1 → strengthen
  )
saveRDS(res2, file = file.path(dirOut, "breast_posterior_sig.Rds"))

# This seems a tat too complicated
# synergy_rows <- compute_synergy_block(
#   outcomes = outcomes,
#   dat_base = dat_base,
#   dat_KD   = dat_KD,
#   dat_MIF  = dat_MIF,
#   dat_KD_MIF = dat_KD_MIF,
#   prior_sd = prior_sd,
#   orient_by_base = TRUE,   # set to FALSE for raw (unoriented) synergy
#   ndraws = NULL            # set to e.g., 4000 for draw-based CIs
# ) 
# saveRDS(synergy_rows, file = file.path(dirOut, "breast_posterior_synergy_sig.Rds"))# This is very hard to interpret, leave out?

# Combine results and annotate #----

res1 <- readRDS(file.path(dirOut, "breast_limma_sig.Rds")) %>%
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
res2 <- readRDS(file = file.path(dirOut, "breast_posterior_sig.Rds")) %>%
  as.data.frame()
res2$modifier <- factor(gsub("KD_MIF","KD + MIF",res2$modifier), levels = c("MIF", "KD", "KD + MIF"))

res <- left_join(res2,res1, by = c("outcome","modifier" = "cell"))
res$Name <- gsub("_logit","",res$outcome)

saveRDS(res, file = file.path(dirOut, "breast_results.Rds"))

anno <- anno %>% dplyr::filter(Name %in% res$Name)
saveRDS(anno, file = file.path(dirOut, "breast_sig_anno.Rds"))

