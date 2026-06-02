
library(here)
library(tidyverse)
library(patchwork)
library(limma)
library(lme4)

source(here("src/helpers_exposure-modifier-synergy.R"))

results_dir <- here("results")
dirOut <- file.path(results_dir,"analysis1/2-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

# Annotation #----
anno <- readRDS(file.path(results_dir,"0-output/anno_12.v1.mm10.Rds"))
cpg_island <- anno %>% dplyr::filter(Relation_to_Island == "Island") %>% rownames()
cpg_nonisland <- anno %>% dplyr::filter(Relation_to_Island != "Island") %>% rownames()
promotor_cg <- anno %>% dplyr::filter(promoter == TRUE) %>% pull(Name)
body_cg <- anno %>% dplyr::filter(primary_feature == "tss_body") %>% pull(Name)

# Data load #----

load(here('data/beta_final.Rdata'))
load(here("data/pheno.Rdata"))

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

outcomes <- intersect(cpg_island,rownames(beta_sub))

# Linear contrasts for different diets #----
beta_selected <- t(beta_sub[rownames(beta_sub) %in% outcomes, ]) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "basename") %>%
  logit_transform()
outcomes <- paste0(outcomes,"_logit")
beta_selected <- beta_selected %>% dplyr::select(basename,all_of(outcomes))
dat <- full_join(dat,beta_selected)

## Run linear model ###----

res1 <- summarize_outcomes_limma(
  dat = dat,
  outcome_cols = outcomes,
  moderate = T
)

# Per CpG synergy results
res_synergy <- res1 %>%
  dplyr::filter(contrast == "Synergy")

sig_synergy <- res_synergy %>%
  dplyr::filter(p.value < 0.05) #509, but none after FDR

t.test(res_synergy$estimate)

# One Sample t-test
# 
# data:  res_synergy$estimate
# t = -42.804, df = 30297, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.07952370 -0.07255963
# sample estimates:
#   mean of x 
# -0.07604167 

saveRDS(res1, file = file.path(dirOut, "breast_limma_synergy.Rds"))
