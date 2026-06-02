
library(tidyverse)

# Inverse-Variance Weighted Gene-Level Meta-Analysis 
# Calculate weighted effect size
# Optionally correct for CpG count by scaling SE

gene_meth_ipw <- function(dmr_result, gene_col = "gene_primary", group_cols = gene_col, adjust_n_cpg = FALSE){
  
  gene_meta <- dmr_result %>%
    filter(!is.na(.data[[gene_col]]), !is.na(estimate), !is.na(SE)) %>%
    mutate(weight = 1 / (SE^2)) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_cpg = n(),
      # IVW estimate
      mean_deltaM_gene = sum(weight * estimate) / sum(weight),
      # Gene-level SE
      SE_gene_raw = sqrt(1 / sum(weight)),
      .groups = "drop"
    )
  
  if(adjust_n_cpg == TRUE) {
    gene_meta <- gene_meta %>%
      mutate(
        # CpG-count adjustment (penalizes genes with many CpGs)
        SE_gene = SE_gene_raw * sqrt(n_cpg),
        z = mean_deltaM_gene / SE_gene,
        p_gene = 2 * pnorm(-abs(z))) %>%
      dplyr::select(-SE_gene_raw)
  } else{
    gene_meta <- gene_meta %>%
      dplyr::rename(SE_gene = SE_gene_raw) %>%
      mutate(
        z = mean_deltaM_gene / SE_gene,
        p_gene = 2 * pnorm(-abs(z)))
  }
  
  return(gene_meta)
  
}


# Analytical posterior for gene-level summaries #---------

posterior_delta_gene <- function(dat_base, dat_mod, prior_sd = 0.5) {
  
  # Join baseline and modifier by gene_primary
  dat <- dat_base %>%
    dplyr::select(gene_primary, mean_deltaM_gene, SE_gene, promoter) %>%
    rename(base_mean = mean_deltaM_gene,
           base_se   = SE_gene) %>%
    inner_join(
      dat_mod %>%
        select(gene_primary, mean_deltaM_gene, SE_gene, promoter) %>%
        rename(mod_mean = mean_deltaM_gene,
               mod_se   = SE_gene),
      by = c("gene_primary", "promoter") # Would have to adjust if different regions examined
    )
  
  # Analytical posterior per gene
  dat %>%
    rowwise() %>%
    mutate(
      # Posterior variance & mean for baseline
      post_var_base  = 1 / (1/prior_sd^2 + 1/base_se^2),
      post_mean_base = post_var_base * (0/prior_sd^2 + base_mean/base_se^2),
      
      # Posterior variance & mean for modifier
      post_var_mod   = 1 / (1/prior_sd^2 + 1/mod_se^2),
      post_mean_mod  = post_var_mod * (0/prior_sd^2 + mod_mean/mod_se^2),
      
      # Delta: modifier − baseline
      direction      = sign(post_mean_base),
      delta_mean     = direction * (post_mean_mod - post_mean_base),
      delta_sd       = sqrt(post_var_base + post_var_mod),
      
      # Probability of attenuation vs strengthening
      prob_attenuate  = pnorm(0, mean = delta_mean, sd = delta_sd),
      prob_strengthen = 1 - prob_attenuate,
      
      delta_l95 = delta_mean - 1.96 * delta_sd,
      delta_u95 = delta_mean + 1.96 * delta_sd
    ) %>%
    ungroup() %>%
    dplyr::select(
      gene_primary, promoter,
      delta_mean, delta_sd,
      delta_l95, delta_u95,
      prob_attenuate, prob_strengthen
    )
}
