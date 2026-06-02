library(tidyverse)
library(patchwork)
library(emmeans)
library(limma)
library(scales)
library(tidytext)

# Helper for M-values
logit_transform <- function(df) {
  df %>%
    mutate(
      across(
        .cols = starts_with("cg") & where(is.numeric),
        .fns = ~ log2(.x / (1 - .x)),
        .names = "{str_replace(col, '_beta$', '')}_logit"
      )
    )
}

# Helper to compute exposure contrasts per outcome
summarize_outcome <- function(y,dat, outcomes, covars = NULL) {

  # fit lm
  if (!is.null(covars)) {
    covar_formula <- paste(covars, collapse = " + ")
    f <- as.formula(
      paste0("`", y, "` ~ exposure * keto * mifepristone + ", covar_formula)
    )
  } else {
    f <- as.formula(
      paste0("`", y, "` ~ exposure * keto * mifepristone")
    )
  }
  fit <- lm(f, data = dat)

  # EMMs for exposure within each keto × mife cell
  emm_E_by_KM <- emmeans(fit, ~ exposure | keto * mifepristone)

  # Exposure effect Δ = (P/D+) − (P/D-) per cell
  delta <- contrast(emm_E_by_KM, method = "revpairwise",
                    by = c("keto","mifepristone"))

  # Tidy summary (unadjusted p-values)
  s <- as.data.frame(summary(delta, infer = c(TRUE, TRUE)))
  s$outcome <- y
  s$sigma <- sigma(fit)
  s$std_effect_d <- s$estimate / s$sigma   # quick Cohen's d (using model σ)
  s # Cancer delta per treatment
}

summarize_outcome_comp <- function(y, dat_list, celltypes = c("Epi","Imm","Stro"), covars = NULL) {
  
  res_list <- lapply(celltypes, function(ct) {
    dat <- dat_list[[ct]]
    
    # for multiple outcomes
    map_dfr(y, function(outc) {
      
      # build formula
      if (!is.null(covars)) {
        covar_formula <- paste(covars, collapse = " + ")
        f <- as.formula(paste0("`", outc, "` ~ exposure + ", covar_formula))
      } else {
        f <- as.formula(paste0("`", outc, "` ~ exposure"))
      }
      
      fit <- lm(f, data = dat)
      
      emm_E <- emmeans(fit, ~ exposure)
      delta <- contrast(emm_E, method = "revpairwise")
      
      s <- as.data.frame(summary(delta, infer = c(TRUE, TRUE)))
      s$outcome <- outc
      s$celltype <- ct
      s$sigma <- sigma(fit)
      s$std_effect_d <- s$estimate / s$sigma
      
      s
    })
  })
  
  do.call(rbind, res_list)
}

# Heatmap of Δ
p_to_star <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

summarize_outcomes_limma <- function(dat,
                                     outcome_cols,
                                     moderate = FALSE,
                                     conf.level = 0.95,
                                     covars = NULL) {
  
  library(limma)
  library(dplyr)
  
  # ---- 1. Build outcome matrix (features x samples) ----
  dat$sample_id <- seq_len(nrow(dat))
  z_mat <- as.matrix(dat[, outcome_cols])
  rownames(z_mat) <- dat$sample_id
  z_mat <- t(z_mat)
  
  # ---- 2. Design matrix (exact same model as lm) ----
  if (!is.null(covars)) {
    covar_formula <- paste(covars, collapse = " + ")
    form <- as.formula(
      paste("~ exposure * keto * mifepristone +", covar_formula)
    )
  } else {
    form <- ~ exposure * keto * mifepristone
  }
  
  design <- model.matrix(form, data = dat)
  rownames(design) <- dat$sample_id
  
  fit <- lmFit(z_mat, design)
  
  coef_names <- colnames(design)
  
  # ---- 3. Build contrast matrix numerically (SAFE) ----
  C <- matrix(0, nrow = length(coef_names), ncol = 4)
  rownames(C) <- coef_names
  colnames(C) <- c("KDminus_MIFminus",
                   "KDplus_MIFminus",
                   "KDminus_MIFplus",
                   "KDplus_MIFplus")
  
  # Reference cell: KD- MIF-
  C["exposureP/D+", "KDminus_MIFminus"] <- 1
  
  # KD+ MIF-
  C["exposureP/D+", "KDplus_MIFminus"] <- 1
  C["exposureP/D+:ketoKD+", "KDplus_MIFminus"] <- 1
  
  # KD- MIF+
  C["exposureP/D+", "KDminus_MIFplus"] <- 1
  C["exposureP/D+:mifepristoneMIF+", "KDminus_MIFplus"] <- 1
  
  # KD+ MIF+
  C["exposureP/D+", "KDplus_MIFplus"] <- 1
  C["exposureP/D+:ketoKD+", "KDplus_MIFplus"] <- 1
  C["exposureP/D+:mifepristoneMIF+", "KDplus_MIFplus"] <- 1
  C["exposureP/D+:ketoKD+:mifepristoneMIF+", "KDplus_MIFplus"] <- 1
  
  fit2 <- contrasts.fit(fit, C)
  
  if (moderate) {
    fit2 <- eBayes(fit2)
  }
  
  df <- if (moderate) fit2$df.total else fit2$df.residual
  
  alpha <- 1 - conf.level
  tcrit <- qt(1 - alpha/2, df)
  
  # ---- 4. Identify reference samples ----
  beta_fun <- function(M) 2^M / (1 + 2^M)
  ref_samples <- dat$sample_id[
    dat$exposure == "P/D-" & dat$keto == "KD-" & dat$mifepristone == "MIF-"
  ]
  M_ref <- rowMeans(z_mat[, ref_samples])       # true baseline M
  beta_ref <- beta_fun(M_ref)                   # true baseline β
  
  # ---- 5. Build output ----
  results <- lapply(seq_len(ncol(C)), function(i) {
    
    est <- fit2$coefficients[, i]
    se  <- fit2$stdev.unscaled[, i] * fit2$sigma
    t   <- est / se
    p   <- 2 * pt(abs(t), df, lower.tail = FALSE)
    
    lwr <- est - tcrit * se
    upr <- est + tcrit * se
    
    delta_beta <- beta_fun(M_ref + est) - beta_ref  # biologically interpretable Δβ
    
    data.frame(
      outcome = rownames(z_mat),
      keto = sub("_.*", "", colnames(C)[i]),
      mifepristone = sub(".*_", "", colnames(C)[i]),
      estimate = est, # limma modeled deltaM
      SE = se,
      df = df,
      t.ratio = t,
      p.value = p,
      lower.CL = lwr,
      upper.CL = upr,
      sigma = fit2$sigma,
      std_effect_d = est / fit2$sigma,
      delta_beta = delta_beta,
      beta_ref = beta_ref,
      row.names = NULL
    )
    
  })
  
  bind_rows(results)
}

plot_heat_delta <- function(res){
  
  res <- res %>%
    mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    mutate(cell = case_when(
        cell == "KD- | MIF-" ~ "Ctrl",
        cell == "KD+ | MIF-" ~ "KD",
        cell == "KD- | MIF+" ~ "MIF",
        cell == "KD+ | MIF+" ~ "KD + MIF"))
  res$cell <- factor(res$cell, levels = c("Ctrl", "KD", "MIF",  "KD + MIF"))
  
  p <- ggplot(res, aes(x = cell, y = outcome, fill = estimate)) +
    geom_tile(color = "grey85") +
    scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
      breaks = pretty_breaks(n = 3)
    ) +
    geom_text(aes(label = p_to_star(p.value)), size = 5) +
    theme_minimal(base_size = 10) +
    labs(title = "ΔExposure",
         fill = "LE") +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(vjust = 0.8),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    xlab("") +
    ylab("")
     
  return(p)
}

# Analytical posterior helper ----
posterior_delta <- function(y, dat_base, dat_mod, prior_sd = 0.5, covars = NULL) {
  
  # Build formulas
  if (!is.null(covars)) {
    covar_formula <- paste(covars, collapse = " + ")
    f <- as.formula(
      paste0("`", y, "` ~ exposure + ", covar_formula)
    )
  } else {
    f <- as.formula(paste0("`", y, "` ~ exposure"))
  }
  
  # Fit OLS to get estimates
  lm_base <- lm(f, data = dat_base)
  lm_mod  <- lm(f, data = dat_mod)
  
  # Extract point estimate and standard error
  b_base  <- coef(lm_base)["exposureP/D+"]
  se_base <- summary(lm_base)$coefficients["exposureP/D+", "Std. Error"]
  
  b_mod   <- coef(lm_mod)["exposureP/D+"]
  se_mod  <- summary(lm_mod)$coefficients["exposureP/D+", "Std. Error"]
  
  # Analytical posterior: Normal conjugate with prior N(0, prior_sd^2)
  post_var_base <- 1 / (1/prior_sd^2 + 1/se_base^2)
  post_mean_base <- post_var_base * (0/prior_sd^2 + b_base/se_base^2)
  
  post_var_mod <- 1 / (1/prior_sd^2 + 1/se_mod^2)
  post_mean_mod <- post_var_mod * (0/prior_sd^2 + b_mod/se_mod^2)
  
  # δ = modifier − baseline
  direction <- sign(post_mean_base) # direction of baseline
  delta_mean <- direction * (post_mean_mod - post_mean_base)
  delta_sd   <- sqrt(post_var_mod + post_var_base)  # sum of independent Normals
  
  prob_attenuate  <- pnorm(0, mean = delta_mean, sd = delta_sd)
  prob_strengthen <- 1 - prob_attenuate
  
  tibble(
    outcome = y,
    delta_mean = delta_mean,
    delta_l95  = delta_mean - 1.96*delta_sd,
    delta_u95  = delta_mean + 1.96*delta_sd,
    prob_attenuate = prob_attenuate,
    prob_strengthen = prob_strengthen
  )
}

# Heatmap of posterior probablities #----

prob_to_arrow <- function(p, direction = c("up","down")) {
  direction <- match.arg(direction)
  
  n_arrow <- case_when(
    p <= 0.7 ~ 0,       # corresponds to P = 0.85, strong certainty
    p <= 0.85 ~ 1, # P = 0.925, very strong certainty
    p <= 0.95 ~ 2, # P = 0.975, near-certain
    TRUE ~ 3
  )
  
  arrow_char <- if(direction=="up") "\u2191" else "\u2193"  # ↑ or ↓
  
  # Vectorized: paste arrows, 0 → empty string
  sapply(n_arrow, function(n) if(n==0) "" else paste0(rep(arrow_char, n), collapse=""))
}

plot_heat_bayes <- function(res){
  
  res$cell <- factor(res$modifier, levels = c("KD","MIF", "KD_MIF"),
                     labels = c("KD","MIF",  "KD + MIF"))
  
  p <- ggplot(res, aes(x=cell, y=outcome, fill=color_value)) + # color = 2P(delta > 0) -1
    geom_tile(color="grey85") +
    geom_text(aes(label=paste0(
      prob_to_arrow(prob_attenuate, "down"), #P(delta < 0)
      prob_to_arrow(prob_strengthen, "up") # P(delta > 0)
    )), size=3) +
    scale_fill_gradient2(
      low="#6a00a8", mid="white", high="#ff7f0e",
      midpoint=0,
      breaks = pretty_breaks(n = 3)
    ) +
    theme_minimal(base_size=10) +
    labs(title = "ΔTreatment",
         fill = "PDC") +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(vjust = 0.8),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    xlab("") +
    ylab("") 

  return(p)
  
}

# Combined heatmap #----

plot_heat_delta_bayes <- function (res1, res2){
  
  res1 <- res1 %>%
    mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    mutate(cell = case_when(
      cell == "KD- | MIF-" ~ "Ctrl", TRUE ~ cell))
  res1$cell <- factor(res1$cell, levels = c("Ctrl", "KD+ | MIF-","KD- | MIF+", "KD+ | MIF+"))
  
  res2$cell <- factor(res2$modifier, levels = c("KD","MIF",  "KD_MIF"),
                      labels = c("KD","MIF",  "KD + MIF"))
  
  # Only show bayes prob when baseline modifier effect is significant
  baseline_sig <- res1 %>%
    filter(cell == "Ctrl") %>%
    transmute(outcome,
              baseline_significant = p.value < 0.05)
  res2 <- res2 %>%
    left_join(baseline_sig, by = "outcome") %>%
    mutate(
      color_value = ifelse(!baseline_significant, 0, color_value),
      prob_attenuate = ifelse(!baseline_significant, 0.5, prob_attenuate),
      prob_strengthen = ifelse(!baseline_significant, 0.5, prob_strengthen)
    )
  
  p_linear <- plot_heat_delta(res1)
  p_bayes <- plot_heat_bayes(res2) + theme(axis.text.y = element_blank())
  
  p_heat <- p_linear + p_bayes +
    plot_layout(widths = c(4, 3))
  
  return(p_heat)
}

# Wrappers #----

wrap_summary_plot <- function(pdat, outcomes, outcome_labels = NULL) {
  
  # ---- Scale outcomes for visualization ----
  pdat_scaled <- pdat %>%
    mutate(across(all_of(outcomes), ~ scale(.x, center = TRUE, scale = TRUE)))
  
  # ---- Linear summaries ----
  res1 <- map_dfr(outcomes, ~ summarize_outcome(.x, dat = pdat_scaled, outcomes = outcomes))
  
  # ---- Analytical posterior calculation ----
  modifiers <- c("KD","MIF","KD_MIF")
  
  dat_base <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF-")
  dat_KD   <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF-")
  dat_MIF  <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF+")
  dat_KD_MIF <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF+")
  
  dat_mod_list <- list(KD = dat_KD, MIF = dat_MIF, KD_MIF = dat_KD_MIF)
  
  # Define prior SD
  deltaM_baseline <- res1 %>% dplyr::filter(keto == "KD-" & mifepristone == "MIF-") %>% pull(estimate) %>% na.omit()
  prior_sd <- mad(deltaM_baseline, constant = 1)  # robust SD
  
  res2 <- map_dfr(modifiers, function(mod) {
    dat_mod <- dat_mod_list[[mod]]
    map_dfr(outcomes, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd))
  }) %>%
    mutate(modifier = rep(modifiers, each = length(outcomes))) %>%
    mutate(color_value = prob_strengthen - prob_attenuate)
  
  # ---- Handle outcome labels exactly like the multifacet function ----
  if (is.null(outcome_labels)) {
    outcome_labels <- outcomes
  }
  
  # relabel factors for plotting
  res1$outcome <- factor(res1$outcome,
                         levels = rev(outcomes),
                         labels = rev(outcome_labels))
  
  res2$outcome <- factor(res2$outcome,
                         levels = rev(outcomes),
                         labels = rev(outcome_labels))
  
  # ---- Define cells as factors ----
  res1 <- res1 %>%
    mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    mutate(cell = case_when(
      cell == "KD- | MIF-" ~ "Ctrl",
      cell == "KD+ | MIF-" ~ "KD",
      cell == "KD- | MIF+" ~ "MIF",
      cell == "KD+ | MIF+" ~ "KD + MIF"))
  
  res1$cell <- factor(res1$cell, levels = c("Ctrl", "KD", "MIF", "KD + MIF"))
  
  res2$cell <- factor(res2$modifier, levels = c("KD","MIF", "KD_MIF"),
                      labels = c("KD","MIF","KD + MIF"))
  
  # ---- Baseline filtering ----
  baseline_sig <- res1 %>%
    filter(cell == "Ctrl") %>%
    transmute(outcome,
              baseline_significant = p.value < 0.05)
  
  res2 <- res2 %>%
    left_join(baseline_sig, by = "outcome") %>%
    mutate(
      color_value = ifelse(!baseline_significant, 0, color_value),
      prob_attenuate = ifelse(!baseline_significant, 0.5, prob_attenuate),
      prob_strengthen = ifelse(!baseline_significant, 0.5, prob_strengthen)
    )
  
  # ---- Plot ----
  p_linear <- plot_heat_delta(res1) +
    theme(
      axis.text.y = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      strip.background.y = element_rect(fill = "#F5F5F5", colour = "black", linewidth = 0.5)
    )
  
  p_bayes <- plot_heat_bayes(res2) +
    scale_y_discrete(position = "right") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_blank(),
          axis.text.y.right = element_text(hjust = 0))
  
  # Combine linear and bayesian plots side-by-side
  p_final <- p_linear + p_bayes + plot_layout(widths = c(4,3))
  
  return(list(p_final, res1))
}

plot_heat_comp <- function(res){
  
  p <- ggplot(res, aes(x = celltype, y = outcome, fill = estimate)) +
    geom_tile(color = "grey85") +
    scale_fill_gradient2(
      low = "#005F73", mid = "white", high = "#E31A1C", midpoint = 0,
      breaks = pretty_breaks(n = 3) 
    ) +
    geom_text(aes(label = p_to_star(p.value)), size = 5) +
    theme_minimal(base_size = 10) +
    labs(title = "Celltype specificity",
         fill = "ΔLSC") +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(vjust = 0.8),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    xlab("") +
    ylab("")
  
  return(p)
  
}

wrap_summary_plot_comp <- function(pdat, pdat_comp, outcomes, outcome_labels = NULL) {
  
  # ---- Scale outcomes for visualization ----
  pdat_scaled <- pdat %>%
    mutate(across(all_of(outcomes), ~ scale(.x, center = TRUE, scale = TRUE)))
  
  # ---- Linear summaries ----
  res1 <- map_dfr(outcomes, ~ summarize_outcome(.x, dat = pdat_scaled, outcomes = outcomes))
  
  # ---- Analytical posterior calculation ----
  modifiers <- c("KD","MIF","KD_MIF")
  
  dat_base <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF-")
  dat_KD   <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF-")
  dat_MIF  <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF+")
  dat_KD_MIF <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF+")
  
  dat_mod_list <- list(KD = dat_KD, MIF = dat_MIF, KD_MIF = dat_KD_MIF)
  
  res2 <- map_dfr(modifiers, function(mod) {
    dat_mod <- dat_mod_list[[mod]]
    map_dfr(outcomes, ~ posterior_delta(.x, dat_base, dat_mod))
  }) %>%
    mutate(modifier = rep(modifiers, each = length(outcomes))) %>%
    mutate(color_value = prob_strengthen - prob_attenuate)
  
  # ---- Outcome labels ----
  if (is.null(outcome_labels)) {
    outcome_labels <- outcomes
  }
  
  res1$outcome <- factor(res1$outcome,
                         levels = rev(outcomes),
                         labels = rev(outcome_labels))
  
  res2$outcome <- factor(res2$outcome,
                         levels = rev(outcomes),
                         labels = rev(outcome_labels))
  
  # ---- Define treatment cells ----
  res1 <- res1 %>%
    mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    mutate(cell = case_when(
      cell == "KD- | MIF-" ~ "Ctrl",
      cell == "KD+ | MIF-" ~ "KD",
      cell == "KD- | MIF+" ~ "MIF",
      cell == "KD+ | MIF+" ~ "KD + MIF"))
  
  res1$cell <- factor(res1$cell, levels = c("Ctrl", "KD", "MIF", "KD + MIF"))
  
  res2$cell <- factor(res2$modifier, levels = c("KD","MIF", "KD_MIF"),
                      labels = c("KD","MIF","KD + MIF"))
  
  # ---- Baseline significance ----
  baseline_sig <- res1 %>%
    filter(cell == "Ctrl") %>%
    transmute(outcome,
              baseline_significant = p.value < 0.05)
  
  res2 <- res2 %>%
    left_join(baseline_sig, by = "outcome") %>%
    mutate(
      color_value = ifelse(!baseline_significant, 0, color_value),
      prob_attenuate = ifelse(!baseline_significant, 0.5, prob_attenuate),
      prob_strengthen = ifelse(!baseline_significant, 0.5, prob_strengthen)
    )
  
  # Baseline estimates per cell type ----
  # res3: column = cell type (Epi, Imm, Stro)
  celltypes <- c("Epi","Imm","Stro")
  
  res3 <- map_dfr(celltypes, function(ct) {
    summarize_outcome_comp(y = outcomes, dat_list = pdat_comp, celltypes = ct)
  })
  
  res3$outcome <- factor(res3$outcome,
                         levels = rev(outcomes),
                         labels = rev(outcome_labels))
  
  res3$celltype <- factor(res3$celltype, levels = celltypes)
  
  # ---- Plots ----
  p_linear <- plot_heat_delta(res1) +
    theme(
      axis.text.y = element_blank()
      # strip.placement = "outside",
      # strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      # strip.background.y = element_rect(fill = "#F5F5F5", colour = "black", linewidth = 0.5)
    )
  
  p_bayes <- plot_heat_bayes(res2) +
    scale_y_discrete(position = "right") +
    theme(axis.text.y = element_blank())
  
  # Add baseline per celltype heatmap
  p_baseline_ct <- plot_heat_comp(res3) +
    scale_y_discrete(position = "right") +
    theme(strip.placement = "outside",
          strip.text.y.left = element_blank(),
          axis.text.y.right = element_text(hjust = 0))
    
  
  # Combine all three side-by-side
  p_final <- p_linear + p_bayes + p_baseline_ct  + plot_layout(widths = c(4,3,3))
  
  return(p_final)
}


wrap_summary_plot_multi <- function(pdat_list, tissues, outcomes, outcome_labels = NULL) {
  
  stopifnot(length(pdat_list) == length(tissues))
  
  modifiers <- c("KD","MIF","KD_MIF")
  
  all_res1 <- list()
  all_res2 <- list()
  
  for(i in seq_along(pdat_list)) {
    
    pdat <- pdat_list[[i]]
    tissue_name <- tissues[i]
    
    # ---- Scale outcomes ----
    pdat_scaled <- pdat %>%
      mutate(across(all_of(outcomes), ~ scale(.x, center = TRUE, scale = TRUE)))
    
    # ---- Linear summaries ----
    res1 <- map_dfr(outcomes, 
                    ~ summarize_outcome(.x, dat = pdat_scaled, outcomes = outcomes)) %>%
      mutate(tissue = tissue_name)
    
    # ---- Split modifier datasets ----
    dat_base <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF-")
    dat_KD   <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF-")
    dat_MIF  <- pdat %>% filter(keto=="KD-" & mifepristone=="MIF+")
    dat_KD_MIF <- pdat %>% filter(keto=="KD+" & mifepristone=="MIF+")
    
    dat_mod_list <- list(KD = dat_KD, 
                         MIF = dat_MIF, 
                         KD_MIF = dat_KD_MIF)
    
    # ---- Bayesian summaries ----
    # Define prior SD
    deltaM_baseline <- res1 %>% dplyr::filter(keto == "KD-" & mifepristone == "MIF-") %>% pull(estimate) %>% na.omit()
    prior_sd <- mad(deltaM_baseline, constant = 1)  # robust SD
    
    res2 <- map_dfr(modifiers, function(mod) {
      dat_mod <- dat_mod_list[[mod]]
      map_dfr(outcomes, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd))
    }) %>%
      mutate(modifier = rep(modifiers, each = length(outcomes))) %>%
      mutate(
        color_value = prob_strengthen - prob_attenuate,
        tissue = tissue_name
      )
    
    all_res1[[i]] <- res1
    all_res2[[i]] <- res2
  }
  
  res1_all <- bind_rows(all_res1)
  res2_all <- bind_rows(all_res2)
  
  # ---- Cell, outcome and tissue levels and labels #----
  res1_all <- res1_all %>%
    mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    mutate(cell = case_when(
      cell == "KD- | MIF-" ~ "Ctrl",
      cell == "KD+ | MIF-" ~ "KD",
      cell == "KD- | MIF+" ~ "MIF",
      cell == "KD+ | MIF+" ~ "KD + MIF"))
  res1_all$cell <- factor(res1_all$cell, 
                          levels = c("Ctrl", "KD", "MIF", "KD + MIF"))
  
  res2_all$cell <- factor(res2_all$modifier, levels = c("KD","MIF", "KD_MIF"),
                      labels = c("KD", "MIF", "KD + MIF"))
  
  if (is.null(outcome_labels)) {
    outcome_labels <- outcomes
  }
  
  res1_all$outcome <- factor(res1_all$outcome,
                             levels = rev(outcomes),
                             labels = rev(outcome_labels))
  res2_all$outcome <- factor(res2_all$outcome,
                             levels = rev(outcomes),
                             labels = rev(outcome_labels))
  
  res1_all$tissue <- factor(res1_all$tissue, levels = tissues)
  res2_all$tissue <- factor(res2_all$tissue, levels = tissues)
  
  # ---- Baseline filtering ----
  baseline_sig <- res1_all %>%
    filter(cell == "Ctrl") %>%
    transmute(outcome,
              tissue,
              baseline_significant = p.value < 0.05)
  
  res2_all <- res2_all %>%
    left_join(baseline_sig, by = c("outcome","tissue")) %>%
    mutate(
      color_value = ifelse(!baseline_significant, 0, color_value),
      prob_attenuate = ifelse(!baseline_significant, 0.5, prob_attenuate),
      prob_strengthen = ifelse(!baseline_significant, 0.5, prob_strengthen)
    )
  
  # ---- Plot ----
  p_linear <- plot_heat_delta(res1_all) +
    facet_grid(tissue ~ ., switch = "y") +
    theme(
      axis.text.y = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5
      ),
      strip.background.y = element_rect(
        fill = "#F5F5F5",
        colour = "black",
        linewidth = 0.5
      )
    )
  
  p_bayes <- plot_heat_bayes(res2_all) +
    facet_grid(tissue ~ ., switch = "y") +
    scale_y_discrete(position = "right") +
    theme(strip.placement = "outside", 
          strip.text.y.left = element_blank(),
          axis.text.y.right = element_text(hjust = 0))
  
  p_final <- p_linear + p_bayes +
    plot_layout(widths = c(4, 3))
  
  return(p_final)
}

# This version can deal with different outcomes per facet ("tissue")
wrap_summary_plot_multi2 <- function(pdat_list, tissues, outcomes, outcome_labels = NULL, drop_ns = TRUE) {
  stopifnot(length(pdat_list) == length(tissues))
  if (is.null(outcome_labels)) outcome_labels <- outcomes
  
  modifiers <- c("KD","MIF","KD_MIF")
  all_res1 <- list()
  all_res2 <- list()
  
  for (i in seq_along(pdat_list)) {
    pdat <- pdat_list[[i]]
    tissue_name <- tissues[i]
    
    # Outcomes available in this tissue
    outcomes_i <- intersect(outcomes, names(pdat))
    if (length(outcomes_i) == 0) {
      warning(sprintf("No specified outcomes present for tissue '%s'; skipping.", tissue_name))
      next
    }
    
    # Scale only available outcome columns
    pdat_scaled <- pdat %>%
      dplyr::mutate(across(dplyr::any_of(outcomes_i), ~ scale(.x, center = TRUE, scale = TRUE)))
    
    # Linear summaries (available outcomes only)
    res1 <- purrr::map_dfr(outcomes_i,
                           ~ summarize_outcome(.x, dat = pdat_scaled, outcomes = outcomes_i)) %>%
      dplyr::mutate(tissue = tissue_name)
    
    # Split datasets
    dat_base   <- pdat %>% dplyr::filter(keto=="KD-" & mifepristone=="MIF-")
    dat_KD     <- pdat %>% dplyr::filter(keto=="KD+" & mifepristone=="MIF-")
    dat_MIF    <- pdat %>% dplyr::filter(keto=="KD-" & mifepristone=="MIF+")
    dat_KD_MIF <- pdat %>% dplyr::filter(keto=="KD+" & mifepristone=="MIF+")
    
    dat_mod_list <- list(KD = dat_KD, MIF = dat_MIF, KD_MIF = dat_KD_MIF)
    
    # Bayesian summaries (available outcomes only)
    # Define prior SD
    deltaM_baseline <- res1 %>% dplyr::filter(keto == "KD-" & mifepristone == "MIF-") %>% pull(estimate) %>% na.omit()
    prior_sd <- mad(deltaM_baseline, constant = 1)  # robust SD
    
    res2 <- purrr::map_dfr(modifiers, function(mod) {
      dat_mod <- dat_mod_list[[mod]]
      purrr::map_dfr(outcomes_i, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd))
    }) %>%
      dplyr::mutate(
        modifier    = rep(modifiers, each = length(outcomes_i)),
        color_value = prob_strengthen - prob_attenuate,
        tissue      = tissue_name
      )
    
    all_res1[[i]] <- res1
    all_res2[[i]] <- res2
  }
  
  # Drop tissues that had no outcomes
  all_res1 <- Filter(Negate(is.null), all_res1)
  all_res2 <- Filter(Negate(is.null), all_res2)
  if (length(all_res1) == 0 || length(all_res2) == 0) {
    stop("No tissues with available outcomes to plot.")
  }
  
  res1_all <- dplyr::bind_rows(all_res1)
  res2_all <- dplyr::bind_rows(all_res2)
  
  # Cell labels
  res1_all <- res1_all %>%
    dplyr::mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    dplyr::mutate(cell = dplyr::case_when(
      cell == "KD- | MIF-" ~ "Ctrl",
      cell == "KD+ | MIF-" ~ "KD",
      cell == "KD- | MIF+" ~ "MIF",
      cell == "KD+ | MIF+" ~ "KD + MIF"
    ))
  res1_all$cell <- factor(res1_all$cell, levels = c("Ctrl", "KD", "MIF", "KD + MIF"))
  
  res2_all$cell <- factor(res2_all$modifier,
                          levels = c("KD","MIF","KD_MIF"),
                          labels = c("KD","MIF","KD + MIF"))
  
  # Outcome labels/order (global); facets will drop unused levels
  res1_all$outcome <- factor(res1_all$outcome,
                             levels = rev(outcomes),
                             labels = rev(outcome_labels))
  res2_all$outcome <- factor(res2_all$outcome,
                             levels = rev(outcomes),
                             labels = rev(outcome_labels))
  
  # Tissue order
  res1_all$tissue <- factor(res1_all$tissue, levels = tissues)
  res2_all$tissue <- factor(res2_all$tissue, levels = tissues)
  
  # Optional: filter to baseline-significant (p < 0.05) outcome × tissue pairs
  if (isTRUE(drop_ns)) {
    baseline_keep <- res1_all %>%
      dplyr::filter(cell == "Ctrl", !is.na(p.value)) %>%
      dplyr::group_by(outcome, tissue) %>%
      dplyr::summarise(baseline_significant = any(p.value < 0.05), .groups = "drop") %>%
      dplyr::filter(baseline_significant) %>%
      dplyr::select(outcome, tissue)
    
    res1_all <- dplyr::semi_join(res1_all, baseline_keep, by = c("outcome","tissue"))
    res2_all <- dplyr::semi_join(res2_all, baseline_keep, by = c("outcome","tissue"))
  }
  
  # Plots: each facet shows only the outcomes retained (no NA rows)
  p_linear <- plot_heat_delta(res1_all) +
    ggplot2::facet_grid(tissue ~ ., switch = "y", scales = "free_y", drop = TRUE) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      strip.background.y = ggplot2::element_rect(fill = "#F5F5F5", colour = "black", linewidth = 0.5)
    )
  
  p_bayes <- plot_heat_bayes(res2_all) +
    ggplot2::facet_grid(tissue ~ ., switch = "y", scales = "free_y", drop = TRUE) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_text(hjust = 0)
    )
  
  p_final <- p_linear + p_bayes + patchwork::plot_layout(widths = c(4, 3))
  return(p_final)
}

# This version can deal with different outcomes per facet ("tissue")
# Baseline delta in different celltypes ("compartment") are added (TCA decomposed meth values)
# Values ordered by bayes in KD cell
wrap_summary_plot_multi3 <- function(pdat_list, pdat_comp_list, tissues, outcomes, outcome_labels = NULL, covars = NULL, drop_ns = TRUE, min_rows = NULL) {
  # library(tidytext) # ensure this is available
  
  stopifnot(length(pdat_list) == length(tissues))
  stopifnot(length(pdat_comp_list) == length(tissues))
  if (is.null(outcome_labels)) outcome_labels <- outcomes
  
  modifiers <- c("KD","MIF","KD_MIF")
  celltypes <- c("Epi","Imm","Stro")
  all_res1 <- list()
  all_res2 <- list()
  all_res3 <- list()
  
  for (i in seq_along(pdat_list)) {
    pdat <- pdat_list[[i]]
    pdat2 <- pdat_comp_list[[i]]
    tissue_name <- tissues[i]
    
    # Outcomes available in this tissue
    outcomes_i <- intersect(outcomes, names(pdat))
    if (length(outcomes_i) == 0) {
      warning(sprintf("No specified outcomes present for tissue '%s'; skipping.", tissue_name))
      next
    }
    
    # Scale only available outcome columns
    pdat_scaled <- pdat %>%
      dplyr::mutate(across(dplyr::any_of(outcomes_i), ~ scale(.x, center = TRUE, scale = TRUE)))
    
    # Linear summaries (available outcomes only)
    res1 <- purrr::map_dfr(outcomes_i,
                           ~ summarize_outcome(.x, dat = pdat_scaled, outcomes = outcomes_i, covars = covars)) %>%
      dplyr::mutate(tissue = tissue_name)
    
    # Split datasets
    dat_base   <- pdat %>% dplyr::filter(keto=="KD-" & mifepristone=="MIF-")
    dat_KD     <- pdat %>% dplyr::filter(keto=="KD+" & mifepristone=="MIF-")
    dat_MIF    <- pdat %>% dplyr::filter(keto=="KD-" & mifepristone=="MIF+")
    dat_KD_MIF <- pdat %>% dplyr::filter(keto=="KD+" & mifepristone=="MIF+")
    
    dat_mod_list <- list(KD = dat_KD, MIF = dat_MIF, KD_MIF = dat_KD_MIF)
    
    # Bayesian summaries (available outcomes only)
    # Define prior SD
    deltaM_baseline <- res1 %>% dplyr::filter(keto == "KD-" & mifepristone == "MIF-") %>% pull(estimate) %>% na.omit()
    prior_sd <- mad(deltaM_baseline, constant = 1)  # robust SD
    
    res2 <- purrr::map_dfr(modifiers, function(mod) {
      dat_mod <- dat_mod_list[[mod]]
      purrr::map_dfr(outcomes_i, ~ posterior_delta(.x, dat_base, dat_mod, prior_sd = prior_sd, covars = covars))
    }) %>%
      dplyr::mutate(
        modifier    = rep(modifiers, each = length(outcomes_i)),
        color_value = prob_strengthen - prob_attenuate,
        tissue      = tissue_name
      )
    
    # Baseline per-cell-type summaries
    res3 <- purrr::map_dfr(celltypes, function(ct) {
      summarize_outcome_comp(y = outcomes_i, dat_list = pdat2, celltypes = ct)
    }) %>%
      dplyr::mutate(tissue = tissue_name)
    
    all_res1[[i]] <- res1
    all_res2[[i]] <- res2
    all_res3[[i]] <- res3
  }
  
  # Drop tissues that had no outcomes
  all_res1 <- Filter(Negate(is.null), all_res1)
  all_res2 <- Filter(Negate(is.null), all_res2)
  all_res3 <- Filter(Negate(is.null), all_res3)
  if (length(all_res1) == 0 || length(all_res2) == 0 || length(all_res3) == 0) {
    stop("No tissues with available outcomes to plot.")
  }
  
  res1_all <- dplyr::bind_rows(all_res1)
  res2_all <- dplyr::bind_rows(all_res2)
  res3_all <- dplyr::bind_rows(all_res3)
  
  # Cell labels
  res1_all <- res1_all %>%
    dplyr::mutate(cell = paste(keto, mifepristone, sep = " | ")) %>%
    dplyr::mutate(cell = dplyr::case_when(
      cell == "KD- | MIF-" ~ "Ctrl",
      cell == "KD+ | MIF-" ~ "KD",
      cell == "KD- | MIF+" ~ "MIF",
      cell == "KD+ | MIF+" ~ "KD + MIF"
    ))
  res1_all$cell <- factor(res1_all$cell, levels = c("Ctrl", "KD", "MIF", "KD + MIF"))
  
  res2_all$cell <- factor(res2_all$modifier,
                          levels = c("KD","MIF","KD_MIF"),
                          labels = c("KD","MIF","KD + MIF"))
  
  # Tissue and celltype order
  res1_all$tissue <- factor(res1_all$tissue, levels = tissues)
  res2_all$tissue <- factor(res2_all$tissue, levels = tissues)
  res3_all$tissue <- factor(res3_all$tissue, levels = tissues)
  res3_all$celltype <- factor(res3_all$celltype, levels = c("Epi","Imm","Stro"))
  
  # Optional: filter to baseline-significant (p < 0.05) outcome × tissue pairs
  if (isTRUE(drop_ns)) {
    baseline_keep <- res1_all %>%
      dplyr::filter(cell == "Ctrl", !is.na(p.value)) %>%
      dplyr::group_by(outcome, tissue) %>%
      dplyr::summarise(baseline_significant = any(p.value < 0.05), .groups = "drop") %>%
      dplyr::filter(baseline_significant) %>%
      dplyr::select(outcome, tissue)
    
    res1_all <- dplyr::semi_join(res1_all, baseline_keep, by = c("outcome","tissue"))
    res2_all <- dplyr::semi_join(res2_all, baseline_keep, by = c("outcome","tissue"))
    res3_all <- dplyr::semi_join(res3_all, baseline_keep, by = c("outcome","tissue"))
  }
  
  # Per-tissue outcome ordering by KD effect (descending)
  kd_order <- res2_all %>%
    dplyr::filter(cell == "KD") %>%
    dplyr::group_by(tissue) %>%
    dplyr::arrange(dplyr::desc(color_value), .by_group = TRUE) %>%  # DESC for "stronger at top"
    dplyr::mutate(order_idx = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(tissue, outcome, order_idx)
  
  # Apply per-facet ordering
  res1_plot <- res1_all %>%
    dplyr::left_join(kd_order, by = c("tissue","outcome")) %>%
    dplyr::mutate(outcome_chr = as.character(outcome),
                  order_idx   = ifelse(is.na(order_idx), 0L, order_idx),
                  outcome     = tidytext::reorder_within(outcome_chr, order_idx, tissue))
  
  res2_plot <- res2_all %>%
    dplyr::left_join(kd_order, by = c("tissue","outcome")) %>%
    dplyr::mutate(outcome_chr = as.character(outcome),
                  order_idx   = ifelse(is.na(order_idx), 0L, order_idx),
                  outcome     = tidytext::reorder_within(outcome_chr, order_idx, tissue))
  
  res3_plot <- res3_all %>%
    dplyr::left_join(kd_order, by = c("tissue","outcome")) %>%
    dplyr::mutate(outcome_chr = as.character(outcome),
                  order_idx   = ifelse(is.na(order_idx), 0L, order_idx),
                  outcome     = tidytext::reorder_within(outcome_chr, order_idx, tissue))
  
  # Optionally drop facets with fewer than min_rows outcomes
  if (!is.null(min_rows)) {
    tissue_keep <- res1_plot %>%
      dplyr::distinct(tissue, outcome_chr) %>%
      dplyr::count(tissue, name = "n_rows") %>%
      dplyr::filter(n_rows >= min_rows) %>%
      dplyr::pull(tissue) %>%
      unique()
    
    
    res1_plot <- res1_plot %>% dplyr::filter(tissue %in% tissue_keep)
    res2_plot <- res2_plot %>% dplyr::filter(tissue %in% tissue_keep)
    res3_plot <- res3_plot %>% dplyr::filter(tissue %in% tissue_keep)
    
    # Optional: fail early if nothing is left to plot
    if (nrow(res1_plot) == 0 || nrow(res2_plot) == 0 || nrow(res3_plot) == 0) {
      stop("All tissues had fewer than 'min_rows' unique outcomes; nothing to plot.")
    }
    
  }
  
  # Label mapping function
  label_fun <- function(x) {
    # strip reorder suffix "___<tissue>"
    base <- sub("___.*$", "", x)
    lbl_map <- setNames(outcome_labels, outcomes)
    out <- ifelse(base %in% names(lbl_map), lbl_map[base], base)
    unname(out)
  }
  
  # Plots using the reordered data
  p_linear <- plot_heat_delta(res1_plot) +
    ggplot2::facet_grid(tissue ~ ., switch = "y", scales = "free_y", space = "free_y", drop = TRUE) +
    tidytext::scale_y_reordered(labels = label_fun) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      strip.background.y = ggplot2::element_rect(fill = "#F5F5F5", colour = "black", linewidth = 0.5)
    )
  
  p_bayes <- plot_heat_bayes(res2_plot) +
    ggplot2::facet_grid(tissue ~ ., switch = "y", scales = "free_y", space = "free_y",  drop = TRUE) +
    tidytext::scale_y_reordered(labels = label_fun, position = "right") +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_blank()
    )
  
  p_baseline_ct <- plot_heat_comp(res3_plot) +
    ggplot2::facet_grid(tissue ~ ., switch = "y", scales = "free_y",space = "free_y",  drop = TRUE) +
    tidytext::scale_y_reordered(labels = label_fun, position = "right") +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_text(hjust = 0)
    )
  
  p_final <- p_linear + p_bayes + p_baseline_ct + patchwork::plot_layout(widths = c(4, 3, 3))
  return(p_final)
}



