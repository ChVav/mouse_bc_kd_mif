library(tidyverse)
library(patchwork)
library(limma)

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

summarize_outcomes_limma_simple <- function(dat,
                                     outcome_cols,
                                     moderate = FALSE,
                                     conf.level = 0.95,
                                     covars = NULL) {
  
  dat$sample_id <- seq_len(nrow(dat))
  
  # set reference level
  dat$diet4 <- relevel(factor(dat$diet4), "KD- MIF-")
  
  # outcome matrix
  z_mat <- as.matrix(dat[, outcome_cols])
  rownames(z_mat) <- dat$sample_id
  z_mat <- t(z_mat)
  
  # design matrix
  if (!is.null(covars)) {
    covar_formula <- paste(covars, collapse = " + ")
    form <- as.formula(
      paste("~ diet4 +", covar_formula)
    )
  } else {
    form <- ~ diet4
  }
  
  design <- model.matrix(form, data = dat)
  rownames(design) <- dat$sample_id
  
  fit <- lmFit(z_mat, design)
  
  if (moderate) {
    fit <- eBayes(fit)
  }
  
  # Degrees of freedom
  df <- if (moderate) fit$df.total else fit$df.residual
  alpha <- 1 - conf.level
  tcrit <- qt(1 - alpha/2, df)
  
  # Extract coefficients
  # Coefficients are:
  # (Intercept)          -> mean of KD- MIF-
  # diet4KD- MIF+        -> KD- MIF+ vs KD- MIF-
  # diet4KD+ MIF-        -> KD+ MIF- vs KD- MIF-
  # diet4KD+ MIF+        -> KD+ MIF+ vs KD- MIF-
  
  beta_fun <- function(M) 2^M / (1 + 2^M)
  ref_samples <- dat$sample_id[
    dat$keto == "KD-" & dat$mifepristone == "MIF-"
  ]
  M_ref <- rowMeans(z_mat[, ref_samples])       # true baseline M
  beta_ref <- beta_fun(M_ref)                   # true baseline β
  
  coef_names <- colnames(fit$coefficients)[-1]  # skip intercept, only want comparisons
  
  results <- lapply(coef_names, function(coef_name) {
    est <- fit$coefficients[, coef_name]
    se  <- fit$stdev.unscaled[, coef_name] * fit$sigma
    tval <- est / se
    pval <- 2 * pt(abs(tval), df, lower.tail = FALSE)
    
    lwr <- est - tcrit * se
    upr <- est + tcrit * se
    
    delta_beta <- beta_fun(M_ref + est) - beta_ref  # biologically interpretable Δβ
    
    data.frame(
      outcome = rownames(z_mat),
      group = coef_name,
      estimate = est,
      SE = se,
      df = df,
      t.ratio = tval,
      p.value = pval,
      lower.CL = lwr,
      upper.CL = upr,
      sigma = fit$sigma,
      std_effect_d = est / fit$sigma,
      delta_beta = delta_beta,
      beta_ref = beta_ref,
      row.names = NULL
    )
  })
  
  bind_rows(results)
}

p_to_star <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

plot_heat_delta_simple <- function(res) {
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  # Keep only what we need for plotting
  res_plot <- res %>%
    dplyr::select(outcome, group, estimate, p.value) %>%
    # Drop the "diet4" prefix and trim whitespace
    mutate(group = sub("^diet4", "", group),
           group = trimws(group))  %>%
    # Map to desired cells
    mutate(cell = dplyr::case_when(
      group == "KD+ MIF-" ~ "KD",
      group == "KD- MIF+" ~ "MIF",
      group == "KD+ MIF+" ~ "KD + MIF",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(cell)) %>%
    mutate(cell = factor(cell, levels = c("KD", "MIF", "KD + MIF")))

  p <- ggplot(res_plot, aes(x = cell, y = outcome, fill = estimate)) +
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

# Helper for p-value stars
p_to_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

# Internal plotting function: only KD, MIF, KD + MIF
plot_heat_delta_simple <- function(res) {
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  res_plot <- res %>%
    dplyr::select(outcome, group, estimate, p.value) %>%
    mutate(group = sub("^diet4", "", group),
           group = trimws(group)) %>%
    mutate(cell = dplyr::case_when(
      group == "KD+ MIF-" ~ "KD",
      group == "KD- MIF+" ~ "MIF",
      group == "KD+ MIF+" ~ "KD + MIF",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(cell)) %>%
    mutate(cell = factor(cell, levels = c("KD", "MIF", "KD + MIF")))
  
  ggplot(res_plot, aes(x = cell, y = outcome, fill = estimate)) +
    geom_tile(color = "grey85") +
    scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
      breaks = pretty_breaks(n = 3)
    ) +
    geom_text(aes(label = p_to_star(p.value)), size = 5) +
    theme_minimal(base_size = 10) +
    labs(title = "ΔTreatment", fill = "LE") +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(vjust = 0.8),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    xlab("") + ylab("")
}

# Wrapper: run simplified limma and return one heatmap
wrap_summary_plot_simple <- function(pdat,
                                     outcomes,
                                     outcome_labels = NULL,
                                     moderate = FALSE,
                                     conf.level = 0.95,
                                     covars = NULL,
                                     scale_outcomes = TRUE) {
  library(dplyr)
  library(purrr)
  
  # Optionally scale outcomes for visualization (no effect on inference intent)
  dat_for_fit <- if (scale_outcomes) {
    pdat %>%
      mutate(across(all_of(outcomes), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))
  } else {
    pdat
  }
  
  # Run simplified limma summary (expects diet4 with "KD- MIF-" as reference)
  res <- summarize_outcomes_limma_simple(
    dat = dat_for_fit,
    outcome_cols = outcomes,
    moderate = moderate,
    conf.level = conf.level,
    covars = covars
  )
  
  # Relabel and order outcomes for plotting
  if (is.null(outcome_labels)) outcome_labels <- outcomes
  res$outcome <- factor(
    res$outcome,
    levels = rev(outcomes),
    labels = rev(outcome_labels)
  )
  
  # Produce the heatmap (only KD, MIF, KD + MIF)
  p <- plot_heat_delta_simple(res)
  return(p)
}

summarize_outcomes_limma_simple2 <- function(dat,
                                            outcome_cols,
                                            moderate = FALSE,
                                            conf.level = 0.95,
                                            covars = NULL) {
  
  dat <- dat %>% dplyr::filter(diet4 %in% c("KD- MIF+", "KD+ MIF+")) %>% droplevels()
  
  dat$sample_id <- seq_len(nrow(dat))
  
  # set reference level
  dat$diet4 <- relevel(factor(dat$diet4), "KD- MIF+")
  
  # outcome matrix
  z_mat <- as.matrix(dat[, outcome_cols])
  rownames(z_mat) <- dat$sample_id
  z_mat <- t(z_mat)
  
  # design matrix
  if (!is.null(covars)) {
    covar_formula <- paste(covars, collapse = " + ")
    form <- as.formula(
      paste("~ diet4 +", covar_formula)
    )
  } else {
    form <- ~ diet4
  }
  
  design <- model.matrix(form, data = dat)
  rownames(design) <- dat$sample_id
  
  fit <- lmFit(z_mat, design)
  
  if (moderate) {
    fit <- eBayes(fit)
  }
  
  # Degrees of freedom
  df <- if (moderate) fit$df.total else fit$df.residual
  alpha <- 1 - conf.level
  tcrit <- qt(1 - alpha/2, df)
  
  # Extract coefficients
  # Coefficients are:
  # (Intercept)          -> mean of KD- MIF+
  # diet4KD+ MIF+        -> KD+ MIF+ vs KD- MIF+
  
  beta_fun <- function(M) 2^M / (1 + 2^M)
  ref_samples <- dat$sample_id[
    dat$keto == "KD+" & dat$mifepristone == "MIF-"
  ]
  M_ref <- rowMeans(z_mat[, ref_samples])       # true baseline M
  beta_ref <- beta_fun(M_ref)                   # true baseline β
  
  coef_names <- colnames(fit$coefficients)[-1]  # skip intercept, only want comparisons
  
  results <- lapply(coef_names, function(coef_name) {
    est <- fit$coefficients[, coef_name]
    se  <- fit$stdev.unscaled[, coef_name] * fit$sigma
    tval <- est / se
    pval <- 2 * pt(abs(tval), df, lower.tail = FALSE)
    
    lwr <- est - tcrit * se
    upr <- est + tcrit * se
    
    delta_beta <- beta_fun(M_ref + est) - beta_ref  # biologically interpretable Δβ
    
    data.frame(
      outcome = rownames(z_mat),
      group = coef_name,
      estimate = est,
      SE = se,
      df = df,
      t.ratio = tval,
      p.value = pval,
      lower.CL = lwr,
      upper.CL = upr,
      sigma = fit$sigma,
      std_effect_d = est / fit$sigma,
      delta_beta = delta_beta,
      beta_ref = beta_ref,
      row.names = NULL
    )
  })
  
  bind_rows(results)
}
