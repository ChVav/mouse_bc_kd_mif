
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


# Gives also 3 way interaction, to investigate synergy at the CpG level
summarize_outcomes_limma <- function(dat,
                                     outcome_cols,
                                     moderate = FALSE,
                                     conf.level = 0.95,
                                     covars = NULL) {
  
  library(limma)
  library(dplyr)
  
  # ---------------------------
  # 1. Outcome matrix
  # ---------------------------
  dat$sample_id <- seq_len(nrow(dat))
  
  z_mat <- as.matrix(dat[, outcome_cols])
  rownames(z_mat) <- dat$sample_id
  z_mat <- t(z_mat)   # CpGs x samples
  
  # ---------------------------
  # 2. Design matrix
  # ---------------------------
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
  
  # ---------------------------
  # 3. Exposure contrasts 
  # ---------------------------
  C <- matrix(0, nrow = length(coef_names), ncol = 4)
  rownames(C) <- coef_names
  colnames(C) <- c(
    "KDminus_MIFminus",
    "KDplus_MIFminus",
    "KDminus_MIFplus",
    "KDplus_MIFplus"
  )
  
  # KD- MIF-
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
  
  # ---------------------------
  # 4. Synergy contrast
  # Δ11 - Δ10 - Δ01 + Δ00
  # equals 3-way interaction exactly
  # ---------------------------
  S <- matrix(0, nrow = length(coef_names), ncol = 1)
  rownames(S) <- coef_names
  colnames(S) <- "Synergy"
  
  S["exposureP/D+:ketoKD+:mifepristoneMIF+", "Synergy"] <- 1
  
  # combine all contrasts
  C_all <- cbind(C, S)
  
  fit2 <- contrasts.fit(fit, C_all)
  
  if (moderate) {
    fit2 <- eBayes(fit2)
  }
  
  df <- if (moderate) fit2$df.total else fit2$df.residual
  
  alpha <- 1 - conf.level
  tcrit <- qt(1 - alpha/2, df)
  
  # ---------------------------
  # 5. Reference methylation
  # ---------------------------
  beta_fun <- function(M) 2^M / (1 + 2^M)
  
  ref_samples <- dat$sample_id[
    dat$exposure == "P/D-" &
      dat$keto == "KD-" &
      dat$mifepristone == "MIF-"
  ]
  
  M_ref <- rowMeans(z_mat[, ref_samples, drop = FALSE])
  beta_ref <- beta_fun(M_ref)
  
  # ---------------------------
  # 6. Build output
  # ---------------------------
  results <- lapply(seq_len(ncol(C_all)), function(i) {
    
    est <- fit2$coefficients[, i]
    se  <- fit2$stdev.unscaled[, i] * fit2$sigma
    t   <- est / se
    p   <- 2 * pt(abs(t), df, lower.tail = FALSE)
    
    lwr <- est - tcrit * se
    upr <- est + tcrit * se
    
    nm <- colnames(C_all)[i]
    
    if (nm == "Synergy") {
      keto_lab <- "KDx"
      mif_lab  <- "MIFx"
      delta_beta <- NA_real_
    } else {
      keto_lab <- sub("_.*", "", nm)
      mif_lab  <- sub(".*_", "", nm)
      delta_beta <- beta_fun(M_ref + est) - beta_ref
    }
    
    data.frame(
      outcome = rownames(z_mat),
      contrast = nm,
      keto = keto_lab,
      mifepristone = mif_lab,
      estimate = est,
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
  
  out <- bind_rows(results)
  
  # FDR within each contrast
  out <- out %>%
    group_by(contrast) %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
  
  return(out)
}