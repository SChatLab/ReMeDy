#### Get formulas ####

build_formulas <- function(expression_var = "expr",
                           expVar,
                           coVars = NULL,
                           rand_effect = "(1|Sample)"){
  
  # RHS terms for mean: exposure + covariates + random effect
  rhs_mean_terms <- c(expVar, coVars, rand_effect)
  
  mean.formula <- as.formula(paste(expression_var, "~", paste(rhs_mean_terms, collapse = " + ")))
  
  # Dispersion model: only exposure
  disp.formula <- as.formula(paste("~", expVar))
  
  list(mean = mean.formula,
       disp = disp.formula)
}


#### Region-Based Modeling Function ####

perRegion.mod <- function(reg, metadata, Region, expVar, coVars, rand_effect){
  
  #### Converting data to long format ####
  
  df_long <- reg %>%
    pivot_longer(cols = all_of(colnames(reg)), names_to = "Sample", values_to = "expr")
  merged_df <- merge(df_long, metadata, by = "Sample")
  merged_df <- merged_df[gtools::mixedorder(merged_df$Sample), ]
  
  #### Get formulas ####
  fm <- build_formulas(expression_var = "expr",
                       expVar,
                       coVars,
                       rand_effect)
  
  #### Fitting hglm ####
  
  fit <- tryCatch({
    fit.m <- hglm::hglm2(fm$mean, data = merged_df, family = gaussian(link = identity), 
                         rand.family = gaussian(link = identity), method = "EQL", conv = 1e-6, maxit = 50, 
                         startval = NULL, X.disp = NULL, disp = fm$disp, link.disp = "log", 
                         weights = NULL, fix.disp = NULL, offset = NULL, sparse = TRUE, vcovmat = FALSE, 
                         calc.like = TRUE, RandC = NULL, bigRR = FALSE, verbose = FALSE)
  }, error = function(err){
    fit.m <- NULL
    return(fit.m)
  })
  
  #### Collecting Outputs ####
  
  if(!(is.null(fit) | is.null(fit$SeFe))){
    summary_fit <- summary(fit)
    mu_pval <- summary_fit$FixCoefMat[2,4]
    if(!is.null(fit$SummVC1)){
      disp_pval <- 2 * pnorm(-abs(summary_fit$SummVC1[2,1]/summary_fit$SummVC1[2,2]))
      if(mu_pval == 1){mu_pval = 0.9999}
      if(disp_pval == 1){disp_pval = 0.9999}
      if(mu_pval == 0){mu_pval = 0.0001}
      if(disp_pval == 0){disp_pval = 0.0001}
      mudisp_pval <- STAAR::CCT(c(mu_pval, disp_pval))
    }else{
      disp_pval <- NA
      mudisp_pval <- NA
    }
  }else{
    mu_pval <- NA
    disp_pval <- NA
    mudisp_pval <- NA
  }
  
  #### Effect sizes (log2FC) + SEs for mean and dispersion ####
  
  expPattern <- paste0("^", expVar)
  
  ## Mean model
  if(!is.null(fit) && !is.null(fit$fixef) && !is.null(fit$SeFe)){
    fe_names <- names(fit$fixef)
    idx_mu <- grep(expPattern, fe_names)
    est_mu <- unname(fit$fixef[idx_mu])
    SE_mu <- unname(fit$SeFe[idx_mu])
  }else{
    est_mu <- NA
    SE_mu <- NA
  }
  
  ## Dispersion model
  if(!is.null(fit) && !is.null(fit$SummVC1)){
    rn <- rownames(fit$SummVC1)
    idx_disp <- grep(expPattern, rn)
    beta_disp <- fit$SummVC1[idx_disp, "Estimate"]
    se_beta_disp <- fit$SummVC1[idx_disp, "Std. Error"]
    est_disp <- beta_disp/log(2)
    SE_disp <- se_beta_disp/log(2)
  }else{
    est_disp <- NA
    SE_disp <- NA
  }
  
  #### Summarizing Outputs ####
  
  est.df <- data.frame(Region,
                       Pval_mu = mu_pval,
                       Pval_disp = disp_pval,
                       Pval_mudisp = mudisp_pval,
                       est_mu = est_mu,
                       SE_mu = SE_mu,
                       est_disp = est_disp,
                       SE_disp = SE_disp)
  
  return(est.df)
}