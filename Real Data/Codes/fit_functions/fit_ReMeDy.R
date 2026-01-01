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

#### Main Function ####

fit.ReMeDy <- function(region_list, 
                       metadata, 
                       expVar, 
                       coVars,
                       rand_effect,
                       parallel = FALSE,
                       verbose = FALSE,
                       numCores = 4){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard hglm pipeline for 2 groups ####
  
  mem <- peakRAM::peakRAM({
    
    #### Per-Region Fitting ####
    if(parallel){
      if(verbose){
        cl <- parallel::makeCluster(numCores)
        doSNOW::registerDoSNOW(cl)
        packages <- c('hglm', 'tidyr')
        exports <- c('perRegion.mod', 'metadata', 'expVar', 'coVars', 'build_formulas', 'rand_effect')
        pb <- txtProgressBar(max = length(region_list), style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        res <- foreach(j = 1:length(region_list), .combine = rbind, .packages = packages, .options.snow = opts, .export = exports) %dopar% {
          reg <- region_list[[j]]
          tmpfit <- perRegion.mod(reg,
                                  metadata,
                                  Region = paste0('Region', j),
                                  expVar, coVars, rand_effect)
          return(tmpfit)
        }
        close(pb)
        stopCluster(cl)
      }else{
        cl <- parallel::makeCluster(numCores)
        doParallel::registerDoParallel(cl)
        packages <- c('hglm', 'tidyr')
        exports <- c('perRegion.mod', 'metadata', 'expVar', 'coVars', 'build_formulas', 'rand_effect')
        res <- foreach(j = 1:length(region_list), .combine = rbind, .packages = packages, .export = exports) %dopar% {
          reg <- region_list[[j]]
          tmpfit <- perRegion.mod(reg,
                                  metadata,
                                  Region = paste0('Region', j),
                                  expVar, coVars, rand_effect)
          return(tmpfit)
        }
        stopCluster(cl)
      }
    }else{
      res <- pblapply(seq_along(region_list), function(j){
        x <- region_list[[j]]
        perRegion.mod(x, metadata, Region = paste0("Region", j), expVar, coVars, rand_effect)
      })
      res <- do.call('rbind', res)
    }
    
    #### Extracting Results ####
    
    # Adjust p-values (BH method)
    res$adjPval_mu <- p.adjust(res$Pval_mu, method = 'BH')
    res$adjPval_disp <- p.adjust(res$Pval_disp, method = 'BH')
    res$adjPval_mudisp <- p.adjust(res$Pval_mudisp, method = 'BH')
    
    # Create an empty list to hold CpG-level dataframes for each region
    df_list <- list()
    
    # Loop over each region in list
    for (region_name in names(region_list)){
      # Extract CpG identifiers (rows) for this region
      cpgs <- rownames(region_list[[region_name]])
      
      # Extract the corresponding p-value results for this region
      pval_info <- res[res$Region == region_name, ]
      
      # Build a dataframe for CpGs in this region
      df_temp <- data.frame(cpgs = cpgs,
                            Region = region_name,
                            Pval_mu = pval_info$Pval_mu,
                            adjPval_mu = pval_info$adjPval_mu,
                            est_mu = pval_info$est_mu,
                            SE_mu = pval_info$SE_mu,
                            Pval_disp = pval_info$Pval_disp,
                            adjPval_disp = pval_info$adjPval_disp,
                            est_disp = pval_info$est_disp,
                            SE_disp = pval_info$SE_disp,
                            Pval_mudisp = pval_info$Pval_mudisp,
                            adjPval_mudisp = pval_info$adjPval_mudisp)
      
      # Store in list
      df_list[[region_name]] <- df_temp
    }
    
    # Combine all region-specific dataframes into one CpG-level dataframe
    df_cpg_info <- do.call(rbind, df_list)
    rownames(df_cpg_info) <- NULL
    
  })
  stop.time <- Sys.time()
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'ReMeDy',
                res = df_cpg_info,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}