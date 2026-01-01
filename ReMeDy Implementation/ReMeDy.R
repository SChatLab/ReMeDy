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