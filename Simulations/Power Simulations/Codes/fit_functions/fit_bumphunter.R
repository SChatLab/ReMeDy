#### Fit bumphunter To A Dataset ####

fit.bumphunter <- function(beta_vals, metadata, expVar, coVars, data_type){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
    #### Standard bumphunter pipeline for 2 groups ####
    
    beta_vals_na <- na.omit(beta_vals)
    beta_vals_mat <- as.matrix(beta_vals_na)
    metadata$Exposure <- as.factor(metadata$Exposure)
    manifest.new <- manifest[rownames(manifest) %in% rownames(beta_vals_na), ]
    
    fit <- tryCatch({
      fitt <- ChAMP::champ.DMR(beta = ENmix::M2B(beta_vals_mat), 
                               pheno = metadata$Exposure, 
                               arraytype = data_type, 
                               method = "Bumphunter",
                               minProbes = 2, 
                               adjPvalDmr = 0.05, 
                               maxGap = 200,
                               pickCutoff = TRUE,
                               smooth = FALSE,
                               smoothFunction = loessByCluster,
                               useWeights = FALSE,
                               B = 10,
                               nullMethod = 'permutation',
                               cores = 1)
    },
    error = function(e){
      return(NULL)
    })
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  
  #### Extracting Results ####
  if(is.null(fit)){
    paras <- data.frame(cpgs = row.names(beta_vals_na), Pval = 1, adjPval = 1)
  }else{
    res <- fit$BumphunterDMR
    res$p_value <-res$p.valueArea
    
    if(sum(res$p_value < 0.05) != 0){
      temp <- filter_cpgs(res, manifest.new)
      paras <- temp[,c('Probe_ID', 'p_value')]
      colnames(paras) <- c('cpgs', 'adjPval')
      cpgs <- data.frame(cpgs = row.names(beta_vals_na))
      paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
      paras$adjPval[is.na(paras$adjPval)] <- 1
      paras <- paras[order(paras$adjPval, decreasing = FALSE),]
      paras$Pval <- paras$adjPval
    }else{
      paras <- data.frame(cpgs = row.names(beta_vals_na), Pval = 1, adjPval = 1)
    }
  }
  
  finres = list(Method = 'bumphunter',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}