#### Fit combp To A Dataset ####

fit.combp <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
    
    #### Standard combp pipeline for 2 groups ####
    
    if(is.null(coVars)){
      regData <- metadata[, c(expVar), drop = FALSE]
    }else{
      regData <- metadata[, c(expVar, coVars)]
    }
    formula <- as.formula(paste("~ ", paste(colnames(regData), collapse = "+")))
    design <- model.matrix(object = formula, data = regData)
    
    beta_vals_na <- na.omit(beta_vals)
    manifest.new <- manifest[rownames(manifest) %in% rownames(beta_vals_na), ]
    groups <- factor(metadata$Exposure)
    
    lm.fit <- lmFit(as.matrix(beta_vals_na), design)
    lm.fit <- eBayes(lm.fit)
    
    stats <- data.frame(estimate = lm.fit$coefficients[, 2],
                        se = sqrt(lm.fit$s2.post) * lm.fit$stdev.unscaled[, 2],
                        p = lm.fit$p.value[, 2])
    stats$probe <- rownames(stats)
    
    acombp.annot <- merge(stats, manifest.new, by.x = 'probe', by.y = 'Probe_ID')
    acombp.annot <- acombp.annot[,c('seqnames', 'start', 'end', 'p', 'probe')]
    names(acombp.annot)[1] <- 'chr'
    acombp.annot <- acombp.annot[match(rownames(beta_vals_na), acombp.annot$probe), ]
    
    fit <- combp(data = acombp.annot,
                 dist.cutoff = 200,
                 bin.size = 100,
                 seed = 0.05,
                 region_plot = FALSE,
                 mht_plot = FALSE,
                 nCores = 1,
                 verbose = FALSE)
    fit <- fit[fit$sidak < 0.05 & fit$nprobe > 2, ]
    
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  
  #### Extracting Results ####
  
  res <- fit
  res$p_value <- res$sidak
  res$seqnames <- res$chr
  
  if(!is.data.frame(res)){
    paras <- data.frame(cpgs = row.names(beta_vals_na), adjPval = 1, Pval = 1)
  }else{
    if(nrow(res) == 0){
      paras <- data.frame(cpgs = row.names(beta_vals_na),
                          Pval = 1,
                          adjPval = 1)
    }else{
      temp <- filter_cpgs(res, manifest.new)
      paras <- temp[,c('Probe_ID', 'p_value')]
      colnames(paras) <- c('cpgs', 'adjPval')
      cpgs <- data.frame(cpgs = row.names(beta_vals_na))
      paras <- merge(cpgs, paras, by = 'cpgs', all = TRUE)
      paras$adjPval[is.na(paras$adjPval)] <- 1
      paras <- paras[order(paras$adjPval, decreasing = FALSE),]
      paras$Pval <- paras$adjPval
    }
  }
  finres = list(Method = 'comb-p',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}