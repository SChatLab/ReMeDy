#### Fit DMRcate To A Dataset ####

fit.DMRcate <- function(beta_vals, metadata, expVar, coVars, data_type){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
    #### Standard DMRcate pipeline for 2 groups ####
    
    if(is.null(coVars)){
      regData <- metadata[, c(expVar), drop = FALSE]
    }else{
      regData <- metadata[, c(expVar, coVars)]
    }
    formula <- as.formula(paste("~ ", paste(colnames(regData), collapse = "+")))
    design <- model.matrix(object = formula, data = regData)
    beta_vals_na <- na.omit(beta_vals)
    manifest.new <- manifest[rownames(manifest) %in% rownames(beta_vals_na), ]
    
    myAnnotation <- DMRcate::cpg.annotate(object = as.matrix(beta_vals_na),
                                          datatype = "array",
                                          what = "M",
                                          analysis.type = "differential",
                                          design = design,
                                          coef = 2,
                                          arraytype = data_type,
                                          contrasts = FALSE, 
                                          cont.matrix = NULL, 
                                          fdr = 0.05, 
                                          varFitcoef = NULL, 
                                          topVarcoef = NULL)
    
    if(sum(myAnnotation@ranges@elementMetadata@listData$is.sig) != 0){
      dmr.object <- DMRcate::dmrcate(myAnnotation, lambda = 1000, C = 2, min.cpgs = 3)
      
      if(length(dmr.object@coord) != 0){
        res <- data.frame(DMRcate::extractRanges(dmr.object, genome = "hg19"))
        res$p_value <- res$Stouffer
        
        #### Extracting Results ####
        
        if(sum(res$Stouffer < 0.05) != 0){
          temp <- filter_cpgs(res, manifest.new)
          paras <- temp[,c('Probe_ID', 'p_value')]
          colnames(paras) <- c('cpgs', 'adjPval')
          cpgs <- data.frame(cpgs = row.names(beta_vals_na))
          paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
          paras$adjPval[is.na(paras$adjPval)] <- 1
          paras <- paras[order(paras$adjPval, decreasing = FALSE),]
          paras$Pval <- paras$adjPval
        }else{
          paras <- data.frame(cpgs = row.names(beta_vals_na), adjPval = 1, Pval = 1)
        }
      }else{
        paras <- data.frame(cpgs = row.names(beta_vals_na), adjPval = 1, Pval = 1)
      }
    }else{
      paras <- data.frame(cpgs = row.names(beta_vals_na), adjPval = 1, Pval = 1)
    }
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  
  finres = list(Method = 'DMRcate',
                res = paras, 
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}