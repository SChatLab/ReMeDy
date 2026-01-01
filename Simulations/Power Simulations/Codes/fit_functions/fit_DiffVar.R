#### Fit DiffVar To A Dataset ####

fit.DiffVar <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard DiffVar pipeline for 2 groups ####
  
  mem <- peakRAM::peakRAM({
    if(is.null(coVars)){
      regData <- metadata[, c(expVar), drop = FALSE]
    }else{
      regData <- metadata[, c(expVar, coVars)]
    }
    formula <- as.formula(paste("~ ", paste(colnames(regData), collapse = "+")))
    design <- model.matrix(object = formula, data = regData)
    beta_vals_na <- na.omit(beta_vals)
    expVar.new <- colnames(design)[grep(expVar, colnames(design))]
    fit <- missMethyl::varFit(beta_vals_na, design = design, coef = c(1, 2), type = 'AD', trend = TRUE, robust = TRUE, weights = NULL)
    res <- missMethyl::topVar(fit, coef = 2, number = length(rownames(beta_vals_na)))
    
    #### Extracting Results ####
    
    paras <- data.frame(cpgs = rownames(res),
                        Pval = res$P.Value,
                        adjPval = res$Adj.P.Value)
    paras <- paras[order(paras$adjPval, decreasing = FALSE),]

  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  cpgs <- data.frame(cpgs = row.names(beta_vals_na))
  paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
  finres = list(Method = 'DiffVar',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}