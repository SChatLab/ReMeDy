#### Fit JLSsc To A Dataset ####

fit.JLSsc <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard JLSsc pipeline for 2 groups ####
  
  mem <- peakRAM::peakRAM({
    conds <- unique(metadata[, expVar])
    if(is.data.frame(conds)){
      conds <- conds$Exposure
    }
    group0indexes <- which(metadata[, expVar] == conds[1])
    group1indexes <- which(metadata[, expVar] == conds[2])
    spe <- function(x){
      g1 <- x[group0indexes]
      g2 <- x[group1indexes]
      names(g1) <- paste('g1_', 1:length(g1), sep = '')
      names(g2) <- paste('g2_', 1:length(g2), sep = '')
      g <- c(g1, g2)
      grpvar <- as.factor(c(rep(1, length(g1)), rep(2, length(g2))))
      fit <- jlst::jlssc(y = x, x = grpvar, type = 1, covar = NULL)
      resTab <- data.frame(Pval = fit$P)
      return(resTab)
    }
    
    #### Per-gene Fitting ####
    
    beta_vals_na <- na.omit(beta_vals)
    res <- apply(beta_vals_na, 1, function(x) spe(x))
    res <- do.call('rbind', res)
    
    #### Extracting Results ####
    
    paras <- data.frame(cpgs = rownames(res),
                        Pval = res$Pval)
    paras$adjPval <- p.adjust(paras$Pval, method = 'BH')
    paras <- paras[order(paras$adjPval, decreasing = FALSE),]

  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  cpgs <- data.frame(cpgs = row.names(beta_vals_na))
  paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
  finres = list(Method = 'JLSsc',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}