#### Fit Bartlett-test To A Dataset ####

fit.Bartlett <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard Bartlett-test pipeline for 2 groups ####
  
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
      q1 <- as.numeric(quantile(x, probs = 0.25, na.rm = TRUE))
      q3 <- as.numeric(quantile(x, probs = 0.75, na.rm = TRUE))
      IQR <- q3 - q1
      keep <- ifelse(x > q3 + (IQR * 1.5) | x < q1 - (IQR * 1.5), FALSE, TRUE)
      g.new <- g[keep]
      g1.new <- g.new[grep('g1', names(g.new))]
      g2.new <- g.new[grep('g2', names(g.new))]
      g12.new <- c(g1.new, g2.new)
      grpvar <- as.factor(c(rep(1, length(g1.new)), rep(2, length(g2.new))))
      fit <- tryCatch({
        fit.b <- stats::bartlett.test(x = g12.new, g = grpvar)
      }, error = function(err){
        fit.b <- NULL
        return(fit.b)
      })
      if(is.null(fit)){
        resTab <- data.frame(statistic = NA, Pval = NA)
      }else{
        resTab <- data.frame(statistic = as.numeric(fit$statistic), Pval = fit$p.value)
      }
      return(resTab)
    }
    
    #### Per-cpg Fitting ####
    
    beta_vals_na <- na.omit(beta_vals)
    res <- apply(beta_vals_na, 1, function(x) spe(x))
    res <- do.call('rbind', res)
    
    #### Extracting Results ####
    
    paras <- data.frame(cpgs = rownames(res),
                        Coef = res$statistic,
                        Pval = res$Pval)
    paras$adjPval <- p.adjust(paras$Pval, method = 'BH')
    paras <- paras[order(paras$adjPval, decreasing = FALSE),]

  })
  stop.time <- Sys.time()
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  cpgs <- data.frame(cpgs = row.names(beta_vals_na))
  paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
  finres = list(Method = 'Bartlett-Test',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}
