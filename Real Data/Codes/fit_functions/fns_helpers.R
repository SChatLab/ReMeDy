#### Obtain CpG-regions based on maximum base pair threshold ####

get.Genoregs <- function(manifest, 
                         maxGap = 1000, 
                         minCpgs = 3, 
                         ceiling = c('equal', 'gr.equal'), 
                         intergenic = TRUE){
  if(intergenic == FALSE){
    manifest <- filter(manifest, manifest$UCSC_RefGene_Group != "")
  }
  manifest <- split(manifest, manifest$seqnames)
  manifest <- lapply(manifest, function(x) x[order(x$probeTarget), ])
  manifest <- do.call('rbind', manifest)
  manifest$bin <- bumphunter::clusterMaker(manifest$seqnames, manifest$probeTarget, assumeSorted = TRUE, maxGap = maxGap)
  manifest <- split(manifest, manifest$bin)
  if(ceiling %in% 'equal'){
    minCpgs.regions <- Filter(function(x) nrow(x) == minCpgs, manifest)
  }else{
    minCpgs.regions <- Filter(function(x) nrow(x) >= minCpgs, manifest)
  }
  minCpgs.names <- lapply(minCpgs.regions, function(x) paste(unique(x$seqnames), min(x$probeTarget), max(x$probeTarget), sep = ':'))
  names(minCpgs.regions) <- minCpgs.names
  return(minCpgs.regions)
}

#### Function to calculate evaluation metrics ####

perfMetrics <- function(fit, dm.cpgs, simulation, alpha, effect.type, all.metrics = TRUE){
  Method <- fit$Method
  resi <- fit$res
  if(Method %in% c('ReMeDy')){
    if(effect.type %in% 'mean'){
      resi$adjPval <- resi$adjPval_mu
    }else if(effect.type %in% 'var'){
      resi$adjPval <- resi$adjPval_disp
    }else{
      resi$adjPval <- resi$adjPval_mudisp
    }
  }
  resi[is.na(resi[, "adjPval"]), "adjPval"] <- 1
  adjPval <- resi$adjPval
  names(adjPval) <- resi$cpgs
  response <- rep(FALSE, length(adjPval))
  names(response) <- resi$cpgs
  response[names(response) %in% dm.cpgs] <- TRUE
  response <- factor(response, levels = c("TRUE", "FALSE"))
  predictor <- factor(adjPval <= 0.05, levels = c("TRUE","FALSE"))
  xtab <- table(predictor, response)
  cm <- caret::confusionMatrix(xtab)
  
  #### Getting Metrics ####
  
  if(all.metrics){
    Sensitivity <- format(cm$byClass['Sensitivity'], scientific = FALSE)
    Sensitivity <- as.numeric(substr(as.character(Sensitivity), 1, 4))
    Specificity <- format(cm$byClass['Specificity'], scientific = FALSE)
    Specificity <- as.numeric(substr(as.character(Specificity), 1, 4))
    FDR <- format(1 - cm$byClass['Precision'], scientific = FALSE)
    FDR <- as.numeric(substr(as.character(FDR), 1, 4))
    F1 <- format(cm$byClass['F1'], scientific = FALSE)
    F1 <- as.numeric(substr(as.character(F1), 1, 4))
    MCC <- format(round(mltools::mcc(preds = predictor, actuals = response), 2), scientific = FALSE)
    MCC <- as.numeric(substr(as.character(MCC), 1, 4))
    
    #### ROC ####
    
    response <- rep(FALSE, length(adjPval))
    names(response) <- resi$cpgs
    response[names(response) %in% dm.cpgs] <- TRUE
    response <- factor(response)
    predictor <- adjPval
    proc_out1 <- pROC::roc(response = response, predictor = predictor, direction = ">", quiet = TRUE)
    ROC <- as.numeric(proc_out1$auc)
    
    #### pAUROC ####
    
    proc_out2 <- roc(response = response, predictor = predictor,
                     direction = ">",quiet = TRUE, partial.auc = c(1, 0.85), 
                     partial.auc.correct = TRUE)
    pAUROC <- as.numeric(proc_out2$auc)
    
    #### Summarizing Metrics ####
    
    Measures <- data.frame(Method = Method,
                           Power = Sensitivity,
                           FPR = 1 - Specificity,
                           FDR = FDR,
                           F1 = F1,
                           MCC = MCC,
                           ROC = ROC,
                           pAUROC = pAUROC,
                           alpha = alpha,
                           replicate = simulation,
                           Time = fit$Time.min,
                           PeakMemory_GB = fit$peak.memory.gb)
  }else{
    Sensitivity <- round(as.numeric(cm$byClass['Sensitivity']), 2)
    Specificity <- round(as.numeric(cm$byClass['Specificity']), 2)
    FDR <- round(as.numeric(1 - cm$byClass['Precision']), 2)
    F1 <- round(as.numeric(1 - cm$byClass['F1']), 2)
    MCC <- round(mltools::mcc(preds = predictor, actuals = response), 2)
    
    #### Summarizing Metrics ####
    
    Measures <- data.frame(Method = Method,
                           Power = Sensitivity,
                           FPR = 1 - Specificity,
                           FDR = FDR,
                           F1 = F1,
                           MCC = MCC,
                           alpha = alpha,
                           replicate = simulation,
                           Time = fit$Time.min,
                           PeakMemory_GB = fit$peak.memory.gb)
  }
  return(Measures)
}

#### Function to Benchmark Models for Null Analysis ####

permMetrics <- function(fit, simulation, alpha, effect.type){
  Method <- fit$Method
  resi <- fit$res
  if(Method %in% c('ReMeDy')){
    if(effect.type %in% 'mean'){
      resi$Pval <- resi$Pval_mu
      resi$adjPval <- resi$adjPval_mu
    }else if(effect.type %in% 'var'){
      resi$Pval <- resi$Pval_disp
      resi$adjPval <- resi$adjPval_disp
    }else{
      resi$Pval <- resi$Pval_mudisp
      resi$adjPval <- resi$adjPval_mudisp
    }
  }
  resi[is.na(resi[, "Pval"]), "Pval"] <- 1
  resi[is.na(resi[, "adjPval"]), "adjPval"] <- 1
  DMcpgs.p <- resi[resi$Pval < alpha, ]
  DMcpgs.fdr <- resi[resi$adjPval < alpha, ]
  if(dim(DMcpgs.p)[1] == 0 & dim(DMcpgs.fdr)[1] == 0){
    permfit <- data.frame(Method = Method,
                          nfeatures = nrow(resi),
                          nDMcpgs.p = 0,
                          nDMcpgs.fdr = 0,
                          pDMcpgs.p = 0,
                          pDMcpgs.fdr = 0,
                          simulation = simulation,
                          Time = fit$Time.min,
                          PeakMemory_GB = fit$peak.memory.gb)
  }else{
    permfit <- data.frame(Method = Method,
                          nfeatures = nrow(resi),
                          nDMcpgs.p = nrow(DMcpgs.p),
                          nDMcpgs.fdr = nrow(DMcpgs.fdr),
                          pDMcpgs.p = nrow(DMcpgs.p)/nrow(resi),
                          pDMcpgs.fdr = nrow(DMcpgs.fdr)/nrow(resi),
                          simulation = simulation,
                          Time = fit$Time.min,
                          PeakMemory_GB = fit$peak.memory.gb)
  }
  return(permfit)
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

filter_cpgs <- function(df1, df2){
  filtered_rows <- do.call(rbind, lapply(1:nrow(df1), function(i) {
    df2_filtered <- df2 %>%
      filter(seqnames == df1$seqnames[i],
             probeTarget >= df1$start[i],
             probeTarget <= df1$end[i]) %>%
      mutate(
        p_value = df1$p_value[i],
        region = paste0("region_", i)
      )
    return(df2_filtered)
  }))
  
  return(filtered_rows)
}

#### Obtain correlation of CpG-regions based on maximum base pair threshold ####

get.GenoregsCorr <- function(minCpgs.regions, 
                             bval_dat,
                             method = 'spearman',
                             minCorr = 0.5,
                             maxCorr = 0.8,
                             numcores = 8){
  cl <- parallel::makeCluster(numcores)
  doParallel::registerDoParallel(cl)
  exports <- c('method', "bval_dat", 'minCorr', 'maxCorr')
  regsCorr <- foreach(l = 1:length(minCpgs.regions), .export = exports) %dopar%{
    x <- minCpgs.regions[[l]]
    bval.tmp <- bval_dat[rownames(bval_dat) %in% x$Probe_ID, ]
    cor.tmp <- cor(t(bval.tmp), method = method)
    cor.tmp <- abs(cor.tmp[upper.tri(cor.tmp, diag = FALSE)])
    cor.tmp1 <- median(cor.tmp)
    Truth <- ifelse(cor.tmp1 >= minCorr & cor.tmp1 <= maxCorr, 'signal', 'noise')
    x$cor.stats <- cor.tmp1
    x$Truth <- Truth
    return(x)
  }
  stopCluster(cl)
  names(regsCorr) <- names(minCpgs.regions)
  regsCorr.sub <- Filter(function(x) unique(x$Truth) %in% 'signal' , regsCorr)
  return(regsCorr.sub)
}
