#### iEVORA Auxiliary Function 1 ####

doDV <- function(tmp.v,pheno.v){
  co.idx <- which(pheno.v==0)
  ca.idx <- which(pheno.v==1)
  bt.o <- stats::bartlett.test(x=tmp.v,g=pheno.v)
  pv <- bt.o$p.value
  logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]))
  avCA <- mean(tmp.v[ca.idx])
  avCO <- mean(tmp.v[co.idx])
  out.v <- c(logR,pv,avCA,avCO)
  names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0")
  return(out.v)
}

#### iEVORA Auxiliary Function 2 ####

doTT <- function(tmp.v,pheno.v){
  tt.o <- stats::t.test(tmp.v ~ pheno.v)
  out.v <- c(-tt.o$stat,tt.o$p.val)
  names(out.v) <- c("t","P")
  return(out.v)
}

#### iEVORA Main Function ####

iEVORA <- function(data.m = beta_vals, pheno.v, thDV = 0.05,thDM = 0.05){
  statDVC.m <- t(apply(data.m, 1, doDV, pheno.v))
  print("Estimated DV statistics")
  qvDVC.v <- qvalue(statDVC.m[,2])$qval
  dvc.idx <- which(qvDVC.v < thDV)
  nDVC <- length(dvc.idx)
  if( nDVC > 0 ){
    statDMC.m <- t(apply(data.m[dvc.idx,],1,doTT,pheno.v))
    print("Preparing output")
    tmp.s <- sort(statDMC.m[,2],decreasing=FALSE,index.return=TRUE)
    pvDMC.v <- tmp.s$x
    ntop <- length(which(pvDMC.v < thDM))
    if(ntop > 0){
      if(ntop == 1){
        topDVMC.m <- c(statDMC.m[tmp.s$ix[1:ntop],],
                       statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],
                       qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]])
        topDVMC.m <- data.frame(as.list(topDVMC.m))
      }else{
        topDVMC.m <- cbind(statDMC.m[tmp.s$ix[1:ntop],],
                           statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],
                           qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]])
      }
      colnames(topDVMC.m) <- c("t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)")
      rownames(topDVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]]
    }
    else {
      topDVMC.m <- data.frame(t = rep(0, nrow(statDVC.m)),
                              PTT = rep(0, nrow(statDVC.m)),
                              Av1 = statDVC.m[,3],
                              Av0 = statDVC.m[,4],
                              V1overV0 = statDVC.m[,1],
                              PBT = statDVC.m[,2],
                              qBT = qvDVC.v)
      names(topDVMC.m) <- c("t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)")
    }
  }
  else {
    topDVMC.m <- data.frame(t = rep(0, nrow(statDVC.m)),
                            PTT = rep(0, nrow(statDVC.m)),
                            Av1 = statDVC.m[,3],
                            Av0 = statDVC.m[,4],
                            V1overV0 = statDVC.m[,1],
                            PBT = statDVC.m[,2],
                            qBT = rep(1, nrow(statDVC.m)))
    names(topDVMC.m) <- c("t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)")
  }
  return(topDVMC.m)
}

#### Fit iEVORA To A Dataset ####

fit.iEVORA <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard iEVORA pipeline for 2 groups ####
  
  mem <- peakRAM::peakRAM({
    conds <- unique(metadata[, expVar])
    if(is.data.frame(conds)){
      conds <- conds$Exposure
    }
    group0indexes <- which(metadata[, expVar] == conds[1])
    group1indexes <- which(metadata[, expVar] == conds[2])
    grpvar <-c(rep(0, length(group0indexes)), rep(1, length(group1indexes)))
    beta_vals_na <- na.omit(beta_vals)
    res <- as.data.frame(iEVORA(data.m = beta_vals_na, pheno.v = grpvar, thDV = 0.05, thDM = 0.05))
    
    #### Extracting Results ####
    
    paras <- data.frame(cpgs = rownames(res), Pval = res$`P(BT)`, adjPval_BT = res$`q(BT)`)
    paras$adjPval <- p.adjust(paras$Pval, method = 'BH')
    cpgs <- data.frame(cpgs = row.names(beta_vals_na))
    paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
    paras$adjPval <- ifelse(is.na(paras$adjPval), 1, paras$adjPval)
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'iEVORA',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}