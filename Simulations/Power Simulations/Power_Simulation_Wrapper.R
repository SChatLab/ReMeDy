#### Load Generic Libraries ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'MASS', 'data.table', 'lmerTest',
              'matrixStats', 'missMethyl', 'evora', 'caret', 'pROC', 'jlst', 'qvalue', 'car', 'gee',
              'tidyverse', 'readr', 'pETM', 'peakRAM', 'onewaytests', 'DMRcate', 'limma', 'tidyr',
              'coMethDMR', 'ClustGeo', 'dendextend', 'ENmix', 'ChAMP', 'bumphunter', 'pbapply',
              'doSNOW', 'gtools', 'dplyr', 'ChIPpeakAnno', 'sesameData', 'dglm', 'hglm', 'robust',
              'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'rngtools', 'STAAR', 'dmrff')

for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}

#### Set Path ####
dir <- 'path_to_directory/Simulations/Power Simulations'

#### Load functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Load Data and manifest files ####

load(paste(dir, "Data", "Preeclampsia_450K.RData", sep = '/'))
bval_dat <- data$beta_vals
meta <- data$metadata
manifest <- data.frame(fread(paste(dir, "Data", 'HM450.h19.manifest.txt', sep = '/')))
rownames(manifest) <- manifest$Probe_ID

if(nrow(manifest) != nrow(bval_dat)){
  manifest <- manifest[rownames(manifest) %in% rownames(bval_dat), ]
}

#### Get Regions with required number of CpGs ####

minCpgs.regions <- get.Genoregs(manifest, 
                                maxGap = 200, 
                                minCpgs = ncpgs.regs, 
                                ceiling = 'equal', 
                                intergenic = FALSE)

reg.list <- do.call(rbind, minCpgs.regions)

#### Generate Simulated data with 100 replicates (in-parallel) ####

ncores <- 10
nreps <- 100
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
exports <- c('acf.table', 'combp', 'doDV', 'doTT', 'filter_cpgs', 'fit.bumphunter', 'fit.combp', 'fit.DGLM', 'fit.SacRobScore',
             'fit.DiffVar', 'fit.DMRcate', 'fit.DMRcateDVMR', 'fit.DMRcateVar', 'fit.dmrff', 'fit.iEVORA', 'fit.SacRobBeta', 
             'trigger.Aclust2', 'fit.ipdmr', 'fit.JLSp', 'fit.JLSsc', 'fit.LR', 'fit.LRRobust', 'fit.pETM', 'fit.probelasso', 
             'fit.RobMethBeta', 'perRegion.mod.beta', 'GEE.clusters', 'fit.RobMeth', 'fit.RobMethBH', 'fit.SacRob', 'fit.SacRobBH', 
             'get.acf', 'get.cpgList', 'fit.Levene', 'score_test_hglm', 'get.Genoregs', 'get.GenoregsCorr','iEVORA', 'ipdmr', 'Dbp.merge',
             'perfMetrics', 'permMetrics', 'perRegion.mod', 'fit.RobMethScore', 'regplot', 'Sacoma', 'simArray', 'fit.aclust2gee',
             'unregister_dopar', 'manifest', 'dir', 'fit.Bartlett', 'fit.BrownForsythe', 'perRegion.mod.score', 'find_cluster_list',
             'calc.dist.d.neighbor', 'update.clust.indicator', 'acluster', 'calc.dist.clusters', 'calc.dist.clust.point',
             'organize.island.repeated', 'perRegion.mod.lmm', 'fit.lmm', 'perRegion.mod.gee', 'fit.gee', 'fit.randcoef',
             'perRegion.mod.randcoef')
restmp <- foreach(r = 1:nreps, .combine = rbind, .packages = packages, .export = exports) %dopar% {
  seedR <- 1233 + r
  
  simData <- simArray(nregions = nregions,
                      cpgs.per.region = ncpgs.regs,
                      nsamples = nsamples,
                      var.effect = var.effect,
                      mean.effect = mean.effect,
                      grp.split.pc = grp.split.pc,
                      pDM = pDM,
                      effect.type = effect.type,
                      outpc = outpc,
                      seed = seedR,
                      hi_corr_prop = 0.5,
                      hi_corr_range = c(0.6, 0.95),
                      lo_corr_range = c(0.05, 0.2),
                      corr.type = corr.type)
  
  #### Gathering Simulated Data ####
  
  mvals <- simData$mval
  rownames(mvals) <- reg.list$Probe_ID[1:nrow(mvals)]
  metadata <- simData$metadata
  cpgs_siminfo <- simData$CpG.Info
  cpgs_siminfo$CpGs <- reg.list$Probe_ID[1:nrow(mvals)]
  dm.cpgs <- cpgs_siminfo$CpGs[cpgs_siminfo$Truth != 'None']
  expVar <- 'Exposure'
  coVars <- NULL
  
  #### Fit models ####
  if(Method %in% c('ReMeDy')){
    
    sacoma.fit <- Sacoma(bval.dnam = mvals,
                         data_type = "450k",
                         manifest.dir = 'path_to_manifest',
                         minCpGs = ncpgs.regs,
                         minCpGs.thres.type = "equal",
                         maxGap = 200,
                         rthresold = 0.5, 
                         method = "spearman",
                         ncores = 1)
    
    # Get significant regions
    mvals.sub <- mvals[sacoma.fit$Probe_ID[sacoma.fit$prediction == 'signal'], ]
    mvals.sub$CpGs <- rownames(mvals.sub)
    merged <- merge(sacoma.fit, mvals.sub, by.x = 'Probe_ID', by.y = 'CpGs')
    merged <- merged[mixedorder(merged$Region),]
    region_list <- split(merged, merged$Region)
    
    region_list_clean <- lapply(region_list, function(x) {
      rownames(x) <- x$Probe_ID
      x <- x[, !(names(x) %in% c("CpGs", "Probe_ID", "CHR", "POS", "bin", "Region", "prediction", "Time.min", "peak.memory.gb"))]
      return(x)
    })
    
    names(region_list_clean) <- paste0('Region', 1:length(region_list_clean))
    
    fit <- fit.ReMeDy(region_list_clean, 
                      metadata, 
                      expVar, 
                      coVars,
                      rand_effect = "(1 | Sample)",
                      parallel = FALSE,
                      verbose = FALSE,
                      numCores = 1)
    
    #### Collect performance metrics ####
    
    if(pDM == 0){
      performance <- permMetrics(fit,
                                 simulation = r, 
                                 alpha = 0.05,
                                 effect.type)
    }else{
      performance <- perfMetrics(fit, 
                                 dm.cpgs, 
                                 simulation = r, 
                                 alpha = 0.05,
                                 effect.type,
                                 all.metrics = TRUE)
    }
  }else{
    
    if(Method %in% c('bumphunter', 'DMRcate', 'DMRcateDVMR', 'DMRcateVar')){
      data_type <- '450K'
      fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars,data_type)", sep = '')
    }else{
      fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars)", sep = '')
    }
    
    fit <- eval(parse(text = fun))
    
    #### Collect performance metrics ####
    
    if(pDM == 0){
      performance <- permMetrics(fit,
                                 simulation = r, 
                                 alpha = 0.05,
                                 effect.type = 'Nope')
    }else{
      performance <- perfMetrics(fit, 
                                 dm.cpgs, 
                                 simulation = r, 
                                 alpha = 0.05,
                                 effect.type = 'Nope',
                                 all.metrics = TRUE)
    }
  }
  return(performance)
}
stopCluster(cl)

#### Save Results ####

sim.scene <- paste(Method, nsamples, mean.effect, var.effect, effect.type, nregions, ncpgs.regs,
                   grp.split.pc, pDM, outpc, corr.type, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')


if(effect.type == 'mean'){
  save(restmp, file = paste(dir, "Results/Power/DMR", savenam, sep = '/'))
}else if(effect.type == 'var'){
  save(restmp, file = paste(dir, "Results/Power/VMR", savenam, sep = '/'))
}else{
  save(restmp, file = paste(dir, "Results/Power/DVMR", savenam, sep = '/'))
}
