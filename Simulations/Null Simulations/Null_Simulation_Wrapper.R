#### Load Generic Libraries ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'MASS', 'data.table',
              'matrixStats', 'missMethyl', 'evora', 'caret', 'pROC', 'jlst', 'qvalue', 'car',
              'tidyverse', 'readr', 'pETM', 'peakRAM', 'onewaytests', 'DMRcate', 'limma',
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
dir <- 'path_to_directory/Simulations/Null Simulations'

#### Load functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Load Data and manifest files ####

if(datnam == 'GSE80970_450K'){
  load(paste(dir, "Data", "GSE80970_450K.RData", sep = '/'))
}else if (datnam == 'Preeclampsia_450K'){
  load(paste(dir, "Data", "Preeclampsia_450K.RData", sep = '/'))
}

bval_dat <- data$beta_vals
meta <- data$metadata
manifest <- data.frame(fread(paste(dir, "Data", 'HM450.h19.manifest.txt', sep = '/')))
rownames(manifest) <- manifest$Probe_ID

nsamples <- nrow(meta)

if(nrow(manifest) != nrow(bval_dat)){
  manifest <- manifest[rownames(manifest) %in% rownames(bval_dat), ]
}

meta$Exposure <- factor(meta$Exposure)

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
             'organize.island.repeated', 'meta', 'bval_dat')
restmp <- foreach(r = 1:nreps, .combine = rbind, .packages = packages, .export = exports) %dopar% {
  seedR <- 1233 + r
  set.seed(seedR)
  
  # Remove any effect by permutating columns
  beta_vals.permuted <- bval_dat[, sample(ncol(bval_dat))]
  manifest.new <- manifest[rownames(manifest) %in% rownames(beta_vals.permuted), ]
  
  # Get regions
  minCpgs.regions <- get.Genoregs(manifest.new, 
                                  maxGap = 200, 
                                  minCpgs = 3, 
                                  ceiling = 'equal', 
                                  intergenic = FALSE)
  
  # Get Null Data
  samp.regs <- sample(minCpgs.regions, 10000)
  reg.list <- do.call(rbind, samp.regs)
  beta_vals <- beta_vals.permuted[reg.list$Probe_ID,]
  names(beta_vals) <- paste0('Sample', 1:ncol(beta_vals))
  metadata <- data.frame(Sample = paste0('Sample', 1:ncol(beta_vals)),
                         Exposure = rep(c(1,0), times = c(ncol(beta_vals)/2, ncol(beta_vals)/2)))
  mvals <- ENmix::B2M(beta_vals)
  expVar = 'Exposure'
  coVars = NULL
  
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
    
    performance <- permMetrics(fit,
                               simulation = r, 
                               alpha = 0.05,
                               effect.type = 'both')
  }else{
    
    if(Method %in% c('bumphunter', 'DMRcate', 'DMRcateDVMR', 'DMRcateVar')){
      data_type <- '450K'
      fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars,data_type)", sep = '')
    }else{
      fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars)", sep = '')
    }
    
    fit <- eval(parse(text = fun))
    
    #### Collect performance metrics ####
    
    performance <- permMetrics(fit,
                               simulation = r, 
                               alpha = 0.05,
                               effect.type = 'Nope')
  }
  return(performance)
}
stopCluster(cl)

#### Save Results ####

sim.scene <- paste(Method, nsamples, 'Type1', sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')
save(restmp, file = paste(dir, "Results/Null", savenam, sep = '/'))