#### Load packages ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'MASS', 'data.table',
              'matrixStats', 'missMethyl', 'evora', 'caret', 'pROC', 'jlst', 'qvalue',
              'tidyverse', 'readr', 'pETM', 'peakRAM', 'onewaytests', 'DMRcate', 'limma',
              'coMethDMR', 'ClustGeo', 'dendextend', 'ENmix', 'ChAMP', 'bumphunter', 'pbapply',
              'doSNOW', 'gtools', 'dplyr', 'ChIPpeakAnno', 'sesameData', 'dglm', 'hglm',
              'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'rngtools', 'STAAR', 'dmrff')

for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}

#### Set Path ####

dir <- 'path_to_directory/Real Data'

#### Load required functions ####

pkgmaker::source_files(paste(dir, "Codes/fit_functions", sep = '/'),'*.R')

# Load Data
load(paste(dir, 'Data', paste0(datnam, '.RData'), sep = '/'))

bval_dat <- data$beta_vals
mvals <- ENmix::B2M(bval_dat)
metadata <- data$metadata
expVar <- 'Exposure'
coVars <- NULL

metadata$Exposure <- factor(metadata$Exposure)

# Get Manifest
manifest_name <- c('HM450.h19.manifest.txt', 'EPIC.hg38.manifest.txt')
if(data_type %in% '450K'){
  manifest <- data.frame(fread(paste(dir, 'Data', manifest_name[1], sep = '/')))
}else{
  manifest <- data.frame(fread(paste(dir, 'Data', manifest_name[2], sep = '/')))
  manifest <- manifest[, c(1:3, ncol(manifest), 4:(ncol(manifest)-1))]
}
manifest <- manifest[manifest$Probe_ID %in% rownames(bval_dat), ]
rownames(manifest) <- manifest$Probe_ID

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
  
}else{
  
  if(Method %in% c('bumphunter', 'DMRcate', 'DMRcateDVMR', 'DMRcateVar')){
    fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars,data_type)", sep = '')
  }else{
    fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars)", sep = '')
  }
  
  fit <- eval(parse(text = fun))
}

#### Save Results ####

sim.scene <- paste(Method, datnam, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')
save(fit, file = paste(dir, "Results", savenam, sep = '/'))
