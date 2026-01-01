#### Function to simulate Methylation Microarray Data ####

simArray <- function(nregions = 10000, # number of regions
                     cpgs.per.region = 3, # number of cpgs per region
                     nsamples = 50, # number of samples
                     var.effect = 0.96, # variance effect
                     mean.effect = 2, # mean effect
                     grp.split.pc = 70, # group proportion
                     pDM = 0.1, # proportion of regions with methylation dysregulation
                     effect.type = "var",  # options: "mean", "var", "both", "None"
                     seed = 1234, # set seed
                     hi_corr_prop = 0.5, # fraction of regions that are "high correlation"
                     hi_corr_range = c(0.6, 0.95), # draw rho ~ Uniform(a,b) for high-corr regions
                     lo_corr_range = c(0.05, 0.2), # draw rho ~ Uniform(a,b) for low-corr regions
                     corr.type = "cs" # correlation structure within region
) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  Ncpgs <- nregions * cpgs.per.region
  
  group1_size <- round(nsamples * (grp.split.pc / 100))
  group2_size <- nsamples - group1_size
  
  # Prior parameters
  d0 <- 20
  s0 <- 0.8
  
  # Number of differential regions (DMR/VMR/DVMR)
  ndm_regions <- if (effect.type == "None" || pDM == 0) 0 else round(nregions * pDM)
  
  # Region IDs
  region_ids_unique <- paste0("Region", 1:nregions)
  
  # Assign regions to low/high correlation
  n_hi <- round(nregions * hi_corr_prop)
  n_lo <- nregions - n_hi
  
  # Randomly choose which regions are high-corr (rest are low-corr)
  hi_idx <- sort(sample.int(nregions, n_hi, replace = FALSE))
  lo_idx <- setdiff(seq_len(nregions), hi_idx)
  
  # Sample a rho per region within its range
  stopifnot(hi_corr_range[1] < hi_corr_range[2], lo_corr_range[1] < lo_corr_range[2])
  rho_vec <- numeric(nregions)
  if (length(hi_idx) > 0) rho_vec[hi_idx] <- runif(length(hi_idx), min = hi_corr_range[1], max = hi_corr_range[2])
  if (length(lo_idx) > 0) rho_vec[lo_idx] <- runif(length(lo_idx), min = lo_corr_range[1], max = lo_corr_range[2])
  
  if (ndm_regions > 0) {
    ndm_regions <- min(ndm_regions, length(hi_idx))
    dm_regions_idx <- sample(hi_idx, ndm_regions, replace = FALSE)
    dm_regions <- region_ids_unique[dm_regions_idx]
  } else {
    dm_regions <- character(0)
  }
  
  # Assign regions to CpGs
  region_ids <- rep(region_ids_unique, each = cpgs.per.region)
  
  # Initialize region variances (per group)
  sigma2_region1 <- d0 * s0^2 / stats::rchisq(nregions, d0)
  sigma2_region2 <- sigma2_region1
  
  # Apply variance effects to DM regions if needed
  if (effect.type %in% c("var", "both") && length(dm_regions) > 0) {
    dm_idx <- which(region_ids_unique %in% dm_regions)
    sigma2_region2[dm_idx] <- var.effect * d0 / stats::rchisq(length(dm_idx), d0)
  }
  
  # Expand per-CpG variances
  sigma2 <- cbind(rep(sigma2_region1, each = cpgs.per.region),
                  rep(sigma2_region2, each = cpgs.per.region))
  
  # Generate M values
  y <- matrix(0, ncol = nsamples, nrow = Ncpgs)
  
  current_cpg <- 1
  for (region_idx in seq_len(nregions)) {
    region_name <- region_ids_unique[region_idx]
    
    # Mean shift for group 2
    if (effect.type == "None") {
      mean_shift_group2 <- 0
    } else if (effect.type == "var") {
      # tiny offset to avoid exactly identical means under pure variance effects
      mean_shift_group2 <- 0.00000011985
    } else if (effect.type %in% c("mean", "both") && region_name %in% dm_regions) {
      mean_shift_group2 <- mean.effect
    } else {
      mean_shift_group2 <- 0
    }
    
    # Build region-specific correlation matrix
    rho <- rho_vec[region_idx]
    
    # Compound symmetric (exchangeable)
    corr_mat <- matrix(rho, nrow = cpgs.per.region, ncol = cpgs.per.region)
    diag(corr_mat) <- 1
    
    base_mean <- if (runif(1) < 0.5) -2 else 2
    for (sample_idx in 1:nsamples) {
      add_mean <- if (sample_idx > group1_size) mean_shift_group2 else 0
      mvn_vals <- MASS::mvrnorm(
        1,
        mu = rep(base_mean + add_mean, cpgs.per.region),
        Sigma = corr_mat
      )
      y[current_cpg:(current_cpg + cpgs.per.region - 1), sample_idx] <- mvn_vals
    }
    current_cpg <- current_cpg + cpgs.per.region
  }
  
  # Apply group-specific variances
  y[, 1:group1_size] <- sqrt(sigma2[, 1]) * y[, 1:group1_size]
  y[, (group1_size + 1):nsamples] <- sqrt(sigma2[, 2]) * y[, (group1_size + 1):nsamples]
  
  # Metadata
  ids <- paste0("Sample", 1:nsamples)
  group <- as.factor(c(rep(0, group1_size), rep(1, group2_size)))
  metadata <- data.frame(Sample = ids, Exposure = group)
  mvals <- as.data.frame(y)
  rownames(mvals) <- paste0("cpg", 1:Ncpgs)
  colnames(mvals) <- ids
  
  # CpG Truth
  CpG.Truth <- ifelse(region_ids %in% dm_regions,
                      effect.type,
                      "None")
  
  # Compute effects
  sim.mu.eff <- rowMeans(mvals[, (group1_size + 1):nsamples, drop = FALSE]) -
    rowMeans(mvals[, 1:group1_size, drop = FALSE])
  sim.var.eff <- matrixStats::rowVars(as.matrix(mvals[, (group1_size + 1):nsamples, drop = FALSE])) /
    matrixStats::rowVars(as.matrix(mvals[, 1:group1_size, drop = FALSE]))
  
  CpG.Info <- data.frame(CpGs = rownames(mvals),
                         Region = region_ids,
                         Mean.Effect = sim.mu.eff,
                         Var.Effect  = sim.var.eff,
                         Truth = CpG.Truth,
                         stringsAsFactors = FALSE)
  
  # Region-wise empirical correlations (average pairwise CpG correlation)
  region_cors <- CpG.Info %>%
    dplyr::select(CpGs, Region) %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(
      Correlation = {
        cpgs_in_region <- CpGs
        if (length(cpgs_in_region) < 2) {
          NA_real_
        } else {
          cor_mat <- stats::cor(t(mvals[cpgs_in_region, , drop = FALSE]))
          mean(cor_mat[upper.tri(cor_mat)])
        }
      },
      .groups = "drop"
    )
  
  # Attach the intended rho used per region
  rho_df <- data.frame(Region = region_ids_unique, IntendedRho = rho_vec, stringsAsFactors = FALSE)
  
  CpG.Info <- CpG.Info %>%
    dplyr::left_join(region_cors, by = "Region") %>%
    dplyr::left_join(rho_df, by = "Region")
  
  # True CpGs
  dm.cpgs <- CpG.Info$CpGs[CpG.Info$Truth != "None"]
  
  list(mval = mvals,
       metadata = metadata,
       CpG.Info = CpG.Info,
       dm.cpgs = dm.cpgs,
       hi_regions = region_ids_unique[hi_idx],
       lo_regions = region_ids_unique[lo_idx],
       dm_regions = dm_regions)
}
