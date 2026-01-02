# ReMeDy

`ReMeDy` A flexible statistical framework for region-based detection of DNA methylation dysregulation.

## ReMeDy Model Pipeline

<figure>
<img src="ReMeDy%20Pipeline.png" alt="Pipeline" />
</figure>

## Installation

For the installation, the R package devtools is needed.

``` r
install.packages('devtools')
library(devtools)
```

For the implementation of ReMeDy, the following R packages are needed to be installed beforehand.

For CRAN packages:

``` r
install.packages(c("data.table", "foreach", "doParallel", "ClustGeo", "readr", "dplyr", "dendextend", "peakRAM", "hglm", "tidyr", "pbapply"))
```

For Bioconductor packages:

``` r
install.packages("BiocManager")
BiocManager::install("bumphunter")
```

After installing these packages, ReMeDy can be installed by using devtools.

``` r
devtools::install_github("SChatLab/ReMeDy")
```

## Arguments

`ReMeDy` has 8 main arguments which the user needs to supply/define.

### The table below details the required arguments:

| Parameter   | Default  | Description                                                                                                                                                                                                                                                         |
|:------------|:--------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| region_list    |          | A list of co-methylated regions of methylation levels with CpGs in the row and samples in the column.|
| metadata    |          | A dataframe with information on the samples such as disease status, gender and etc.|
| expVar | Exposure | The name of the variable on which the differential methylation, variable methylation or both will be evaluated. If the user provides no name then the metadata should have a column named as exposure which should have information on things such as disease status, treatment conditions or etc.|
| coVars | NULL | The names of the covariates/confounders that needs to adjusted for in the DNA methylation analysis.|
| rand_effect |     | This parameter specifies the random-effects structure of the model.|
| parallel | FALSE | 	If true, the analysis will be performed on multiple cores with faster runtimes.|
| verbose | FALSE | 	If true, it will print the progress report.|
| numCores | 1 | The number of cores on which the analysis will be serially performed.|

## Values

`ReMeDy` outputs has 12 value arguments.

### The table below details the values returned by Sacoma:

| Value   | Description                                                                                       |
|:--------|:--------------------------------------------------------------------------------------------------|
| cpgs | Identifier of the CpG probe included in the analysis.|
| Region | Identifier of the region to which the CpG belongs.|
| Pval_mu | Region-level p-value for differential methylation (mean effect) derived from the mean model.|
| adjPval_mu | Benjamini–Hochberg adjusted p-value for differential methylation (DMR).|
| est_mu | Estimated exposure effect of mean methylation for the region (on the M-value scale).|
| SE_mu | Standard error of the estimated mean effect.|
| Pval_disp | Region-level p-value for variable methylation derived from the dispersion model.|
| adjPval_disp | Benjamini–Hochberg adjusted p-value for variable methylation (VMR).|
| est_disp | Estimated exposure effect of methylation variability for the region.|
| SE_disp | Standard error of the estimated variance effect.|
| Pval_mudisp | Joint p-value for combined mean and variance effects obtained using the Cauchy Combination Test.|
| adjPval_mudisp | Benjamini–Hochberg adjusted joint p-value for identifying DVMRs.|

## Working Example

### Importing libraries

``` r
packages <- c('foreach', 'doParallel', 'data.table', 'tidyverse', 'hglm', 'bumphunter', 'ENmix',
              'readr', 'ClustGeo', 'dendextend', 'dplyr', 'peakRAM', 'tidyr', 'gtools', 'ReMeDy')
for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}
```

### Loading Example Data

Loading a sample data

``` r
load("~/Data/sample_data.RData")
beta_vals <- data$beta_vals
metadata <- data$metadata
```

### Snapshot of data

A typical DNA methylation beta value data frame looks something
like below. Note, the CpGs should be in the rows and the samples in
columns.

``` r
beta_vals[1:5, 1:5]
```

    ##              GSM1088602 GSM1088603 GSM1088604 GSM1088605 GSM1088606
    ##  cg00000029  0.1042731 0.09970492 0.08676946 0.09811605  0.1202126
    ##  cg00000108  0.9521152 0.94776479 0.94435024 0.95042081  0.9413818
    ##  cg00000109  0.9071243 0.89921650 0.90713672 0.91270799  0.9044580
    ##  cg00000165  0.5675178 0.69421567 0.76857651 0.76790992  0.6425491
    ##  cg00000236  0.8681640 0.85634641 0.88599162 0.86816904  0.8576151

### Identifying co-methylated regions using SACOMA

To get co-methylated regions as input for ReMeDy you need to run SACOMA which can be found at **[SACOMA](https://github.com/SChatLab/SACOMA/tree/main/SACOMA%20Implementation)**.
Additionally, you will also need array manifest files for the Infinium 450k or EPIC platforms.  
These can be downloaded from following links:

- **[450k manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EUDB4zlBM7NPt6BLLVwiq8sBuT1pdhKPcsWlXuD_-jNc2A?e=UlUx8U)**
- **[EPIC manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EWOBny4uHHpDhef2_0OQcisBxvtm0kHPMBX9zzNK4Y-puQ?e=d4dJGZ)**  

Download the files and point the `manifest.dir` argument to the folder where they are stored when running the function.

To obtain co-methylated regions using SACOMA please use the code below.

``` r
sacoma.regions <- Sacoma(bval.dnam = beta_vals,
                         data_type = "450k",
                         manifest.dir = paste(dir, "Data", sep = '/'),
                         minCpGs = 3,
                         minCpGs.thres.type = "gr.equal",
                         maxGap = 200,
                         rthresold = 0.5, 
                         method = "spearman",
                         ncores = 1)
```

The results table from SACOMA should look something like the following:

``` r
sacoma.regions[1:5,]
```

    ##        Probe_ID  CHR     POS bin prediction Region Time.min peak.memory.gb
    ##    1 cg01537494 chr1 1102276  53      noise      1    0.847      0.2583691
    ##    2 cg02825344 chr1 1102294  53      noise      1    0.847      0.2583691
    ##    3 cg02331673 chr1 1102389  53      noise      1    0.847      0.2583691
    ##    4 cg03003335 chr1 1229351  86      noise      2    0.847      0.2583691
    ##    5 cg01937247 chr1 1229534  86      noise      2    0.847      0.2583691

Filter the results table by prediction = 'signal' to obtain co-methylated regions and then create a list of these regions, which will be the input for ReMeDy.

``` r
beta_vals.sub <- beta_vals[sacoma.regions$Probe_ID[sacoma.regions$prediction == 'signal'], ]
mvals.sub <- ENmix::B2M(beta_vals.sub)
mvals.sub$CpGs <- rownames(mvals.sub)
merged <- merge(sacoma.regions, mvals.sub, by.x = 'Probe_ID', by.y = 'CpGs')
merged <- merged[mixedorder(merged$Region),]
region_list <- split(merged, merged$Region)

region_list_clean <- lapply(region_list, function(x) {
  rownames(x) <- x$Probe_ID
  x <- x[, !(names(x) %in% c("CpGs", "Probe_ID", "CHR", "POS", "bin", "Region", "prediction", "Time.min", "peak.memory.gb"))]
  return(x)
})

names(region_list_clean) <- paste0('Region', 1:length(region_list_clean))
```

### Identifying DMRs, VMRs and DVMRs using ReMeDy

``` r
fit <- ReMeDy(region_list = region_list_clean,
              metadata = metadata,
              expVar = "Exposure",
              coVars = NULL,
              rand_effect = "(1 | Sample)",
              parallel = FALSE,
              verbose = FALSE,
              numCores = 1)
```

After performing analysis you can extract the results table in the following manner

``` r
results <- fit$res
```

The results table from ReMeDy should look something like the following

``` r
results[1:5,]
```

    ##           cpgs  Region      Pval_mu   adjPval_mu     est_mu     SE_mu Pval_disp adjPval_disp    est_disp   SE_disp  Pval_mudisp adjPval_mudisp
    ##    1 cg01261194 Region1 8.724734e-08 1.145284e-05 -0.2651702 0.0455216 0.1724971    0.8936299  0.58841011 0.4313149 1.744946e-07   2.290118e-05
    ##    2 cg01913436 Region1 8.724734e-08 1.145284e-05 -0.2651702 0.0455216 0.1724971    0.8936299  0.58841011 0.4313149 1.744946e-07   2.290118e-05
    ##    3 cg03748986 Region1 8.724734e-08 1.145284e-05 -0.2651702 0.0455216 0.1724971    0.8936299  0.58841011 0.4313149 1.744946e-07   2.290118e-05
    ##    4 cg01215682 Region2 1.270749e-05 1.836810e-04 -0.3132433 0.0677032 0.9593214    0.9989061 -0.02213983 0.4340699 2.542288e-05   3.368531e-04
    ##    5 cg03448635 Region2 1.270749e-05 1.836810e-04 -0.3132433 0.0677032 0.9593214    0.9989061 -0.02213983 0.4340699 2.542288e-05   3.368531e-04


## Analysis Scripts and Data Availability

All R scripts used for the simulation study and real-data analysis in the manuscript are available in this repository at the following locations:

- **Simulation analysis scripts:** [Simulations](https://github.com/SChatLab/ReMeDy/tree/Manuscript_Materials/Simulations)
- **Real data analysis scripts:** [Real Data Analysis](https://github.com/SChatLab/ReMeDy/tree/Manuscript_Materials/Real%20Data)

The datasets used in the paper were obtained from GEO under the following accession IDs:

- **450k datasets:** GSE44667 (EOPET dataset), GSE80970 (Alzheimer dataset)
- **EPIC dataset:** GSE306095 (B-ALL dataset)

After downloading, place the data files in a preferred directory and update the paths inside the analysis scripts accordingly.
