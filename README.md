# ReMeDy

This repository contains all the R codes for implementing ReMeDy and reproducing the findings that can be found in our research article - ReMeDy: A flexible statistical framework for region-based detection of DNA methylation dysregulation.

## ReMeDy Model Pipeline

<figure>
<img src="ReMeDy%20Pipeline.png" alt="Pipeline" />
</figure>

## Implementation

For the implementation of ReMeDy, the following R packages are needed to be installed beforehand.

For CRAN packages:

``` r
install.packages(c("data.table", "foreach", "doParallel", "ClustGeo", "readr", "dplyr", "dendextend", "peakRAM", "hglm"))
```

For Bioconductor packages:

``` r
install.packages("BiocManager")
BiocManager::install("bumphunter")
```

The main function used in this workflow can be found: **[ReMeDy Function](https://github.com/SChatLab/ReMeDy/tree/main/ReMeDy%20Implementation)** 

To get co-methylated regions as input for ReMeDy you need to run SACOMA which can be found at **[SACOMA](https://github.com/SChatLab/SACOMA/tree/main/SACOMA%20Implementation)**.
Additionally, you will also need array manifest files for the Infinium 450k or EPIC platforms.  
These can be downloaded from following links:

- **[450k manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EUDB4zlBM7NPt6BLLVwiq8sBuT1pdhKPcsWlXuD_-jNc2A?e=UlUx8U)**
- **[EPIC manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EWOBny4uHHpDhef2_0OQcisBxvtm0kHPMBX9zzNK4Y-puQ?e=d4dJGZ)**  

Download the files and point the `manifest.dir` argument to the folder where they are stored when running the function.

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
