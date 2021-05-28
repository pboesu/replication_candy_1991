# [Re] Modeling Insect Phenology Using Ordinal Regression and Continuation Ratio Models
## Replication of Dennis et al. (1986) and Candy (1991)
<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4012772.svg)](https://doi.org/10.5281/zenodo.4012772)-->

This project reimplements a number of ordinal regression models from Dennis et al. (1986) and Candy (1991) in `R` and attempts to replicate their application to insect phenology data. The full references to the original articles are:

> B. Dennis, W.P. Kemp,and R.C. Beckwith (1986): "Stochastic Model of Insect Phenology: Estimation and Testing." Environmental Entomology 15(3):540–546. https://doi.org/10.1093/ee/15.3.540    
> S.G. Candy (1991): "Modeling Insect Phenology Using Ordinal Regression and Continuation Ratio Models." Environmental Entomology 20(1):190–195. https://doi.org/10.1093/ee/20.1.190

The reproduction was successful and has been submitted to [ReScience C](https://rescience.github.io/) 

### Data
The `data/` subdirectory contains the dataset that was used in both original publications in its original form, as well as three reformatted versions that fulfill the input requirements for the replication code.   

### Model fitting 
All `R` scripts required to reformat the original data and reproduce the statistical models are in the `code/` subdirectory.
They assume that the project folder is used as the working directory.    
Additional `SAS` and `GLIM` code from the original publications is provided in the `code/` directory for completeness.


#### Initial setup
The project makes use of the the `renv` package for project-local R dependency management. To reproduce the replication, I recommend to follow the following steps, assuming *R v3.6.0 or later* is already installed:

1. Clone this repository
```
git clone https://github.com/pboesu/replication_candy_1991
```
2. If you do not want to use `renv` remove or rename `.Rprofile` and skip to stage 6 (not recommended). 

3. If you want to use `renv` launch `R` in the repository. The local `.Rprofile` file should ensure that the correct version of `renv` is installed. At the time of writing it appears safe to ignore any message about using `renv::status()`, and any message stating that the projected requested an older R version than the currently installed one. 

4. Run the following command to install the exact versions of packages specified in the renv.lock file.
```r
renv::restore()
```
5. Agree to installing the required packages when prompted. Installing the packages to the project library before first time use may take up to 20 minutes depending on the end-user's internet connection and whether or not packages are installed from source.

6. If you do not choose to use `renv` for dependency management you can launch `R` in the repository root and manually install the necessary packages as listed in the session info below:
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] splines   stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
[1] ggplot2_3.3.3      ordinal_2019.12-10 VGAM_1.1-5         kableExtra_1.2.1   dplyr_1.0.2       

loaded via a namespace (and not attached):
 [1] latex2exp_0.4.0     Rcpp_1.0.5          pillar_1.4.6        compiler_3.6.1      tools_3.6.1         digest_0.6.25      
 [7] gtable_0.3.0        lattice_0.20-41     evaluate_0.14       lifecycle_0.2.0     tibble_3.0.3        viridisLite_0.3.0  
[13] ucminf_1.1-4        pkgconfig_2.0.3     rlang_0.4.7         Matrix_1.2-18       cli_2.0.2           rstudioapi_0.11    
[19] xfun_0.17           withr_2.4.2         httr_1.4.2          stringr_1.4.0       knitr_1.29          xml2_1.3.2         
[25] generics_0.0.2      vctrs_0.3.4         gtools_3.8.2        hms_0.5.3           grid_3.6.1          webshot_0.5.2      
[31] tidyselect_1.1.0    glue_1.4.2          R6_2.4.1            fansi_0.4.1         rmarkdown_2.3       tidyr_1.1.2        
[37] readr_1.3.1         purrr_0.3.4         magrittr_1.5        MASS_7.3-51.6       scales_1.1.1        ellipsis_0.3.1     
[43] htmltools_0.5.0     assertthat_0.2.1    rvest_0.3.6         colorspace_1.4-1    numDeriv_2016.8-1.1 renv_0.12.0        
[49] stringi_1.5.3       munsell_0.5.0       crayon_1.3.4   
```


#### Running the models
After setup, the complete analysis can be run to reproduce all figures and tables by sourcing the `R` scripts `01_reformat_budworm_data.R` to `code/07_initial_value_sensitivity.R` in sequence. 
However, each individual analysis script can be run individually and the user is encouraged to step through the individual scripts to learn more about the replication, and in particular to see where particular model implementations fail or produce warnings.
Note that the final script to reproduce the appendix will take about 15 minutes to run on a regular desktop.

```r
# data wrangling script
source('code/01_reformat_budworm_data.R')
# fit the models from Candy (1991) by direct optimisation of the likelihood
source('code/02_fit_candy_1991_models.R')
# fit the models from Dennis et al. (1986) by direct optimisation of the likelihood
source('code/03_fit_dennis_et_al_model.R')
# fit sequential and cumulative models using VGAM
source('code/04_fit_VGAM_models.R')
# fit cumulative and sequential models using ordinal
source('code/05_fit_ordinal_models.R')
# assemble parameter estimates for typesetting
source('code/06_collate_parameter_tables.R')
# conduct simulations to assess sensitivity to initial values
source('code/07_initial_value_sensitivity.R')
```

The reformatted data is written in CSV format to `data/`, the resulting figures are written to the `figures/` subdirectory, and the resulting tables are written as `.tex` files to the `outputs/` subdirectory. The manuscript can be recompiled using these.
In addition, several intermediate outputs are saved in `outputs/` as binary `.RDS` files. This is done to facilitate rerunning individual scripts while minimizing numerical inaccuracies caused by converting between numerical and textual representations of parameter values.

Data formatting and model fitting of the main manuscript are not computationally intensive. On an ordinary laptop, running the analysis scripts 01 - 06 takes less than 1 minute, and uses less than 200MB of memory on a single core. 
Simulations for the appendix, however, take about 15 minutes to execute on a single core.

#### Notes on warnings and errors

```
glm.fit: fitted probabilities numerically 0 or 1 occurred
```
- This warning is expected when running the `code/02_fit_candy_1991_models.R` script. It occurs because predicted probabilities are numerically indistinguishable from 1. This is further discussed in section 4.3 of the manuscript.

```
In checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon) :
  5 diagonal elements of the working weights variable 'wz' have been replaced by 1.819e-12
```

- This warning is expected when running `code/04_fit_VGAM_models.R`. It is related to the sparsity issues that result in the above warning produced by `stats::glm` and notifies the user that some diagonal elements of the weights matrix used in the iterative model fitting procedure have been adjusted to prevent numerical underflow. 
In addition, neither of the `vglm` models converged when using the cloglog link. This behaviour has been reported to the package maintainer as a potential bug. In these cases the replication relied on alternative fitting approaches in `code/02_fit_candy_1991_models.R` and `code/05_fit_ordinal_models.R`.

```
In dpois(data$count, pred_n, log = T) : NaNs produced
```

- This warning is produced during when running `code/07_initial_value_sensitivity.R`, because the code explicitly checks for false model convergence under infinite likelihoods.

### Article reproduction
The article uses the [ReScience C](https://rescience.github.io/) journal template. All elements are in the `article/` subfolder. Instructions to reproduce the article are provided in the subfolder [README](article/README.md).