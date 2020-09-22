# [Re] Modeling Insect Phenology Using Ordinal Regression and Continuation Ratio Models
## Replication of Dennis et al. (1986) and Candy (1991)
<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4012772.svg)](https://doi.org/10.5281/zenodo.4012772)-->

This project reimplements a number of ordinal regression models from Dennis et al. (1986) and Candy (1991) in `R` and attempts to replicate their application to insect phenology data. The full references to the original articles are:

> B. Dennis, W.P. Kemp,and R.C. Beckwith (1986): "Stochastic Model of Insect Phenology: Estimation and Testing." Environmental Entomology 15(3):540–546. https://doi.org/10.1093/ee/15.3.540
> S.G. Candy (1991): "Modeling Insect Phenology Using Ordinal Regression and Continuation Ratio Models." Environmental Entomology 20(1):190–195. https://doi.org/10.1093/ee/20.1.190

The reproduction was successful and has been submitted to [ReScience C](https://rescience.github.io/) 

### Model reproduction

All `R` scripts required to reformat the original data and reproduce the statistical models are in the `code/` subdirectory.
They assume that the project folder is used as the working directory.
Additional `SAS` and `GLIM` code from the original publications is provided in the `code/` directory for completeness.


### Data
The `data/` subdirectory contains the dataset that was used in both original publications in its original form, as well as three reformatted versions that fulfill the input requirements for the replication code.   

#### Initial setup
The project makes use of the the `renv` package for project-local R dependency management. To reproduce the replication, I recommend to follow the following steps, assuming *R v3.6.0 or later* is already installed:

1. Clone this repository
```
git clone https://github.com/pboesu/replication_candy_1991
```
2. Launch *R* in the repository.

3. Run the following commands to install the exact versions of packages (as specified in the Project Environment)
```julia
import Pkg; Pkg.activate(".")
Pkg.instantiate()
```

#### Running the model
After setup, the analysis can be run to reproduce all figures and tables by sourcing the `R` scripts in sequence

```r
profvis::profvis({
source('code/01_reformat_budworm_data.R')
source('code/02_fit_candy_1991_models.R')
source('code/03_fit_dennis_et_al_model.R')
source('code/04_fit_VGAM_models.R')
source('code/05_fit_ordinal_models.R')
source('code/06_collate_parameter_tables.R')
})
```

The resulting figures are written to the `figures/` subdirectory, the resulting tables are written to the `outputs/` subdirectory.


Data formatting and model fitting are not computationally intensive. On an ordinary laptop, running the analysis scripts takes less than 1 minute, and uses less than 200MB of memory on a single core. Installing the packages before first time use may take an additional few minutes depending on the end-user's internet connection and whether or not packages are installed from source.


#### Notes on warnings and errors

```
In dpois(data$count, pred_n, log = T) : NaNs produced
```

- This warning is expected when running the `code/02_fit_candy_1991_models.R` script. It will be repeated several times, as the model likelihood is not finite under some parameter proposals before convergence is achieved.

```
glm.fit: fitted probabilities numerically 0 or 1 occurred
```
- This warning is expected when running the `code/02_fit_candy_1991_models.R` script. It occurs because complete or partial separation occurs for some stages in the continuation ratio model.

```
Error in checkwz(wz, M, trace = trace, wzepsilon = control$wzepsilon) : 
  NAs found in the working weights variable 'wz'
```

- This error is expected when running `code/04_fit_VGAM_models.R` as models fitted to the data with the `VGAM` package did not converge when using the cloglog link. In these cases the replication relied on alternative fitting approaches in `code/02_fit_candy_1991_models.R` and `code/05_fit_ordinal_models.R`.

### Article reproduction

The article uses the [ReScience C](https://rescience.github.io/) journal template. All elements are in the `article/` subfolder. Instructions to reproduce the article are provided in the subfolder [README](article/README.md).