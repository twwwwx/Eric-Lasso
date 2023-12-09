# Outline

```
├── README.md
├── data # for plots and real data analysis
├── simulation
│   ├── eval.R # simulation
│   ├── gereate_data.R # data function
│   ├── analysis4combo.R # real data analysis
├── R
│   ├── eic.R # main function
│   ├── ADMM_proj.R # ADMM projection
│   ├── cv_covariance_matrices.R # calibrated CV
│   ├── cv_covariance_noproj.R # CV
│   ├── lasso_covariance.R # lasso solution
│   ├── lasso_covariance_constrained.R 
            #lasso solution with constraint
├── results
│   ├── results_table.csv # main simulation
│   ├── results_sum.csv # sums in simulation
│   ├── RDAresults.csv # real data results
├── plot.Rmd # plots
├── action.sh # bash script for running simulation

``````
# Eric-Lasso
This is the code for the paper "High-dimensional regression analysis of compositional covariates with measurement errors". It implements four methods for high-dimensional regression analysis of compositional covariates with measurement errors, including Eric, [CoCo](https://arxiv.org/pdf/1510.07123.pdf), [Coda](https://academic.oup.com/biomet/article/101/4/785/1775476) and vanilla lasso. The tuning parameter selection in Eric and CoCo Lasso was done by using the [calibrated cross validation](https://arxiv.org/pdf/1510.07123.pdf) method and that in Coda and Vani Lasso was done by using cross validation. 
This code is built on the [BDcocolasso](https://github.com/celiaescribe/BDcocolasso) package. 


# Run simulation

For the simulation, we only need to run `eval.R`.

Required R package:
```
install.packages("MASS","boot","rlist","emdbook","dirmult")
```
### Parameters in `eval.R`

1. Model parameters:

- `constrain=T, proj=T`: Eric

- `constrain=F, proj=T`: CoCo

- `constrain=T, proj=F`: Coda

- `constrain=F, proj=F`: vanilla lasso

2. data parameters:
- `data_type`: chose from `lognormal` `dirichlet` `dirmult`
- `n`: sample size
- `p`: dimension
- `rho`: correlation
- `tau`: sdandard deviation of measurement error
- `sigma`: sdandard deviation of true covariance matrix


1. simulation parameters:
- `N_sim`: number of simulations
  
### Outputs
Files will be automatedly saved in `results/results_table.csv`.
If you are running it on linux (which is not necessary), create a new directory named `log` then run the bash script
```
bash action.sh
```
The output log will be automatedly saved in `log/month_day_hour.log`.
