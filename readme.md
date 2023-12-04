# Outline

```
├── README.md
├── data # for plots
├── simulation
│   ├── eval.R #evaluation
│   ├── gereate_data.R # data function
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
│   ├── RDAresults.csv # real data
├── plot.Rmd # plots
├── action.sh # bash script for running simulation

``````
# Simulation

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
- `n`: sample size
- `p`: dimension
- `rho`: correlation
- `tau`: sdandard deviation of measurement error
- `sigma`: sdandard deviation of true covariance matrix


3. simulation parameters:
- `N_sim`: number of simulations
- `N_bs`: number of bootstrap size
  
### Outputs
Files will be automatedly saved in `results/results_table.csv`.
If you are running it on linux (which is not necessary), run the bash script
```
bash action.sh
```
The output log will be automatedly saved in `log/month_day_hour.log`.
