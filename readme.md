### Outline

```
├── README.md
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
``````

For the simulation, we only need to run `eval.R`.

required R package:
```
install.packages("MASS","boot","rlist","emdbook")
```
### Parameters in `eval.R`

1. Model parameters:

- `constrain=T, proj=T`: eic

- `constrain=F, proj=T`: coco

- `constrain=T, proj=F`: coda

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
If you are running it on linux, run the bash script
```
bash action.sh
```
The output log will be automatedly saved in `results/month_day_hour.log`.
