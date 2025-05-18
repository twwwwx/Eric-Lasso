#!/bin/bash
#SBATCH --job-name=sim          # Job name
#SBATCH --time=5:00:00             # Time limit hrs:min:sec
#SBATCH --account=one
#SBATCH --output=/dev/null      # Prevent SLURM from handling output directly
#SBATCH --mem=1G                # Memory demand (e.g., 4GB)
#SBATCH --mail-type=END
#SBATCH --mail-user=wenccccy1003@gmail.com
module load r/4.4.2
log_name=$(date +"%m%d%H%M")
log_file="log/${log_name}.log"
# Redirect all output to the log file
exec > "$log_file" 2>&1

Rscript simulation/eval.R 
# Rscript simulation/analysis4combo.R > log/${log_name}_RDA.log 2>&1 &