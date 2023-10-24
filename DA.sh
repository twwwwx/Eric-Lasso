log_name=$(date +"%m%d%H")
Rscript real_DA/RealDataAnalysis.R > log/${log_name}RDA.log 2>&1 &