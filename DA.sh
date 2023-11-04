log_name=$(date +"%m%d%H")
Rscript real_DA/RealDataAnalysis_full.R > RDAlog/${log_name}RDA.log 2>&1 &