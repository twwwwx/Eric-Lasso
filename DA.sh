log_name=$(date +"%m%d%H")
# Rscript real_DA/RealDataAnalysis.R > log/${log_name}RDA.log 2>&1 &
Rscript real_DA/analysis4combo.R > RDAlog/${log_name}.log 2>&1 &