# get sys time month/day/hour
# month=`date +%m`
# day=`date +%d`
# hour=`date +%H`
# echo $month $day $hour
log_name=$(date +"%m%d%H")
Rscript simulation/eval.R > log/${log_name}.log 2>&1 &