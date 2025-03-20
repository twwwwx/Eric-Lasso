# get sys time month/day/hour
# month=`date +%m`
# day=`date +%d`
# hour=`date +%H`
# echo $month $day $hour
module load r/4.4.2
log_name=$(date +"%m%d%H")
Rscript test/eval_test_2.R > test/${log_name}.log 2>&1 &