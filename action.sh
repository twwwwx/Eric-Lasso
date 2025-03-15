# get sys time month/day/hour
# month=`date +%m`
# day=`date +%d`
# hour=`date +%H`
# echo $month $day $hour
module load r/4.4.2
log_name=$(date +"%m%d%H")
# Check for -n flag and its value
# while getopts "n:" opt; do
#   case $opt in
#     n)
#       log_name="$(date +"%m%d%H")_$OPTARG"
#       ;;
#     *)
      
#       ;;
#   esac
# done

Rscript simulation/eval.R > log/${log_name}.log 2>&1 &