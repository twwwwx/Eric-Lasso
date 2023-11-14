log_name=$(date +"%m%d%H")
mkdir -p RDAlog/${log_name}
# file_name = "real_DA/analysis4combo.R"
file_name="real_DA/analysis4combo_whole.R"
cp ${file_name} RDAlog/${log_name}
Rscript ${file_name} > RDAlog/${log_name}/${log_name}.log 2>&1 &