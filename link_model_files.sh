#!/bin/bash -l

#source $PWD/validation_hfr.ini

component=$1
old_name_exp=$2
new_name_exp=$3
time_res_model=$4
path_to_model_files=$5
ini_date=$6
fin_date=$7
out_model=$8

current=$(date -d "$ini_date")
end=$(date -d "$fin_date +1 day")
echo $end
mkdir -p ${out_model}
#loop over all dates
while [ "$(date -d "$end" +%Y%m%d)" != "$(date -d "$current" +%Y%m%d)" ] 
do
    #fileU=${path_to_model_files}/${name_exp}_${time_res_model}_gridV25h_$(date -d "$current" +%Y%m%d)-$(date -d "$current" +%Y%m%d).nc
    file=`find ${path_to_model_files} -name "${old_name_exp}_${time_res_model}_grid${component}25h_$(date -d "$current" +%Y%m%d)-$(date -d "$current" +%Y%m%d).nc"`
    ln -sf $file ${out_model}/${new_name_exp}_${time_res_model}_$(date -d "$current" +%Y%m%d)_grid_${component}.nc
    current=$(date -d "$current +1 day")
    echo $(date -d "$current" +%Y%m%d)
    echo $current
done
