#!/bin/bash -l

source $PWD/validation_hfr.ini

IFS=',' read -r -a name_exp_array <<< "$name_exp"

if [ $want_destag == "yes" ] ; then
    bsub -K -n 1 -q s_long -J HFRVAL -e aderr_0 -o adout_0 -P 0510 "python destaggering_UV.py $ini_date $fin_date $path_to_model_files $path_to_destag_model_files $time_res_model $path_to_mesh_mask ${name_exp_array[@]}" &
fi

wait

ncrcat $path_to_destag_model_files/${name_exp}_${time_res_model}*_grid_U2T.nc $path_to_destag_model_files/${name_exp}_${time_res_model}_grid_U2T_combined.nc
ncrcat $path_to_destag_model_files/${name_exp}_${time_res_model}*_grid_V2T.nc $path_to_destag_model_files/${name_exp}_${time_res_model}_grid_V2T_combined.nc

u_combined=$path_to_destag_model_files/${name_exp}_${time_res_model}*_grid_U2T_combined.nc
v_combined=$path_to_destag_model_files/${name_exp}_${time_res_model}*_grid_V2T_combined.nc

bsub -K -n 1 -q s_long -J HFRVAL -e aderr_1 -o adout_1 -P 0510 "python model_validation_hfr.py $ini_date $fin_date $path_to_hfr_files $path_to_mesh_mask $time_res_to_average $path_to_model_files $path_to_destag_model_files $want_destag $time_res_model $name_exp $output_plot_folder $u_combined $v_combined $work_dir" &
