#!/bin/bash -l

source $PWD/validation_hfr.ini

IFS=',' read -r -a time_res_model_array <<< "$time_res_model"
IFS=',' read -r -a want_destag_array <<< "$want_destag"
IFS=',' read -r -a name_exp_array <<< "$name_exp"
IFS=',' read -r -a model_array <<< "$path_to_model_files"
IFS=',' read -r -a out_model_array <<< "$path_to_out_model_files"
IFS=',' read -r -a mesh_mask_array <<< "$path_to_mesh_mask"
IFS=',' read -r -a label_plot_array <<< "$label_plot"
IFS=',' read -r -a irregular_grid_array <<< "$irregular_grid"
IFS=',' read -r -a output_plot_folder_array <<< "$output_plot_folder"
#IFS=',' read -r -a want_link_array <<< "$want_link"

u_combined_array=()
v_combined_array=()
for i in $(seq 0 $((num_exp-1)))
do
    echo "ciao"
    echo ${irregular_grid_array[$i]}
    echo ${model_array[$i]}
    if [ "${irregular_grid_array[$i]}" = "no" ] ; then
        if [ "${want_destag_array[$i]}" = "yes" ] ; then
            bsub -K -n 1 -q s_long -J HFRVAL -e aderr_0 -o adout_0 -P 0510 "python destaggering_UV.py $ini_date $fin_date ${model_array[$i]} ${out_model_array[$i]} ${time_res_model_array[$i]} ${mesh_mask_array[$i]} ${name_exp_array[$i]}" &
            wait
            ncrcat ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}*_grid_U2T.nc ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_U2T_combined.nc
            ncrcat ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}*_grid_V2T.nc ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_V2T_combined.nc
            
          
        fi
        u_combined_array[$i]=${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_U2T_combined.nc
        v_combined_array[$i]=${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_V2T_combined.nc
    fi
    if [ "${irregular_grid_array[$i]}" = "yes" ] ; then
#        if [ "${want_link_array[$i]}" = "yes" ] ; then
#            sh link_model_files.sh U ${name_exp_array[$i]} ${name_exp_array[$i]} ${time_res_model_array[$i]} ${model_array[$i]} ${ini_date} ${fin_date} ${out_model_array[$i]}
#            sh link_model_files.sh V ${name_exp_array[$i]} ${name_exp_array[$i]} ${time_res_model_array[$i]} ${model_array[$i]} ${ini_date} ${fin_date} ${out_model_array[$i]}
#            ncrcat ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}*_grid_U.nc ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_U_combined.nc
#            ncrcat ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}*_grid_V.nc ${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_V_combined.nc
#        fi
        u_combined_array[$i]=${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_U_combined.nc
        v_combined_array[$i]=${out_model_array[$i]}/${name_exp_array[$i]}_${time_res_model_array[$i]}_grid_V_combined.nc
    fi
done
TMP_DIR_array=()
TMP_DIR_array+=(${out_model_array[@]})
TMP_DIR_array+=(${irregular_grid_array[@]})
TMP_DIR_array+=(${mesh_mask_array[@]})
TMP_DIR_array+=(${name_exp_array[@]})
TMP_DIR_array+=(${label_plot_array[@]})
TMP_DIR_array+=(${output_plot_folder_array[@]})
TMP_DIR_array+=(${u_combined_array[@]})
TMP_DIR_array+=(${v_combined_array[@]})

bsub -K -n 1 -q s_long -J HFRVAL -e aderr_1 -o adout_1 -P 0510 "python model_validation_hfr.py $ini_date $fin_date $path_to_hfr_files $time_res_to_average $num_exp $work_dir $output_plot_folder_comparison ${TMP_DIR_array[@]}" &
