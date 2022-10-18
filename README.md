# HFR-HCURR

This tool provides statistics from comparison between model data and High Frequency radars current velocities. It is necessary to have a python virtual environment with the following packages: sys, warnings, netCDF4, datetime, collections, os, csv, numpy, xarray, pandas, scipy, matplotlib, cartopy.

1. Open validation_hfr.ini

INPUT VARIABLES:

  a. ini_date: initial date  
  b. fin_date: final date  
  c. path_to_hfr_files: path to .nc HFR files  
  d. path_to_mesh_mask: path to .nc mesh mask file  
  e. time_res_to_average: time resolution you want the tool use to average the current velocities data (es: 1M, monthly resolution)  
  f. path_to_model_files: path to .nc model files  
  g. want_destag: if yes the tool executes destaggering on model files, if no the tool skips destaggering (because you already have destaggered velocity .nc files for example and want just to run the plot part)  
  h. time_res_model: time resolution of nc model files (es: 1h)  
  i. name_exp: name of the experiment in model file names. If the model files of the experiment change name over time, write all the names separated by comma (es:   name_exp=mfs1,mfs2)  

OUTPUT VARIABLES:

  a. work_dir: path to working directory  
  b. path_to_destag_model_files: path to the folder where you want the destaggered files provided by the tool  
  c. output_plot_folder: path to the folder where you want to have the plots  

2. run sh main_validation.sh
