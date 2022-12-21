# HFR-HCURR
This tool provides statistics comparing model data and high frequency radars current velocities.

## Table of Contents
* [General Info](#general-info)
* [Setup](#setup)
* [Usage](#usage)
* [Project Status](#project-status)
* [Room for Improvement](#room-for-improvement)
* [Acknowledgements](#acknowledgements)
* [Contact](#contact)
<!-- * [License](#license) -->
## General Info
The tool can manage the comparison of multiple experiments (with regular and irregular grid) with the same observation dataset.

## Setup
It is necessary to have a python virtual environment with the following packages:  
* sys  
* warnings  
* netCDF4  
* datetime  
* collections  
* os  
* csv  
* numpy  
* xarray (version 0.20.1)
* pandas  
* scipy  
* matplotlib  
* cartopy  

To clone the repository: git clone git@github.com:miky21121996/CLASS2-HCUR-HFR.git name_your_repo

## Usage
1. If the model files of the experiment change name over time, you can use *link_model_files.sh* because the tool needs to have an unique name for each experiment: sh link_model_files.sh component name_exp time_res_model path_to_model_files ini_date fin_date out_model  
Where:  
* component: U or V
* old_name_exp: name of the exp you want to change to uniform all the model file names of the experiment  
* new_name_exp: name you have chosen to uniform the model file names
* time_res_model: time resolution of nc model files (es 1h)  
* path_to_model_fles: path to .nc model files  
* ini_date: initial date from which you want to link the files with the new name
* fin_date: final date up to which you want to link the files with the new name
* out_model: name of the folder that will be created to contain the renamed linked nc files  
Note:
* it is important to have in the out_model folder the model files for the entire period of interest.  
* It is likely you need to change the code in *link_model_files.sh* to adapt it to the original path you want to rename linking it. The final linked path must have this format: out_model/name_exp_${time_res_model}_$(date -d "$current" +%Y%m%d)_grid_${component}.nc

2. Open *validation_hfr.ini*

**INPUT VARIABLES**:

* ini_date: initial date  
* fin_date: final date  
* path_to_hfr_files: path to .nc HFR files (example of HFR file: /work/oda/mg28621/prova_destag/hfr_validation/hfr_data/hfr_data/GL_TV_HF_HFR-TirLig-Total.nc)
* path_to_mesh_mask: path to .nc mesh mask file correspondent to the model experiment; if more experiments separate the mesh mask paths with a comma without spaces  
* time_res_to_average: time resolution you want the tool to use to average the current velocities data (es: 1M, monthly resolution)  
* num_exp: number of experiments
* path_to_model_files: path to .nc model files (if more experiments separate the model paths with a comma without spaces)  
* irregular_grid: 'yes' if the experiment is associated to a model irregular grid, 'no' if the experiment is associated to a model regural grid; if more experiments separate the yes/no with a comma without spaces  
* want_destag: if 'yes' (model grid must be regular) the tool executes destaggering on model files, if 'no' the tool skips destaggering (because you already have destaggered velocity .nc files for example and want just to run the plot part or because the model grid is irregular); if more experiments, separate the yes/no with a comma without spaces  
* time_res_model: time resolution of nc model files (es: 1h); if more experiments separate the resolutions with a comma without spaces  
* name_exp: name of the experiment in model file names; if more experiments separate the names with a comma without spaces  
* label_plot: name of the experiment you want to be shown in plots; if more experiments separate the names with a comma without spaces  

**OUTPUT VARIABLES**:

* work_dir: path to working directory  
* path_to_out_model_files: path to the folder where you want the destaggered or linked files provided by the tool to be in; if more experiments separate the paths with a comma without spaces  
* output_plot_folder: path to the folder where you want to have the plots of the experiment; if more experiments separate the paths with a comma without spaces  
* output_plot_folder_comparison: path to the folder where you want to have the plots that directly compare the various experiments (num_exp must be larger than 1)

2. run sh main_validation.sh &
