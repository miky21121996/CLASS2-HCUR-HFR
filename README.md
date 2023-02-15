# HFR-HCURR
This tool provides statistics comparing model data and high frequency radars current velocities.

## Table of Contents
* [General Info](#general-info)
* [Setup](#setup)
* [Usage](#usage)
* [Project Status](#project-status)
* [Contact](#contact)

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

If the model files of the experiment change name over time, you can use *link_model_files.sh* because the tool needs to have an unique name for each experiment:  
*sh link_model_files.sh component name_exp time_res_model path_to_model_files ini_date fin_date out_model*  
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

If your mesh mask lat and lon variables' names are "x" and "y", as the dimension names, it is necessary to rename them in "nav_lon" (depending on dimensions x and y) and "nav_lat" (depending on dimensions x and y). You can use *mesh_mask_and_coord_file.py*:  
*python mesh_mask_and_coord_file.py path_to_original_mesh_mask coordinates_file path_to_output_folder name_exp*  
Where:  
* path_to_original_mesh_mask: path to the mesh mask file from which you want to create the new renamed mesh mask file
* coordinates_file: if the model has an irregular grid you will need a coordinates.nc file to create an angles.nc file, necessary to execute the interpolation from the model irregular grid to the HF radar regular grid. If "no" the run will provide just the renamed mesh mask file. If "yes" the run will provide the renamed mesh mask file and the coordinates file
* path_to_output_folder: path to folder where you want new mesh mask and coordinates file to be saved
* name_exp: name you want your new mesh mask and coordinates file to have as prefix (es: name_exp_mesh_mask_to_save.nc and name_exp_coordinates.nc)

If the model has an irregular grid, you can use *provide_angles.F90*:  
*gfortran -o provide_angles.x provide_angles.F90 -I/zeus/opt/impi19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L/zeus/opt/impi19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1/lib -lnetcdff -lnetcdf*
*./provide_angles.x*  

In *provide_angles.F90* you must have:  
*  CHARACTER(LEN=16), PARAMETER :: conf = 'name_exp' (name_exp is the same of the prefix of coordinates.nc)  
*   INTEGER, PARAMETER :: jpiglo = number of x points and INTEGER, PARAMETER :: jpjglo = number of y points (you can check them in mesh mask or coordinates file dimensions)

The run will provide a name_exp_angle.nc file that you need to include in *zapata.py*.
In *def _resolve_grid(self, ingrid, grid,verbose=False):* you must add something like this:  
if ingrid == 'name_exp':  
    print(f' Tripolar 1/128 Grid -- {ingrid}')  
    grid = xr.open_dataset(self.mdir + 'name_exp_coordinates.nc')  
    angle = xr.open_dataset(self.mdir + '/name_exp_angle.nc')  
    struct={'tmask': grid.tmask, 'umask': grid.umask,'vmask': grid.vmask, 'tangle': angle.tangle, \  
            'lonT': grid.glamt,'latT':grid.gphit,'lonU':grid.glamu, \  
            'latU':grid.gphiu,'lonV':grid.glamv,'latV':grid.gphiv  }  
  
For what concerns the HF radar nc files, in *def _resolve_grid(self, ingrid, grid,verbose=False):* you need to include for each of them something like this:  
elif ingrid == 'GL_TV_HF_HFR-Ibiza-Total':  
    print(f' Regular HFR Lat-Lon Grid -- {ingrid}')  
    tk=np.logical_not(grid['mask'])  
    mask = tk.assign_coords({'lat':grid.nav_lat,'lon':grid.nav_lon})  
    struct={'mask': mask}  

## Usage

1. Open *validation_hfr.ini*

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

## Project Status
Project is: _in progress_ 

## Contact
Created by Michele Giurato (michele.giurato@cmcc.it) - feel free to contact me!
