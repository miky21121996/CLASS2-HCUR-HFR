import sys
import numpy as np
import xarray
from netCDF4 import Dataset

def CreateCoordinatesFile(mesh_mask,name_exp,path_to_output_folder):
    coord_file=xarray.Dataset(data_vars=dict(glamt=(['y','x'], mesh_mask['glamt'].data[0,:,:]),glamu=(['y','x'], mesh_mask['glamu'].data[0,:,:]),glamv=(['y','x'], mesh_mask['glamv'].data[0,:,:]),glamf=(['y','x'], mesh_mask['glamf'].data[0,:,:]),gphit=(['y','x'], mesh_mask['gphit'].data[0,:,:]),gphiu=(['y','x'], mesh_mask['gphiu'].data[0,:,:]),gphiv=(['y','x'], mesh_mask['gphiv'].data[0,:,:]),gphif=(['y','x'], mesh_mask['gphif'].data[0,:,:]),tmask=(['y','x'], mesh_mask['tmask'].data[0,0,:,:]),umask=(['y','x'], mesh_mask['umask'].data[0,0,:,:]),vmask=(['y','x'], mesh_mask['vmask'].data[0,0,:,:])))
    new_filename = path_to_output_folder + '/' + name_exp + '_coordinates.nc'
    print ('saving ', new_filename)
    coord_file.to_netcdf(path=new_filename)

if __name__ == "__main__":

    argv=sys.argv
    print("argv: ", argv)
    path_to_original_mesh_mask=argv[1]
    coordinates_file=argv[2]
    path_to_output_folder=argv[3]
    name_exp=argv[4]

    #domain_cfg_ds = Dataset(path_to_domain_cfg)
    mesh_mask_ds = Dataset(path_to_original_mesh_mask)
    mesh_mask=xarray.Dataset(data_vars=dict(tmask=(['time_counter','z','y','x'], mesh_mask_ds['tmask']), nav_lon=(['y','x'], mesh_mask_ds['x']), nav_lat=(['y','x'], mesh_mask_ds['y']),nav_lev=(['z'], mesh_mask_ds['nav_lev']), umask=(['time_counter','nav_lev','y','x'], mesh_mask_ds['umask']),vmask=(['time_counter','nav_lev','y','x'], mesh_mask_ds['vmask']),glamt=(['time_counter','y','x'], mesh_mask_ds['glamt']),gphit=(['time_counter','y','x'], mesh_mask_ds['gphit']),glamu=(['time_counter','y','x'], mesh_mask_ds['glamu']),gphiu=(['time_counter','y','x'], mesh_mask_ds['gphiu']),glamv=(['time_counter','y','x'], mesh_mask_ds['glamv']),gphiv=(['time_counter','y','x'], mesh_mask_ds['gphiv']),glamf=(['time_counter','y','x'], mesh_mask_ds['glamf']),gphif=(['time_counter','y','x'], mesh_mask_ds['gphif'])))

    print ('saving to ', path_to_output_folder)
    mesh_mask.to_netcdf(path=path_to_output_folder+"/"+name_exp+"_mesh_mask_to_save.nc")

    if coordinates_file=='yes':
        CreateCoordinatesFile(mesh_mask,name_exp,path_to_output_folder)
