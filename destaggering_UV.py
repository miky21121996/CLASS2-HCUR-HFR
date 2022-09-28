import sys
import os.path
from os import listdir
from os.path import isfile, join
from datetime import timedelta, date
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xarray as xr
import matplotlib as mpl
mpl.use('Agg')

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)+1):
        yield start_date + timedelta(n)

def destaggering(date_in, date_fin, path_to_mod_output, path_to_destag_output_folder, name_exp, time_res, path_to_mask):

    isExist = os.path.exists(path_to_destag_output_folder)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path_to_destag_output_folder)
        print("The new directory is created!")

    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(path_to_mod_output,followlinks=True):
        print("dirpath: ", dirpath)
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    print(listOfFiles)
    start_date = date(int(date_in[0:4]),int(date_in[4:6]) , int(date_in[6:8]))
    end_date = date(int(date_fin[0:4]),int(date_fin[4:6]) , int(date_fin[6:8]))

    mesh_mask = xr.open_dataset(path_to_mask)
    u_mask = mesh_mask.umask.values
    v_mask = mesh_mask.vmask.values
    t_mask = mesh_mask.tmask.values
    u_mask = np.squeeze(u_mask[0, 0, :, :])
    v_mask = np.squeeze(v_mask[0, 0, :, :])
    t_mask = np.squeeze(t_mask[0, 0, :, :])

    for single_date in daterange(start_date, end_date):
        print(single_date)
        timetag=single_date.strftime("%Y%m%d")
        counter = 0
        for name in name_exp:
            u_filename = name + "_" + time_res + "_" + timetag + "_grid_U.nc"
            v_filename = name + "_" + time_res + "_" + timetag + "_grid_V.nc"
            t_filename = name + "_" + time_res + "_" + timetag + "_grid_T.nc"
            # print(any(v_filename in r for r in listOfFiles))
            if any(u_filename in s for s in listOfFiles) and any(v_filename in r for r in listOfFiles) and any(t_filename in q for q in listOfFiles):
                matching_u = [u_match for u_match in listOfFiles if u_filename in u_match]
                matching_v = [v_match for v_match in listOfFiles if v_filename in v_match]
                matching_t = [t_match for t_match in listOfFiles if t_filename in t_match]
                U_current = xr.open_dataset(listOfFiles[listOfFiles.index(matching_u[0])])
                V_current = xr.open_dataset(listOfFiles[listOfFiles.index(matching_v[0])])
                T_current = xr.open_dataset(listOfFiles[listOfFiles.index(matching_t[0])])
                experiment=name
                break
            else:
                counter = counter + 1
                continue
        if counter==len(name_exp):
            continue
        
        lon = U_current.nav_lon.values
        lat = U_current.nav_lat.values
        [dim_t, dim_lat, dim_lon] = U_current.ssu.shape

        u = U_current.ssu.values
        v = V_current.ssv.values
        #u = np.squeeze(u[:, :, :])
        #v = np.squeeze(v[:, :, :])

        #destaggering of u
        sum_u_mask = u_mask[:, 1:]+u_mask[:, :(dim_lon-1)]
        sum_u_mask = np.repeat(sum_u_mask[np.newaxis, :, :], dim_t, axis=0)
        sum_u = u[:, :, 1:]+u[:, :, :(dim_lon-1)]
        denominator_u_mask = np.maximum(sum_u_mask, 1)
        destaggered_u = np.zeros(u.shape)
        destaggered_u[:,:, 1:] = sum_u / denominator_u_mask
        destaggered_u=destaggered_u * t_mask
        #destaggered_u[destaggered_u == 0] = 'nan'

        #destaggering of v
        sum_v_mask = v_mask[1:, :]+v_mask[:(dim_lat-1), :]
        sum_v_mask = np.repeat(sum_v_mask[np.newaxis, :, :], dim_t, axis=0)
        sum_v = v[:, 1:, :]+v[:, :(dim_lat-1), :]
        denominator_v_mask = np.maximum(sum_v_mask, 1)
        destaggered_v = np.zeros(v.shape)
        destaggered_v[:, 1:, :] = sum_v / denominator_v_mask
        destaggered_v=destaggered_v*t_mask
        #destaggered_v[destaggered_v == 0] = 'nan'

        # save destaggered u in nc file
        U_for_attributes = U_current
        destaggered_U_current = T_current
        destaggered_U_current = destaggered_U_current.assign(destaggered_u=(('time_counter', 'y', 'x'),destaggered_u))
        destaggered_U_current.destaggered_u.attrs=U_for_attributes.ssu.attrs
        #destaggered_U_current = destaggered_U_current.drop(("ssu"))
        destaggered_U_current.to_netcdf(path_to_destag_output_folder + experiment + "_" + "1h" + "_" + timetag + "_grid_U2T.nc")

        # save destaggered v in nc file
        V_for_attributes = V_current
        destaggered_V_current = T_current
        destaggered_V_current = destaggered_V_current.assign(destaggered_v=(('time_counter', 'y', 'x'),destaggered_v))
        destaggered_V_current.destaggered_v.attrs=V_for_attributes.ssv.attrs
        #destaggered_V_current = destaggered_V_current.drop(("ssu"))
        destaggered_V_current.to_netcdf(path_to_destag_output_folder + experiment + "_" + "1h" + "_" + timetag + "_grid_V2T.nc")

if __name__ == "__main__":

    argv=sys.argv
    print("argv: ", argv)
    ini_date=argv[1]
    fin_date=argv[2]
    path_to_model_files=argv[3]
    path_to_destag_model_files=argv[4]
    time_res_model=argv[5]
    path_to_mesh_mask=argv[6]
    name_exp=argv[7:]
        
    destaggering(ini_date,fin_date,path_to_model_files,path_to_destag_model_files,name_exp,time_res_model,path_to_mesh_mask)
