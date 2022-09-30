import warnings
warnings.filterwarnings("ignore")
import sys
import datetime
from datetime import timedelta, date, datetime
from collections import OrderedDict
import os
import csv
import numpy as np
import xarray
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerBase
from matplotlib import colors
from matplotlib.colors import rgb2hex, Normalize, ListedColormap
import matplotlib.lines as mlines
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.dates as mdates
from matplotlib.colorbar import ColorbarBase
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import numpy.ma as ma
from SeaOverLand import seaoverland
from scipy.stats import linregress, pearsonr, gaussian_kde

def line_A(x, m_A, q_A):
    return (m_A*x+q_A)

def BIAS(data,obs):
    return  np.round((np.nanmean( data-obs)).data, 2)

def RMSE(data,obs):
    return np.round(np.sqrt(np.nanmean((data-obs)**2)),2)

def ScatterIndex(data,obs):
    num=np.sum(((data-np.nanmean(data))-(obs-np.nanmean(obs)))**2)
    denom=np.sum(obs**2)
    return np.round(np.sqrt((num/denom)),2)

def Normalized_std(data,obs):
    data_std=np.std(data)
    data_obs=np.std(obs)
    return np.round(data_std/data_obs,2)

def find_nearest(array, value, pprint=False):
    array = np.asarray(array)
    if pprint:
        print("model: ",array)
        print("obs: ",value)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def getSourceAntennas(ds):
    data = {'antennas':[] ,'bbox':[]}
    try:
        for name,lat,lon in zip(ds['SCDT'][0,:].astype('str').data,ds['SLTT'][0,:].data,ds['SLNT'][0,:].data):
            if np.isnan(lat) == False:
                name,lat,lon = name.strip(),round(lat, 4),round(lon, 4)
                data['antennas'].append({'name': name, 'lat': lat, 'lon': lon})
        data['bbox'] = [float(ds.geospatial_lon_min),float(ds.geospatial_lon_max),float(ds.geospatial_lat_min),float(ds.geospatial_lat_max)]
    except Exception as e:
        data = None
        print('An error ocurred when checking antennas')
    return data

def Get_List_Of_Files(path_to_hfr_files):
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(path_to_hfr_files,followlinks=True):
        print("dirpath: ", dirpath)
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    print(listOfFiles)
    return listOfFiles

def Get_String_Time_Resolution(start_date, end_date, time_res_to_average):
    dates = [start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")]
    start, end = [datetime.strptime(_, "%Y-%m-%d") for _ in dates]
    if time_res_to_average[-1] == 'D':
        string_time_res=list(OrderedDict(((start + timedelta(_)).strftime(r"%d-%b-%y"), None) for _ in range((end - start).days+1)).keys())
    if time_res_to_average[-1] == 'M':
        string_time_res=list(OrderedDict(((start + timedelta(_)).strftime(r"%b-%y"), None) for _ in range((end - start).days+1)).keys())

    return string_time_res

def Get_Closest_Hfr_Time_Range_Index(time_res_to_average,ini_date,fin_date,averaged_ds):
    if time_res_to_average[-1] == 'D':
        timestamp_start = ini_date[0:4]+'-'+ini_date[4:6]+'-'+ini_date[6:8]
        timestamp_end = fin_date[0:4]+'-'+fin_date[4:6]+'-'+fin_date[6:8]
        datetime_obj1 = datetime.strptime(timestamp_start, '%Y-%m-%d')
        datetime_obj2 = datetime.strptime(timestamp_end, '%Y-%m-%d')

    if time_res_to_average[-1] == 'M':
        timestamp_start = ini_date[0:4]+'-'+ini_date[4:6]
        timestamp_end = fin_date[0:4]+'-'+fin_date[4:6]
        datetime_obj1 = datetime.strptime(timestamp_start, '%Y-%m')
        datetime_obj2 = datetime.strptime(timestamp_end, '%Y-%m')
    print("TIME",averaged_ds['TIME'])
    print("datetime_obj1: ",datetime_obj1)
    print("datetime_obj2: ",datetime_obj2)
    closerval1 = averaged_ds['TIME'].sel(TIME=datetime_obj1,method="backfill")
    idx1 = averaged_ds['TIME'].astype(str).values.tolist().index(str(closerval1.data))
    print("nearest start month: ",averaged_ds['TIME'][idx1])
    closerval2 = averaged_ds['TIME'].sel(TIME=datetime_obj2,method="backfill")
    idx2 = averaged_ds['TIME'].astype(str).values.tolist().index(str(closerval2.data))
    print("nearest end month: ",averaged_ds['TIME'][idx2])
    return idx1,idx2

def Get_Max_Min_Interpolated_Model(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr):
    min_value = 0
    max_value = 0
    for time_counter,index in enumerate(range(idx1,idx2+1)):
        U = averaged_ds['EWCT'][index,0].data
        V = averaged_ds['NSCT'][index,0].data
        speed_hfr = (U ** 2 + V ** 2) ** 0.5
        mask_hfr=np.ma.masked_invalid(speed_hfr).mask

        subset_speed_model_instant = seaoverland(masked_subset_speed_model[time_counter],3)

        f=interpolate.interp2d(x_subset_model,y_subset_model,subset_speed_model_instant)
        speed_interpolated=f(lon_hfr,lat_hfr)
        masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_hfr)
        min_interpolated_subset_model = np.nanmin(masked_speed_interpolated.data)
        min_value = min(min_value,min_interpolated_subset_model)
        max_interpolated_subset_model = np.nanmax(masked_speed_interpolated.data)
        max_value = max(max_value,max_interpolated_subset_model)
    return min_value,max_value

def Get_Max_Min_Bias(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr):
    min_value = 0.0
    max_value = 0.0
    for time_counter,index in enumerate(range(idx1,idx2+1)):
        U = averaged_ds['EWCT'][index,0].data
        V = averaged_ds['NSCT'][index,0].data
        speed_hfr = (U ** 2 + V ** 2) ** 0.5
        mask_hfr=np.ma.masked_invalid(speed_hfr).mask

        subset_speed_model_instant = seaoverland(masked_subset_speed_model[time_counter],3)

        f=interpolate.interp2d(x_subset_model,y_subset_model,subset_speed_model_instant)
        speed_interpolated=f(lon_hfr,lat_hfr)
        masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_hfr)
        min_bias = np.nanmin(masked_speed_interpolated.data-speed_hfr.data)
        min_value = min(min_value,min_bias)
        max_bias = np.nanmax(masked_speed_interpolated.data-speed_hfr.data)
        max_value = max(max_value,max_bias)

    return min_value,max_value

def Get_Max_Min_Rmsd(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr):
    min_value = 0.0
    max_value = 0.0
    for time_counter,index in enumerate(range(idx1,idx2+1)):
        U = averaged_ds['EWCT'][index,0].data
        V = averaged_ds['NSCT'][index,0].data
        speed_hfr = (U ** 2 + V ** 2) ** 0.5
        mask_hfr=np.ma.masked_invalid(speed_hfr).mask

        subset_speed_model_instant = seaoverland(masked_subset_speed_model[time_counter],3)

        f=interpolate.interp2d(x_subset_model,y_subset_model,subset_speed_model_instant)
        speed_interpolated=f(lon_hfr,lat_hfr)
        masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_hfr)
        min_rmsd = np.nanmin(np.sqrt((masked_speed_interpolated.data-speed_hfr.data)**2))
        min_value = min(min_value,min_rmsd)
        max_rmsd = np.nanmax(np.sqrt((masked_speed_interpolated.data-speed_hfr.data)**2))
        max_value = max(max_value,max_rmsd)

    return min_value,max_value

def plot_hfr_wind_field(info, extent, min_hfr, min_model, max_hfr, max_model, x, y, speed_hfr, skip, skip_coords, date_str, ds, output_plot_folder):
    plt.figure(num=None, figsize=(18, 13), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.axes(projection=ccrs.Mercator())# Map projection
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--') #adding grid lines
    #plotting antennas
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for antenna in info['antennas']:
        plt.plot(antenna['lon'], antenna['lat'], color=np.random.rand(3,), markeredgecolor=np.random.rand(3,), marker='o',transform=ccrs.Geodetic(),label=antenna['name'])#add points
    if ds.id in ["GL_TV_HF_HFR-TirLig-Total","GL_TV_HF_HFR-NAdr-Total"]:
        ax.legend(loc='lower left')
    else:
        ax.legend()
    #personalized limts
    ax.set_extent(extent)
    # quiver plot: set vectors, colormaps, colorbars
    print("max: ", max(max_hfr,max_model))
    norm = colors.Normalize(vmin=min(min_hfr,min_model), vmax=max(max_hfr,max_model))
    ax.pcolor(x, y, speed_hfr, cmap='viridis', vmin=min(min_hfr,min_model), vmax=max(max_hfr,max_model),
            transform=cartopy.crs.PlateCarree())
    # quiver plot: arrows
    Q=ax.quiver(x[skip_coords], y[skip_coords], U[skip]/max(max_hfr,max_model), V[skip]/max(max_hfr,max_model),
            transform=cartopy.crs.PlateCarree())
    ax.quiverkey(Q, 0.1, 0.9, 0.4, r'$0.4 m/s$',fontproperties={'weight': 'bold'})
    # title and colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis',norm=norm)
    a = np.random.random((10, 20))
    im_ratio = a.shape[0]/a.shape[1]
    plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.15, fraction=0.047*im_ratio).ax.set_xlabel('velocity (m/s)', labelpad=9)
    plt.title(date_str+' | ' + ds.id + '\n surface current velocity', pad=28)
    figure_name = ds.id+'_surface_current_velocity_'+ date_str +'.png'
    plt.savefig(output_plot_folder+figure_name, dpi=300, bbox_inches = "tight")
    plt.clf()

def plot_model_wind_field(info, extent, min_hfr, min_model, max_hfr, max_model, x, y, skip, skip_coords, masked_speed_interpolated, masked_u_interpolated, masked_v_interpolated, date_str, output_plot_folder,title_substring,name_file_substring):

    plt.figure(num=None, figsize=(18, 13), dpi=100, facecolor='w', edgecolor='k')
    # Map projection
    ax = plt.axes(projection=ccrs.Mercator())
    # adding grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                  color='gray', alpha=0.5, linestyle='--')
    #plotting antennas
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for antenna in info['antennas']:
        plt.plot(antenna['lon'], antenna['lat'], color=np.random.rand(3,),
             markeredgecolor=np.random.rand(3,), marker='o',transform=ccrs.Geodetic(),label=antenna['name'])
    if ds.id in ["GL_TV_HF_HFR-TirLig-Total","GL_TV_HF_HFR-NAdr-Total"]:
        ax.legend(loc='lower left')
    else:
        ax.legend()
    # personalized limts
    ax.set_extent(extent)
    # quiver plot: set vectors, colormaps, colorbars
    print("max: ", max(max_hfr,max_model))
    norm = colors.Normalize(vmin=min(min_hfr,min_model), vmax=max(max_hfr,max_model))
    ax.pcolor(x, y, masked_speed_interpolated, cmap='viridis', vmin=min(min_hfr,min_model),
            vmax=max(max_hfr,max_model), transform=cartopy.crs.PlateCarree())
    # quiver plot: arrows
    Q=ax.quiver(x[skip_coords], y[skip_coords], masked_u_interpolated[skip]/max(max_hfr,max_model),
            masked_v_interpolated[skip]/max(max_hfr,max_model), transform=cartopy.crs.PlateCarree())
    ax.quiverkey(Q, 0.1, 0.9, 0.4, r'$0.4 m/s$',fontproperties={'weight': 'bold'})

    # title and colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis',norm=norm)
    a = np.random.random((10, 20))
    im_ratio = a.shape[0]/a.shape[1]
    plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.15, fraction=0.047*im_ratio).ax.set_xlabel('velocity (m/s)', labelpad=9)
    plt.title(date_str+' | ' + ds.id + '\n'+ title_substring, pad=28)
    figure_name = ds.id+name_file_substring+ date_str +'.png'
    plt.savefig(output_plot_folder+figure_name, dpi=300, bbox_inches = "tight")
    plt.clf()

def plot_interpolated_hfr_wind_field(info, extent, min_hfr, min_model_value, max_hfr, max_model_value, x_subset_model, y_subset_model, skip, skip_coords, masked_hfr_speed_interpolated, masked_hfr_u_interpolated, masked_hfr_v_interpolated, date_str, output_plot_folder,title_substring,name_file_substring):

    plt.figure(num=None, figsize=(18, 13), dpi=100, facecolor='w', edgecolor='k')
    # Map projection
    ax = plt.axes(projection=ccrs.Mercator())
    # adding grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                  color='gray', alpha=0.5, linestyle='--')
    #plotting antennas
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for antenna in info['antennas']:
        plt.plot(antenna['lon'], antenna['lat'], color=np.random.rand(3,),
             markeredgecolor=np.random.rand(3,), marker='o',transform=ccrs.Geodetic(),label=antenna['name'])
    if ds.id in ["GL_TV_HF_HFR-TirLig-Total","GL_TV_HF_HFR-NAdr-Total"]:
        ax.legend(loc='lower left')
    else:
        ax.legend()
    # personalized limts
    ax.set_extent(extent)
    # quiver plot: set vectors, colormaps, colorbars
    print("max: ", max(max_hfr,max_model))
    norm = colors.Normalize(vmin=min(min_hfr,min_model), vmax=max(max_hfr,max_model))
    ax.pcolor(x_subset_model, y_subset_model, masked_hfr_speed_interpolated, cmap='viridis', vmin=min(min_hfr,min_model),
            vmax=max(max_hfr,max_model), transform=cartopy.crs.PlateCarree())
    # quiver plot: arrows
    Q=ax.quiver(x_subset_model[skip_coords], y_subset_model[skip_coords], masked_hfr_u_interpolated[skip]/max(max_hfr,max_model),
            masked_hfr_v_interpolated[skip]/max(max_hfr,max_model), transform=cartopy.crs.PlateCarree())
    ax.quiverkey(Q, 0.1, 0.9, 0.4, r'$0.4 m/s$',fontproperties={'weight': 'bold'})

    # title and colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis',norm=norm)
    a = np.random.random((10, 20))
    im_ratio = a.shape[0]/a.shape[1]
    plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.15, fraction=0.047*im_ratio).ax.set_xlabel('velocity (m/s)', labelpad=9)
    plt.title(date_str+' | ' + ds.id + '\n' + title_substring, pad=28)
    figure_name = ds.id+name_file_substring+ date_str +'.png'
    plt.savefig(output_plot_folder+figure_name, dpi=300, bbox_inches = "tight")
    plt.clf()

def plot_bias(info, extent, x, y, min_bias, max_bias, masked_speed_interpolated, speed_hfr, date_str, output_plot_folder,title_substring,name_file_substring):

    plt.figure(num=None, figsize=(18, 13), dpi=100, facecolor='w', edgecolor='k')
    # Map projection
    ax = plt.axes(projection=ccrs.Mercator())
    # adding grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                  color='gray', alpha=0.5, linestyle='--')
    #plotting antennas
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for antenna in info['antennas']:
        plt.plot(antenna['lon'], antenna['lat'], color=np.random.rand(3,),
             markeredgecolor=np.random.rand(3,), marker='o',transform=ccrs.Geodetic(),label=antenna['name'])
    if ds.id in ["GL_TV_HF_HFR-TirLig-Total","GL_TV_HF_HFR-NAdr-Total"]:
        ax.legend(loc='lower left')
    else:
        ax.legend()    
    # personalized limts
    ax.set_extent(extent)
    # quiver plot: set vectors, colormaps, colorbars
    norm = colors.Normalize(vmin=min_bias, vmax=max_bias)
    ax.pcolor(x, y, masked_speed_interpolated-speed_hfr, cmap='RdBu_r', vmin=min_bias, vmax=max_bias, transform=cartopy.crs.PlateCarree())
    # quiver plot: arrows
    #Q=ax.quiver(x[skip_coords], y[skip_coords], masked_u_interpolated[skip],
    #        masked_v_interpolated[skip], transform=cartopy.crs.PlateCarree(), scale=10)
    #ax.quiverkey(Q, 0.1, 0.9, 0.4, r'$0.4 m/s$',fontproperties={'weight': 'bold'})

    # title and colorbar
    sm = plt.cm.ScalarMappable(cmap='RdBu_r',norm=norm)
    a = np.random.random((10, 20))
    im_ratio = a.shape[0]/a.shape[1]
    plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.15, fraction=0.047*im_ratio).ax.set_xlabel('velocity (m/s)', labelpad=9)
    plt.title(date_str+' | ' + ds.id + '\n '+ title_substring, pad=28)
    figure_name = ds.id+ name_file_substring + date_str +'.png'
    plt.savefig(output_plot_folder+figure_name, dpi=300, bbox_inches = "tight")
    plt.clf()

def plot_rmsd(info, extent, x, y, min_rmsd, max_rmsd, masked_speed_interpolated, speed_hfr, date_str, output_plot_folder,title_substring,name_file_substring):

    plt.figure(num=None, figsize=(18, 13), dpi=100, facecolor='w', edgecolor='k')
    # Map projection
    ax = plt.axes(projection=ccrs.Mercator())
    # adding grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                  color='gray', alpha=0.5, linestyle='--')
    #plotting antennas
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    for antenna in info['antennas']:
        plt.plot(antenna['lon'], antenna['lat'], color=np.random.rand(3,),
             markeredgecolor=np.random.rand(3,), marker='o',transform=ccrs.Geodetic(),label=antenna['name'])
    if ds.id in ["GL_TV_HF_HFR-TirLig-Total","GL_TV_HF_HFR-NAdr-Total"]:
        ax.legend(loc='lower left')
    else:
        ax.legend()
    # personalized limts
    ax.set_extent(extent)
    # quiver plot: set vectors, colormaps, colorbars
    norm = colors.Normalize(vmin=min_rmsd, vmax=max_rmsd)
    ax.pcolor(x, y, np.sqrt((masked_speed_interpolated-speed_hfr)**2), cmap='viridis', vmin=min_rmsd, vmax=max_rmsd, transform=cartopy.crs.PlateCarree())
    # quiver plot: arrows
    #Q=ax.quiver(x[skip_coords], y[skip_coords], masked_u_interpolated[skip],
    #        masked_v_interpolated[skip], transform=cartopy.crs.PlateCarree(), scale=10)
    #ax.quiverkey(Q, 0.1, 0.9, 0.4, r'$0.4 m/s$',fontproperties={'weight': 'bold'})

    # title and colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis',norm=norm)
    a = np.random.random((10, 20))
    im_ratio = a.shape[0]/a.shape[1]
    plt.colorbar(sm,ax=ax, orientation='vertical', pad=0.15, fraction=0.047*im_ratio).ax.set_xlabel('velocity (m/s)', labelpad=9)
    plt.title(date_str+' | ' + ds.id + '\n ' + title_substring, pad=28)
    figure_name = ds.id+ name_file_substring + date_str +'.png'
    plt.savefig(output_plot_folder+figure_name, dpi=300, bbox_inches = "tight")
    plt.clf()

def plot_mod_obs_ts_comparison(obs_ts, mod_ts, time_res_to_average, ds, date_in, date_fin, output_plot_folder,timerange,name_exp,title_substring,name_file_substring):
    plotname = ds.id + '_' + date_in + '_' + date_fin + '_' + time_res_to_average + title_substring +'.png'
    fig = plt.figure(figsize=(18,12))
    ax = fig.add_subplot(111)
    plt.rc('font', size=24)
    plt.title(title_substring+': '+ ds.id + '\n Period: '+ date_in + '-' + date_fin, fontsize=29)
    
    mean_vel_mod = round(np.nanmean(np.array(mod_ts)),2)
    print("timerange shape: ",timerange.shape)
    print("mod_ts shape: ", np.array(mod_ts).shape)
    plt.plot(timerange,np.array(mod_ts),label = name_exp + ' : '+str(mean_vel_mod)+' m/s', linewidth=2)
    mean_vel_obs = round(np.nanmean(np.array(obs_ts)),2)
    plt.plot(timerange,np.array(obs_ts),label = 'Observation : '+str(mean_vel_obs)+' m/s', linewidth=2)
    plt.grid()
    ax.tick_params(axis='both', labelsize=26)
    if time_res_to_average[1]=='D':
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=int(time_res_to_average[0])))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
    if time_res_to_average[1]=='M':
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=int(time_res_to_average[0])))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%Y'))
    if time_res_to_average[1]=='Y':
        ax.xaxis.set_major_locator(mdates.YearLocator(interval=int(time_res_to_average[0])))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    fig.autofmt_xdate()
    plt.ylabel('Velocity [m/s]', fontsize=40)
    plt.xlabel('Date', fontsize=40)
    plt.legend(prop={'size': 30}, framealpha=0.2)
    plt.savefig(output_plot_folder + plotname)
    plt.clf()

    return mean_vel_mod,mean_vel_obs

class MarkerHandler(HandlerBase):
    def create_artists(self, legend, tup, xdescent, ydescent,
                        width, height, fontsize,trans):
        return [plt.Line2D([width/2], [height/2.],ls="",
                       marker=tup,color='blue', transform=trans)]

def mscatter(x,y,ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    ax = ax or plt.gca()
    sc = ax.scatter(x,y,clip_on=False,cmap='plasma',**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

def scatterPlot(mod, obs, outname, name, n_stations, n_time, possible_markers, hfr_name, pos_colors, time_string,**kwargs):

    if np.isnan(obs).any() or np.isnan(mod).any():

        obs_no_nan = obs[~np.isnan(obs) & ~np.isnan(mod)]
        mod_no_nan = mod[~np.isnan(obs) & ~np.isnan(mod)]
        xy = np.vstack([obs_no_nan, mod_no_nan])
    else:
        xy = np.vstack([obs, mod])

    #possible_markers=np.array(["o","^","s","P","*","D"])
    if n_stations==1:
        print("prima repeat: ",possible_markers)
        m=np.repeat(possible_markers,n_time)
    if n_stations>1:
        m=np.repeat(possible_markers[:n_stations],n_time)
    c_prova = np.tile(np.arange(0,6*n_time,6),n_stations)
#    z = gaussian_kde(xy)(xy)
#    idx = z.argsort()

    if np.isnan(obs).any() or np.isnan(mod).any():
#        x, y, z = obs_no_nan[idx], mod_no_nan[idx], z[idx]
        x, y = obs_no_nan, mod_no_nan
    else:
        #x, y, z = obs[idx], mod[idx], z[idx]
        x, y = obs, mod

    color_list = pos_colors
    color_list_seq = np.tile(color_list[:n_time],n_stations)
    classes = time_string
    markers_labels = hfr_name
    fig, ax = plt.subplots(figsize=(10,6))
#    im = ax.scatter(x, y, c=z, s=8, edgecolor=None, cmap='jet', clip_on=False)
    #im = mscatter(x, y, cmap_prova(c_prova), ax=ax, m=m, c=c_prova,s=15)
    im = mscatter(x, y, ax=ax, m=m, c=c_prova,s=15)
    marker_array=[]
    if n_stations==1:
        marker_array.append(mlines.Line2D([], [], color='blue', marker=possible_markers, linestyle='None', markersize=5, label=markers_labels))
    else:
        for mark,mark_label in zip(possible_markers,markers_labels):
            print('label: ',mark_label)
            marker_array.append(mlines.Line2D([], [], color='blue', marker=mark, linestyle='None', markersize=5, label=mark_label))

    legend_1=plt.legend(handles=im.legend_elements(num=n_time)[0], labels=classes, loc='right',prop={"size":9},bbox_to_anchor=(1.3, 0.5))
    plt.legend(handles=marker_array,loc='upper left',prop={"size":12})
    plt.gca().add_artist(legend_1)
    #plt.legend(list(tuple(possible_markers[:n_stations])), markers_labels,handler_map={tuple:MarkerHandler()},loc='upper left') 
    #plt.legend((im),(hfr_name),loc='lower right')

    maxVal = np.nanmax((x, y))
    ax.set_ylim(0, maxVal)
    ax.set_xlim(0, maxVal)
    ax.set_aspect(1.0)
    ax.tick_params(axis='both', labelsize=12.5)

    bias = BIAS(y,x)
    corr, _ = pearsonr(x, y)
    rmse=RMSE(y,x)
    nstd=Normalized_std(y,x)
    si=ScatterIndex(y,x)
    slope,intercept, rvalue,pvalue,stderr=linregress(y,x)

    prova = x[:,np.newaxis]
    a, _, _, _ = np.linalg.lstsq(prova, y)
    xseq = np.linspace(0, maxVal, num=100)
    ax.plot(xseq, a*xseq, 'r-')
    plt.text(0.01, 0.7, name, weight='bold',transform=plt.gcf().transFigure,fontsize=16)

    plt.text(0.01, 0.32, 'Entries: %s\n'
             'BIAS: %s m/s\n'
             'RMSD: %s m/s\n'
             'NSTD: %s\n'
             'SI: %s\n'
             'corr:%s\n'
             'Slope: %s\n'
             'STDerr: %s m/s'
             %(len(obs),bias,rmse,nstd,si,np.round(corr,2),
               np.round(a[0],2),np.round(stderr,2)),transform=plt.gcf().transFigure,fontsize=15)

    stat_array=[bias,rmse,si,np.round(corr,2),np.round(stderr,2),len(obs)]

    if 'title' in kwargs:
        plt.title(kwargs['title'], fontsize=15, x=0.5, y=1.01)

    if 'xlabel' in kwargs:
        plt.xlabel(kwargs['xlabel'], fontsize=18)

    if 'ylabel' in kwargs:
        plt.ylabel(kwargs['ylabel'], fontsize=18)


    ax.plot([0,maxVal],[0,maxVal],c='k',linestyle='-.')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #ticks_1 = np.linspace(z.min(), z.max(), 5,endpoint=True)
    #cbar=plt.colorbar(im,fraction=0.02,ticks=ticks_1)
    #cbar.ax.set_yticklabels(['{:.1f}'.format(x) for x in ticks_1], fontsize=13)
    #cbar.set_label('probaility density [%]', rotation=270,size=18,labelpad=15)

    plt.savefig(outname)
    plt.close()
    return stat_array

def interp_hfr_mask_to_mod_mask(x_obs,y_obs,mask_hfr,x_model,y_model,threshold):

    f=interpolate.interp2d(x_obs,y_obs,mask_hfr)
    hfr_mask_interpolated=f(x_model,y_model)
    hfr_mask_interpolated[hfr_mask_interpolated<threshold] = 0
    hfr_mask_interpolated[hfr_mask_interpolated>threshold] = 1
    return hfr_mask_interpolated

def interp_mod_to_obs(x_mod,y_mod,speed_model,u_model,v_model,lon_obs,lat_obs,mask_obs):

    f=interpolate.interp2d(x_mod,y_mod,speed_model)
    speed_interpolated=f(lon_obs,lat_obs)
    masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_obs)
    spatial_mean_model_ts = masked_speed_interpolated.mean()

    f=interpolate.interp2d(x_mod,y_mod,u_model)
    u_interpolated=f(lon_obs,lat_obs)
    masked_u_interpolated = ma.masked_array(u_interpolated, mask=mask_obs)

    f=interpolate.interp2d(x_mod,y_mod,v_model)
    v_interpolated=f(lon_obs,lat_obs)
    masked_v_interpolated = ma.masked_array(v_interpolated, mask=mask_obs)

    return masked_speed_interpolated, masked_u_interpolated, masked_v_interpolated, spatial_mean_model_ts

def interp_obs_to_mod(lon_obs,lat_obs,sol_speed_hfr,sol_u_hfr,sol_v_hfr,x_model,y_model,mask_mod):

    f=interpolate.interp2d(lon_obs,lat_obs,sol_speed_hfr)
    speed_interpolated=f(x_model,y_model)
    masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_mod)
    spatial_mean_hfr_ts = masked_speed_interpolated.mean()

    f=interpolate.interp2d(lon_obs,lat_obs,sol_u_hfr)
    u_interpolated=f(x_model,y_model)
    masked_u_interpolated = ma.masked_array(u_interpolated, mask=mask_mod)

    f=interpolate.interp2d(lon_obs,lat_obs,sol_v_hfr)
    v_interpolated=f(x_model,y_model)
    masked_v_interpolated = ma.masked_array(v_interpolated, mask=mask_mod)

    return masked_speed_interpolated, masked_u_interpolated, masked_v_interpolated, spatial_mean_hfr_ts

if __name__ == "__main__":

    argv=sys.argv
    print("argv: ", argv)
    ini_date=argv[1]
    fin_date=argv[2]
    path_to_hfr_files=argv[3]
    path_to_mesh_mask=argv[4]
    time_res_to_average=argv[5]
    path_to_model_files=argv[6]
    path_to_destag_model_files=argv[7]
    want_destag=argv[8]
    time_res_model=argv[9]
    name_exp=argv[10]
    output_plot_folder=argv[11]
    destaggered_u_combined=argv[12]
    destaggered_v_combined=argv[13]
    work_dir=argv[14]
#ini_date = "20210801"
#fin_date = "20210831"
#path_to_hfr_files = "/work/oda/mg28621/prova_destag/hfr_validation/hfr_data/prova_hfr/"
#path_to_mesh_mask = "/work/oda/mg28621/prova_destag/tool_hfr/2021-sicomar-summer-school-master/mesh_mask.nc"
#time_res_to_average = "1D"
#path_to_model_files = "/work/oda/med_dev/exp/eas6_v8/assw_eas7_v0/output/"
#path_to_destag_model_files = "/work/oda/mg28621/prova_destag/tool_hfr/2021-sicomar-summer-school-master/destag_output_folder/"
#want_destag = "no"
#time_res_model = "1h"
#name_exp = "assw_eas7_v0"

    os.makedirs(path_to_destag_model_files, exist_ok=True)
    print("The new directory is created!")

    os.makedirs(output_plot_folder, exist_ok=True)
    print("The new directory is created!")

    listOfFiles = Get_List_Of_Files(path_to_hfr_files)

    #if want_destag == 'yes':
    #    print("start destaggering")
    #    destaggering(ini_date, fin_date, path_to_model_files, path_to_destag_model_files, 
    #             name_exp, time_res_model, path_to_mesh_mask)

    start_date = date(int(ini_date[0:4]),int(ini_date[4:6]) , int(ini_date[6:8]))
    end_date = date(int(fin_date[0:4]),int(fin_date[4:6]) , int(fin_date[6:8]))
    if time_res_to_average[1]=='M':
        timerange = pd.date_range(start_date, end_date, freq=time_res_to_average[1]) - pd.DateOffset(days=15)
    if time_res_to_average[1]=='D':
        timerange = pd.date_range(start_date, end_date, freq=time_res_to_average[1]) - pd.DateOffset(hours=12)
    mesh_mask = xarray.open_dataset(path_to_mesh_mask)
    u_mask = mesh_mask.umask.values
    v_mask = mesh_mask.vmask.values
    t_mask = mesh_mask.tmask.values

    # read the model input
#    ds_model_u = xarray.open_dataset('/work/oda/mg28621/prova_destag/tool_hfr/2021-sicomar-summer-school-master/right_model_U_combined.nc')
    ds_model_u = xarray.open_dataset(destaggered_u_combined)
#    ds_model_v = xarray.open_dataset('/work/oda/mg28621/prova_destag/tool_hfr/2021-sicomar-summer-school-master/right_model_V_combined.nc')
    ds_model_v = xarray.open_dataset(destaggered_v_combined)

    # average in time the model input
    averaged_model_u = ds_model_u.resample(time_counter=time_res_to_average).mean(skipna=True)
    ds_model_u.close()
    averaged_model_v = ds_model_v.resample(time_counter=time_res_to_average).mean(skipna=True)
    ds_model_v.close()

    x_model = averaged_model_u['nav_lon'].data[0,:]
    y_model = averaged_model_u['nav_lat'].data[:,0]

    u_model = averaged_model_u.variables['destaggered_u'][:,:,:].data
    v_model = averaged_model_v.variables['destaggered_v'][:,:,:].data
    speed_model = np.sqrt(u_model*u_model + v_model*v_model)

    # repeat the mask for all the days
    print("time speed_model shape: ",speed_model.shape[0])
    t_mask = np.repeat(t_mask[:, :, :, :], speed_model.shape[0], axis=0)
    t_mask = t_mask[:, 0, :, :]
    print("t_mask shape: ",t_mask.shape)
    print("t_mask: ",t_mask)

    speed_model = ma.masked_array(speed_model, mask=np.logical_not(t_mask))
    min_model = np.nanmin(speed_model.data)
    max_model = np.nanmax(speed_model.data)
    u_model = ma.masked_array(u_model, mask=np.logical_not(t_mask))
    v_model = ma.masked_array(v_model, mask=np.logical_not(t_mask))

    string_time_res = Get_String_Time_Resolution(start_date, end_date, time_res_to_average)

    skip = (slice(None, None, 3), slice(None, None, 3))
    skip_coords = (slice(None,None,3))

    skip_model = (slice(None,None, 1), slice(None,None, 1))
    skip_coords_model = (slice(None,None,1))

    statistics={}
    statistics_rev_interp={}
    mod_array = np.array([])
    obs_array = np.array([])
    mod_array_rev_interp = np.array([])
    obs_array_rev_interp = np.array([])

    possible_colors=['red', 'blue', 'black','green','purple','orange','brown','pink','grey','olive']
    possible_markers=np.array(["o","^","s","P","*","D"])
    hfr_names=[]
    hfr_names_rev=[]
    for hfr_counter,hfr_file in enumerate(listOfFiles):

        print('loading ' + hfr_file + ' ...')
        ds = xarray.open_dataset(hfr_file)
        #lat_hfr=ds.variables['LATITUDE'][:]
        #lon_hfr=ds.variables['LONGITUDE'][:]
        info = getSourceAntennas(ds)
        ds_1=ds[['QCflag','EWCT','NSCT','LATITUDE','LONGITUDE']]
        print("shape ds_1['EWCT']: ",ds_1['EWCT'].data.shape)
        ds_restricted=ds_1[['EWCT','NSCT','LATITUDE','LONGITUDE']].where((ds.QCflag==0) | (ds.QCflag == 1) | (ds.QCflag==2),drop=True)
        print("shape ds_restricted['EWCT']: ",ds_restricted['EWCT'].data.shape)
        averaged_ds=ds_restricted.resample(TIME=time_res_to_average).mean(skipna=True)
        print("shape averaged_ds['EWCT']: ",averaged_ds['EWCT'].data.shape)
        lat_hfr=averaged_ds.variables['LATITUDE'][:]
        lon_hfr=averaged_ds.variables['LONGITUDE'][:]

        x = averaged_ds['LONGITUDE'].data
        y = averaged_ds['LATITUDE'].data

        timestamp_start = ini_date[0:4]+'-'+ini_date[4:6]+'-'+ini_date[6:8]
        timestamp_end = fin_date[0:4]+'-'+fin_date[4:6]+'-'+fin_date[6:8]

        idx1,idx2 = Get_Closest_Hfr_Time_Range_Index(time_res_to_average,ini_date,fin_date,averaged_ds)

        extent = [info['bbox'][0], info['bbox'][1]+0.2, 
        info['bbox'][2], info['bbox'][3]+0.1] 

        max_hfr = np.nanmax((averaged_ds['EWCT'][idx1:idx2+1,0].data ** 2 + averaged_ds['NSCT'][idx1:idx2+1,0].data ** 2) ** 0.5)
        min_hfr = np.nanmin((averaged_ds['EWCT'][idx1:idx2+1,0].data ** 2 + averaged_ds['NSCT'][idx1:idx2+1,0].data ** 2) ** 0.5)

        lon_model = averaged_model_u['nav_lon'].data[0,:]
        lat_model = averaged_model_u['nav_lat'].data[:,0]

        closer_min_mod_lon = find_nearest(lon_model, np.nanmin(ds['LONGITUDE'].data))
        closer_min_mod_lat = find_nearest(lat_model, np.nanmin(ds['LATITUDE'].data))
        closer_max_mod_lon = find_nearest(lon_model, np.nanmax(ds['LONGITUDE'].data))
        closer_max_mod_lat = find_nearest(lat_model, np.nanmax(ds['LATITUDE'].data))

        coord_min_lon_min_lat = [closer_min_mod_lon, closer_min_mod_lat]
        coord_idx_min_lon_min_lat = np.argwhere((averaged_model_u['nav_lon'].data==coord_min_lon_min_lat[0]) & 
                                        (averaged_model_u['nav_lat'].data==coord_min_lon_min_lat[1]))[0]
        coord_min_lon_max_lat = [closer_min_mod_lon, closer_max_mod_lat]
        coord_idx_min_lon_max_lat = np.argwhere((averaged_model_u['nav_lon'].data==coord_min_lon_max_lat[0]) & 
                                        (averaged_model_u['nav_lat'].data==coord_min_lon_max_lat[1]))[0]
        coord_max_lon_min_lat = [closer_max_mod_lon, closer_min_mod_lat]
        coord_idx_max_lon_min_lat = np.argwhere((averaged_model_u['nav_lon'].data==coord_max_lon_min_lat[0]) & 
                                        (averaged_model_u['nav_lat'].data==coord_max_lon_min_lat[1]))[0]
        coord_max_lon_max_lat = [closer_max_mod_lon, closer_max_mod_lat]
        coord_idx_max_lon_max_lat = np.argwhere((averaged_model_u['nav_lon'].data==coord_max_lon_max_lat[0]) & 
                                        (averaged_model_u['nav_lat'].data==coord_max_lon_max_lat[1]))[0]

        subset_averaged_model_u = averaged_model_u.isel(x=slice(coord_idx_min_lon_min_lat[1]-1,
								coord_idx_max_lon_min_lat[1]+1),
							y=slice(coord_idx_min_lon_min_lat[0]-1,
								coord_idx_min_lon_max_lat[0]+1))
        subset_averaged_model_v = averaged_model_v.isel(x=slice(coord_idx_min_lon_min_lat[1]-1,
								coord_idx_max_lon_min_lat[1]+1),
							y=slice(coord_idx_min_lon_min_lat[0]-1,
								coord_idx_min_lon_max_lat[0]+1))
        x_subset_model = subset_averaged_model_u['nav_lon'].data[0,:]
        y_subset_model = subset_averaged_model_u['nav_lat'].data[:,0]

        subset_u_model = subset_averaged_model_u.variables['destaggered_u'][:,:,:].data
        subset_v_model = subset_averaged_model_v.variables['destaggered_v'][:,:,:].data
        subset_speed_model = np.sqrt(subset_u_model * subset_u_model + subset_v_model * subset_v_model)

        subset_t_mask = t_mask[:, slice(coord_idx_min_lon_min_lat[0]-1, coord_idx_min_lon_max_lat[0]+1),slice(coord_idx_min_lon_min_lat[1]-1, coord_idx_max_lon_min_lat[1]+1)]
        subset_t_mask = np.logical_not(subset_t_mask)

        masked_subset_speed_model = ma.masked_array(subset_speed_model, mask=subset_t_mask)
        min_subset_model = np.nanmin(masked_subset_speed_model.data)
        max_subset_model = np.nanmax(masked_subset_speed_model.data)
        masked_subset_u_model = ma.masked_array(subset_u_model, mask=subset_t_mask)
        masked_subset_v_model = ma.masked_array(subset_v_model, mask=subset_t_mask)
        print(string_time_res)
        print(len(range(idx1,idx2+1)))
        min_model_value,max_model_value = Get_Max_Min_Interpolated_Model(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr)
        min_bias,max_bias = Get_Max_Min_Bias(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr)
        min_rmsd,max_rmsd = Get_Max_Min_Rmsd(idx1,idx2,averaged_ds,masked_subset_speed_model,x_subset_model,y_subset_model,lon_hfr,lat_hfr)
        spatial_mean_hfr_ts=[]
        spatial_mean_model_ts=[]
        spatial_not_interp_mean_model_ts=[]
        spatial_mean_interp_hfr_ts=[]
        for time_counter,index in enumerate(range(idx1,idx2+1)):
            date_str = string_time_res[time_counter]

            U = averaged_ds['EWCT'][index,0].data
            V = averaged_ds['NSCT'][index,0].data
            speed_hfr = (U ** 2 + V ** 2) ** 0.5
            #print("obs type: ", speed_hfr)
            #print("speed_hfr.mean: ", speed_hfr.mean())
            spatial_mean_hfr_ts.append(np.nanmean(speed_hfr))
            
            print("mean_obs: ",spatial_mean_hfr_ts[time_counter])
            mask_hfr=np.ma.masked_invalid(speed_hfr).mask
            plot_hfr_wind_field(info, extent, min_hfr, min_model_value, max_hfr, max_model_value, x, y, speed_hfr, skip, skip_coords, date_str, ds, output_plot_folder)

            #u_model_instant = u_model[time_counter,:,:]
            #v_model_instant = v_model[time_counter,:,:]

            subset_speed_model_instant = seaoverland(masked_subset_speed_model[time_counter],3)
            subset_u_model_instant=seaoverland(masked_subset_u_model[time_counter],3)
            subset_v_model_instant=seaoverland(masked_subset_v_model[time_counter],3)
            spatial_not_interp_mean_model_ts.append(masked_subset_speed_model[time_counter].mean())

            masked_speed_interpolated, masked_u_interpolated, masked_v_interpolated, spatial_mean_model_ts_instant=interp_mod_to_obs(x_subset_model,y_subset_model,subset_speed_model_instant,subset_u_model_instant,subset_v_model_instant,lon_hfr,lat_hfr,mask_hfr)
            spatial_mean_model_ts.append(spatial_mean_model_ts_instant)

            masked_speed_hfr=ma.masked_array(speed_hfr,mask=mask_hfr)
            masked_U=ma.masked_array(U,mask=mask_hfr)
            masked_V=ma.masked_array(V,mask=mask_hfr)
            sol_speed_hfr = seaoverland(masked_speed_hfr,3)
            sol_u_hfr=seaoverland(masked_U,3)
            sol_v_hfr=seaoverland(masked_V,3)
            threshold=0.7
            step_lon=lon_hfr[1]-lon_hfr[0]
            step_lat=lat_hfr[1]-lat_hfr[0]
            X=np.concatenate(([lon_hfr[0]-step_lon], lon_hfr, [lon_hfr[-1]+step_lon]))
            Y=np.concatenate(([lat_hfr[0]-step_lat], lat_hfr, [lat_hfr[-1]+step_lat]))
            mask_hfr_prova=np.pad(np.logical_not(mask_hfr), 1)
            hfr_mask_interpolated=interp_hfr_mask_to_mod_mask(X,Y,np.logical_not(mask_hfr_prova),x_subset_model,y_subset_model,threshold)

            masked_hfr_speed_interpolated, masked_hfr_u_interpolated, masked_hfr_v_interpolated, spatial_mean_hfr_ts_instant=interp_obs_to_mod(lon_hfr,lat_hfr,sol_speed_hfr,sol_u_hfr,sol_v_hfr,x_subset_model,y_subset_model,hfr_mask_interpolated)
            spatial_mean_interp_hfr_ts.append(spatial_mean_hfr_ts_instant)

#            f=interpolate.interp2d(x_subset_model,y_subset_model,subset_speed_model_instant)
#            speed_interpolated=f(lon_hfr,lat_hfr)
#            masked_speed_interpolated = ma.masked_array(speed_interpolated, mask=mask_hfr)
#            spatial_mean_model_ts.append(masked_speed_interpolated.mean())
#            print("mean_mod: ",spatial_mean_model_ts[time_counter])
#            #min_interpolated_subset_model = np.nanmin(masked_speed_interpolated.data)
#            #max_interpolated_subset_model = np.nanmax(masked_speed_interpolated.data)
#            #print("min_interpolated_subset_model: ", min_interpolated_subset_model)
#            #print("max_interpolated_subset_model: ",  max_interpolated_subset_model)

#            f=interpolate.interp2d(x_subset_model,y_subset_model,subset_u_model_instant)
#            u_interpolated=f(lon_hfr,lat_hfr)
#            masked_u_interpolated = ma.masked_array(u_interpolated, mask=mask_hfr)

#            f=interpolate.interp2d(x_subset_model,y_subset_model,subset_v_model_instant)
#            v_interpolated=f(lon_hfr,lat_hfr)
#            masked_v_interpolated = ma.masked_array(v_interpolated, mask=mask_hfr)
            title_substring='interpolated model surface current'
            name_file_substring='model_surface_current_velocity_'
            plot_model_wind_field(info, extent, min_hfr, min_model_value, max_hfr, max_model_value, x, y, skip, skip_coords, masked_speed_interpolated, masked_u_interpolated, masked_v_interpolated, date_str, output_plot_folder,title_substring,name_file_substring)
            title_substring='surface current bias'
            name_file_substring='surface_current_velocity_bias'
            plot_bias(info, extent, x, y, min_bias, max_bias, masked_speed_interpolated, speed_hfr, date_str, output_plot_folder,title_substring,name_file_substring)
            title_substring='surface current rmsd'
            name_file_substring='surface_current_velocity_rmsd'
            plot_rmsd(info, extent, x, y, min_rmsd, max_rmsd, masked_speed_interpolated, speed_hfr, date_str, output_plot_folder,title_substring,name_file_substring)

            title_substring='interpolated hfr surface current'
            name_file_substring='interp_hfr_surface_current_velocity_'
            plot_interpolated_hfr_wind_field(info, extent, min_hfr, min_model_value, max_hfr, max_model_value, x_subset_model, y_subset_model, skip_model, skip_coords_model, masked_hfr_speed_interpolated, masked_hfr_u_interpolated, masked_hfr_v_interpolated, date_str, output_plot_folder,title_substring,name_file_substring)
            title_substring='model surface current'
            name_file_substring='model_surface_current_velocity_'
            plot_model_wind_field(info, extent, min_hfr, min_model_value, max_hfr, max_model_value, x_subset_model, y_subset_model, skip_model, skip_coords_model, masked_subset_speed_model[time_counter], masked_subset_u_model[time_counter], masked_subset_v_model[time_counter], date_str, output_plot_folder,title_substring,name_file_substring)
            title_substring='surface current bias (rev_interp)'
            name_file_substring='surface_current_velocity_bias_rev_interp'
            plot_bias(info, extent, x_subset_model, y_subset_model, min_bias, max_bias, masked_subset_speed_model[time_counter], masked_hfr_speed_interpolated, date_str, output_plot_folder,title_substring,name_file_substring)
            title_substring='surface current rmsd (rev_interp)'
            name_file_substring='surface_current_velocity_rmsd_rev_interp'
            plot_rmsd(info, extent, x_subset_model, y_subset_model, min_rmsd, max_rmsd, masked_subset_speed_model[time_counter], masked_hfr_speed_interpolated, date_str, output_plot_folder,title_substring,name_file_substring)


        print("spatial mean hfr ts: ", spatial_mean_hfr_ts)
        print("spatial mean model ts: ",spatial_mean_model_ts)
        print("sua shape: ",len(spatial_mean_model_ts))
        title_substring='Spatial Surface Current Velocity Mean Comparison'
        name_file_substring='_mod_obs_ts_comparison'
        mean_vel_mod,mean_vel_obs=plot_mod_obs_ts_comparison(spatial_mean_hfr_ts, spatial_mean_model_ts, time_res_to_average, ds, ini_date, fin_date, output_plot_folder,timerange,name_exp,title_substring,name_file_substring)
        tot_mean_stat=[mean_vel_mod,mean_vel_obs]
        title_substring='Spatial Surface Current Velocity Mean Comparison\n (Rev Interp)'
        name_file_substring='_mod_obs_ts_comparison_rev_interp'
        mean_vel_no_interp_mod,mean_vel_interp_obs=plot_mod_obs_ts_comparison(spatial_mean_interp_hfr_ts, spatial_not_interp_mean_model_ts, time_res_to_average, ds, ini_date, fin_date, output_plot_folder,timerange,name_exp,title_substring,name_file_substring)
        tot_mean_stat_rev_interp=[mean_vel_no_interp_mod,mean_vel_interp_obs]
        
        plotname = ds.id + '_' + ini_date + '_' + fin_date + '_' + time_res_to_average +  '_qqPlot.png'
        title = 'Spatial Mean Surface Current Velocity ' + ds.id + '\n Period: ' + ini_date + ' - ' + fin_date
        xlabel = 'Observation Current Velocity [m/s]'
        ylabel = 'Model Current Velocity [m/s]'
        print("spatial_mean_model_ts: ",np.array(spatial_mean_model_ts))
        print("spatial_mean_hfr_ts: ",np.array(spatial_mean_hfr_ts))
        splitted_name=ds.id.split("-")
        print(splitted_name)
        hfr_names.append(splitted_name[1])
        if timerange.shape[0] > 2:
            statistics_array=scatterPlot(np.array(spatial_mean_model_ts),np.array(spatial_mean_hfr_ts),output_plot_folder + plotname,name_exp,1,len(spatial_mean_model_ts),possible_markers[hfr_counter],splitted_name[1],possible_colors,string_time_res,title=title,xlabel=xlabel,ylabel=ylabel)
            row_stat = tot_mean_stat + statistics_array
            statistics[ds.id] = row_stat
        
        mod_array = np.concatenate([mod_array, np.array(spatial_mean_model_ts)])
        obs_array = np.concatenate([obs_array, np.array(spatial_mean_hfr_ts)])


        plotname = ds.id + '_' + ini_date + '_' + fin_date + '_' + time_res_to_average +  '_qqPlot_rev_interp.png'
        title = 'Spatial Mean Surface Current Velocity (Rev Interp)' + ds.id + '\n Period: ' + ini_date + ' - ' + fin_date
        xlabel = 'Observation Current Velocity [m/s]'
        ylabel = 'Model Current Velocity [m/s]'
        print("spatial_not_interp mean_model_ts: ",np.array(spatial_not_interp_mean_model_ts))
        print("spatial_mean_interp_hfr_ts: ",np.array(spatial_mean_interp_hfr_ts))
        splitted_name=ds.id.split("-")
        print(splitted_name)
        hfr_names_rev.append(splitted_name[1])
        if timerange.shape[0] > 2:
            statistics_array_rev_interp=scatterPlot(np.array(spatial_not_interp_mean_model_ts),np.array(spatial_mean_interp_hfr_ts),output_plot_folder + plotname,name_exp,1,len(spatial_not_interp_mean_model_ts),possible_markers[hfr_counter],splitted_name[1],possible_colors,string_time_res,title=title,xlabel=xlabel,ylabel=ylabel)
            row_stat_rev_interp = tot_mean_stat_rev_interp + statistics_array_rev_interp
            statistics_rev_interp[ds.id] = row_stat_rev_interp

        mod_array_rev_interp = np.concatenate([mod_array_rev_interp, np.array(spatial_not_interp_mean_model_ts)])
        obs_array_rev_interp = np.concatenate([obs_array_rev_interp, np.array(spatial_mean_interp_hfr_ts)])

    tot_mean_mod=round(np.nanmean(mod_array),2)
    tot_mean_obs=round(np.nanmean(obs_array),2)
    mean_all=[tot_mean_mod,tot_mean_obs]

    plotname = ini_date + '_' + fin_date + '_' + time_res_to_average +  '_qqPlot.png'
    title = 'Surface Current Velocity -ALL \n Period: ' + ini_date + '-' + fin_date
    xlabel = 'Observation Current Velocity [m/s]'
    ylabel = 'Model Current Velocity [m/s]'
    if timerange.shape[0] > 2:
        statistics_array=scatterPlot(mod_array,obs_array,output_plot_folder + plotname,name_exp,len(listOfFiles),timerange.shape[0],possible_markers,hfr_names,possible_colors,string_time_res,title=title,xlabel=xlabel,ylabel=ylabel)
        row_all = mean_all + statistics_array
        statistics["ALL HFR STATIONS"] = row_all

        a_file = open(work_dir+"statistics_" + name_exp + "_" + ini_date + "_" + fin_date + ".csv", "w")
        writer = csv.writer(a_file)
        writer.writerow(["name_hfr", "mean_mod", "mean_obs", "bias","rmse","si","corr","stderr","number_of_points"])
        for key, value in statistics.items():
            array = [key] + value
            print(array)
            writer.writerow(array)
        a_file.close()

    tot_mean_mod_rev_interp=round(np.nanmean(mod_array_rev_interp),2)
    tot_mean_obs_rev_interp=round(np.nanmean(obs_array_rev_interp),2)
    mean_all_rev_interp=[tot_mean_mod_rev_interp,tot_mean_obs_rev_interp]

    plotname = ini_date + '_' + fin_date + '_' + time_res_to_average +  '_qqPlot_rev_interp.png'
    title = 'Surface Current Velocity -ALL (Rev Interp)\n Period: ' + ini_date + '-' + fin_date
    xlabel = 'Observation Current Velocity [m/s]'
    ylabel = 'Model Current Velocity [m/s]'
    if timerange.shape[0] > 2:
        statistics_array_rev_interp=scatterPlot(mod_array_rev_interp,obs_array_rev_interp,output_plot_folder + plotname,name_exp,len(listOfFiles),timerange.shape[0],possible_markers,hfr_names_rev,possible_colors,string_time_res,title=title,xlabel=xlabel,ylabel=ylabel)
        row_all_rev_interp = mean_all_rev_interp + statistics_array_rev_interp
        statistics_rev_interp["ALL HFR STATIONS"] = row_all_rev_interp

        a_file = open(work_dir+"statistics_rev_interp_" + name_exp + "_" + ini_date + "_" + fin_date + ".csv", "w")
        writer = csv.writer(a_file)
        writer.writerow(["name_hfr", "mean_mod", "mean_obs", "bias","rmse","si","corr","stderr","number_of_points"])
        for key, value in statistics_rev_interp.items():
            array = [key] + value
            print(array)
            writer.writerow(array)
        a_file.close()

