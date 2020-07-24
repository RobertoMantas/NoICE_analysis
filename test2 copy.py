from __future__ import print_function
import os
from netCDF4 import Dataset
from wrf import getvar
from wrf import getvar, to_np, get_cartopy, latlon_coords
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta
import tkinter
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import griddata
import scipy.interpolate as spint
import scipy.interpolate as interpolate
import scipy.spatial.qhull as qhull
import itertools
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

mpl.use('TkAgg')

def WRF_extract_list(in_dir, out_dir):
    for base, dirs, files in os.walk(in_dir):
        for file in files:
            if str(file).startswith('wrfout_d02_'):
                wrf_file = Dataset(in_dir+file)
                with open(out_dir+ '/OUTPUT_'+ str(file)+".dat" ,'w') as output_file:
                    #If we are reading ncfiles:
                    #lons = wrf_file.variables['TH2'][:]
                    # Changing data from xarrays to numpy arrays and reading each variable from WRF file and 
                    # lat: latitude, o
                    lat = getvar(wrf_file, "lat")
                    data_lat = np.asarray(lat)
                    # lon: longitude, o
                    lon = getvar(wrf_file, "lon")
                    data_lon = np.asarray(lon)
                    #TH2: Potential Temperature at 2-meters, K
                    TH2 = getvar(wrf_file, "TH2")
                    data_TH2 = np.asarray(TH2)
                    #T2: Temperature at 2-meters, K
                    T2 = getvar(wrf_file, "T2")
                    data_T2 = np.asarray(T2)
                    #slp: Sea Level Pressure, hPa/mb
                    slp = getvar(wrf_file, "slp")
                    data_slp = np.asarray(slp)
                    #rh2: Relative Humidity at 2-meters, %
                    rh2 = getvar(wrf_file, "rh2")
                    data_rh2 = np.asarray(rh2)
                    #TSK: Surface Skin Temperature, K
                    TSK = getvar(wrf_file, "TSK")
                    data_TSK = np.asarray(TSK)
                    #PSFC: Surface Pressure, Pa
                    PSFC = getvar(wrf_file, "PSFC")
                    data_PSFC = np.asarray(PSFC)
                    #U-component of Wind on Mass Points, m/s
                    U10 = getvar(wrf_file, "U10")
                    data_U10 = np.asarray(U10)
                    #V-component of Wind on Mass Points, m/s
                    V10 = getvar(wrf_file, "V10")
                    data_V10 = np.asarray(V10)
                    wind_speed=np.sqrt(data_U10**2+data_V10**2)
                    #SNOWNC: Accumulated Snow and Ice from the Microphysics Scheme between model output steps, mm
                    SNOWNC = getvar(wrf_file, "SNOWNC")
                    data_SNOWNC = np.asarray(SNOWNC)
                    #Precipitation, Accumulated Precipitation from the Microphysics Scheme between model output steps, mm
                    RAINNC = getvar(wrf_file, "RAINNC")
                    data_RAINNC = np.asarray(RAINNC)
                    #Downward Short-wave Flux at Ground Surface, W m-2
                    SWDNB = getvar(wrf_file, "SWDNB")
                    data_SWDNB = np.asarray(SWDNB)
                    #Upward Short-wave Flux at Ground Surface, W m-2
                    SWUPB = getvar(wrf_file, "SWUPB")
                    data_SWUPB = np.asarray(SWUPB)
                    #Downward Long-wave Flux at Ground Surface, W m-2
                    LWDNB = getvar(wrf_file, "LWDNB")
                    data_LWDNB = np.asarray(LWDNB)
                    #Upward Heat Flux at the Surface, W m-2
                    HFX = getvar(wrf_file, "HFX")
                    data_HFX = np.asarray(HFX)
                    #Latent Heat Flux at the Surface, W m-2
                    LH = getvar(wrf_file, "LH")
                    data_LH = np.asarray(LH)
                    # Ice Depth, m
                    ICEDEPTH = getvar(wrf_file, "ICEDEPTH")
                    data_ICEDEPTH = np.asarray(ICEDEPTH)
                    # Ice Depth, m
                    MMAKO = getvar(wrf_file, "MMAKO")
                    data_MMAKO = np.asarray(MMAKO)
                    # Ice Depth, m
                    RMAKO = getvar(wrf_file, "RMAKO")
                    data_RMAKO = np.asarray(RMAKO)
                    # Ice Depth, m
                    MMAKOW = getvar(wrf_file, "MMAKOW")
                    data_MMAKOW = np.asarray(MMAKOW)
                    # Ice Depth, m
                    RMAKOW = getvar(wrf_file, "RMAKOW")
                    data_RMAKOW = np.asarray(RMAKOW)
                    # Ice Depth, m
                    MJONES = getvar(wrf_file, "MJONES")
                    data_MJONES = np.asarray(MJONES)
                    # Ice Depth, m
                    RJONES = getvar(wrf_file, "RJONES")
                    data_RJONES = np.asarray(RJONES)
                    
                    output_file.write('Latitude, Longitude, Potential_Temperature_2_meters, Temperature_2_meters, Sea_Level_Pressure, RH, Surface_Skin_Temperature, Surface_Pressure, U_10, V10, Wind_Speed, SNOWNC, RAINNC, SWDNB, SWUPB, LWDNB, HFX, LH, ICEDEPTH, MMAKO, RMAKO, MMAKOW, RMAKOW, MJONES, RJONES'+ '\n')
                    for i in range(len(data_lat)):
                        for j in range(len(data_lat[i])):
                            float_lat = data_lat[i][j]
                            float_lon = data_lon[i][j]
                            float_TH2 = data_TH2[i][j]
                            float_T2 = data_T2[i][j]
                            float_slp = data_slp[i][j]
                            float_rh2 = data_rh2[i][j]
                            float_TSK = data_TSK[i][j]
                            float_PSFC = data_PSFC[i][j]
                            float_U10 = data_U10[i][j]
                            float_V10 = data_V10[i][j]
                            float_wind_speed = wind_speed[i][j]
                            float_SNOWNC = data_SNOWNC[i][j]
                            float_RAINNC = data_RAINNC[i][j]
                            float_SWDNB = data_SWDNB[i][j]
                            float_SWUPB = data_SWUPB[i][j]
                            float_LWDNB = data_LWDNB[i][j]
                            float_HFX = data_HFX[i][j]
                            float_LH = data_LH[i][j]
                            float_ICEDEPTH = data_ICEDEPTH[i][j]
                            float_MMAKO = data_MMAKO[4][i][j]
                            float_RMAKO = data_RMAKO[4][i][j]
                            float_MMAKOW = data_MMAKOW[4][i][j]
                            float_RMAKOW = data_RMAKOW[4][i][j]
                            float_MJONES = data_MJONES[4][i][j]
                            float_RJONES = data_RJONES[4][i][j]
                            output_file.write(str(float_lat) +', '+ str(float_lon) +', '+ str(float_TH2) +', '+ str(float_T2) +', '+ str(float_slp) +', '+ str(float_rh2) +', '+ str(float_TSK) +', '+ str(float_PSFC) +', '+ str(float_U10) +', '+ str(float_V10) +', '+ str(float_wind_speed) +', '+ str(float_SNOWNC) +', '+ str(float_RAINNC) +', '+ str(float_SWDNB) 
                            +', '+ str(float_SWUPB) +', '+ str(float_LWDNB) +', '+ str(float_HFX) +', '+ str(float_LH) +', '+ str(float_ICEDEPTH) +', '+ str(float_MMAKO) +', '+ str(float_RMAKO) +', '+ str(float_MMAKOW) +', '+ str(float_RMAKOW) +', '+ str(float_MJONES) +', '+ str(float_RJONES) + '\n')

                    output_file.write("# lat: latitude, o"+ '\n')
                    output_file.write("# lon: longitude, o"+ '\n')
                    output_file.write("#TH2: Potential Temperature at 2-meters, K"+ '\n')
                    output_file.write("#T2: Temperature at 2-meters, K"+ '\n')
                    output_file.write("#slp: Sea Level Pressure, hPa/mb"+ '\n')
                    output_file.write("#rh2: Relative Humidity at 2-meters, %"+ '\n')
                    output_file.write("#TSK: Surface Skin Temperature, K"+ '\n')
                    output_file.write("#PSFC: Surface Pressure, Pa"+ '\n')
                    output_file.write("#U10: U-component of Wind on Mass Points, m/s"+ '\n')
                    output_file.write("#V10: V-component of Wind on Mass Points, m/s"+ '\n')
                    output_file.write("#Wind Speed sqrt(U10^2+V10^2) of Wind on Mass Points, m/s"+ '\n')
                    output_file.write("#SNOWNC: Accumulated Snow and Ice from the Microphysics Scheme between model output steps, mm"+ '\n')
                    output_file.write("#Precipitation, Accumulated Precipitation from the Microphysics Scheme between model output steps, mm"+ '\n')
                    output_file.write("#Downward Short-wave Flux at Ground Surface, W m-2"+ '\n')
                    output_file.write("#Upward Short-wave Flux at Ground Surface, W m-2"+ '\n')
                    output_file.write("#Downward Long-wave Flux at Ground Surface, W m-2"+ '\n')
                    output_file.write("#Upward Heat Flux at the Surface, W m-2"+ '\n')
                    output_file.write("#Latent Heat Flux at the Surface, W m-2"+ '\n')
                    output_file.write("# Ice Depth, m"+ '\n')
                    output_file.write("# MASS OF ACCRETED ICE ON IDEALISED CYLINDER PER UNIT LENGTH USING MAKKONEN MODEL FOR IN-CLOUD ICING" "kg m-1"+ '\n')
                    output_file.write("# RADIAL ICE THICKNESS ON IDEALISED CYLINDER USING MAKKONEN MODEL FOR IN-CLOUD ICING, m"+ '\n')
                    output_file.write("# MASS OF ACCRETED ICE ON IDEALISED CYLINDER PER UNIT LENGTH USING MAKKONEN MODEL FOR WET SNOW kg m-1"+ '\n')
                    output_file.write("# RADIAL ICE THICKNESS ON IDEALISED CYLINDER USING MAKKONEN MODEL FOR WET SNOW, m"+ '\n')
                    output_file.write("# MASS OF ACCRETED ICE ON IDEALISED CYLINDER PER UNIT LENGTH USING JONES MODEL FOR FREEZING RAIN, kg·m-1"+ '\n')
                    output_file.write("#RADIAL ICE THICKNESS ON IDEALISED CYLINDER USING JONES MODEL FOR FREEZING RAIN m"+ '\n')

WRF_extract_list("/Users/robertomantas/Downloads/", os.path.dirname(__file__)+"/Output_WRF_list")
