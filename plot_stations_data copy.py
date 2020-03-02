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
                    
                    output_file.write('Latitude, Longitude, Potential_Temperature_2_meters, Temperature_2_meters, Sea_Level_Pressure, RH, Surface_Skin_Temperature, Surface_Pressure, U_10, V10, Wind_Speed, SNOWNC, RAINNC, SWDNB, SWUPB, LWDNB, HFX, LH, ICEDEPTH'+ '\n')
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
                            output_file.write(str(float_lat) +', '+ str(float_lon) +', '+ str(float_TH2) +', '+ str(float_T2) +', '+ str(float_slp) +', '+ str(float_rh2) +', '+ str(float_TSK) +', '+ str(float_PSFC) +', '+ str(float_U10) +', '+ str(float_V10) +', '+ str(float_wind_speed) +', '+ str(float_SNOWNC) +', '+ str(float_RAINNC) +', '+ str(float_SWDNB) 
                            +', '+ str(float_SWUPB) +', '+ str(float_LWDNB) +', '+ str(float_HFX) +', '+ str(float_LH) +', '+ str(float_ICEDEPTH) + '\n')

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

#WRF_extract_list("/Volumes/Samsung_T5/WRF/WRF_Code_updated/", os.path.dirname(__file__)+"/Output_WRF_list")

def interp_weights(xyz, uvw):
    d = 2
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)

in_dir = os.path.dirname(__file__) 
array_time = []
array_1 = []
array_2 = []
array_3 = []
array_4 = []
array_5 = []
array_6 = []
array_7 = []
array_8 = []
array_9 = []
array_10 = []
array_11 = []
array_12 = []
array_13 = []
array_14 = []
array_15 = []
array_16 = []
kelvin = 273.15

def grid_WRF_OBS(lat_station, long_station, year, month, day):

    array_gridded_Potential_Temperature_2_meters = []
    array_gridded_Temperature_2_meters = []
    array_gridded_Sea_Level_Pressure = []
    array_gridded_RH = []
    array_gridded_Surface_Skin_Temperature = []
    array_gridded_Surface_Pressure = []
    array_gridded_U_10 = []
    array_gridded_V10 = []
    array_gridded_Wind_Speed = []
    array_gridded_SNOWNC = []
    array_gridded_RAINNC = []
    array_gridded_SWDNB = []
    array_gridded_SWUPB = []
    array_gridded_LWDNB = []
    array_gridded_HFX = []
    array_gridded_LH = []
    array_gridded_ICEDEPTH = []
    array_date = []
    col_names=["Latitude", "Longitude", "Potential_Temperature_2_meters", "Temperature_2_meters", "Sea_Level_Pressure", "RH", "Surface_Skin_Temperature", "Surface_Pressure", "U_10", "V10", "Wind_Speed", "SNOWNC", "RAINNC", "SWDNB", "SWUPB", "LWDNB", "HFX", "LH","ICEDEPTH"]
    for base, dirs, files in os.walk(in_dir+ "/Output_WRF_list/"):
        files2 = sorted(files)
        for file in files2:
            if str(file).endswith('.dat') and file[18:22] == year and file[23:25] == month and file[26:28] == day:
                df = pd.read_csv(in_dir+ "/Output_WRF_list/"+file, delimiter=",", skiprows=1, skipfooter=19, names=col_names, engine='python')
                date = file[18:28] + " " + file[29:37]
                Latitude = np.array(df["Latitude"])
                Longitude = np.array(df["Longitude"])
                Potential_Temperature_2_meters = np.array(df["Potential_Temperature_2_meters"])-kelvin
                Temperature_2_meters = np.array(df["Temperature_2_meters"])-kelvin
                Sea_Level_Pressure = np.array(df["Sea_Level_Pressure"])
                RH = np.array(df["RH"])
                Surface_Skin_Temperature = np.array(df["Surface_Skin_Temperature"])
                Surface_Pressure = np.array(df["Surface_Pressure"])
                U_10 = np.array(df["U_10"])
                V10 = np.array(df["V10"])
                Wind_Speed = np.array(df["Wind_Speed"])
                SNOWNC = np.array(df["SNOWNC"])
                RAINNC = np.array(df["RAINNC"])
                SWDNB = np.array(df["SWDNB"])
                SWUPB = np.array(df["SWUPB"])
                LWDNB = np.array(df["LWDNB"])
                HFX = np.array(df["HFX"])
                LH = np.array(df["LH"])
                ICEDEPTH = np.array(df["ICEDEPTH"])

                #Calculating the gridded data 
                d = 2
                all_WRF_coordinates = np.column_stack((Latitude, Longitude))
                station_coordinates = np.column_stack((lat_station, long_station))
                vtx, wts = interp_weights(all_WRF_coordinates, station_coordinates)


                grid_tempPotential_Temperature_2_meters = float(interpolate(Potential_Temperature_2_meters, vtx, wts))
                grid_tempTemperature_2_meters = float(interpolate(Temperature_2_meters, vtx, wts))
                grid_tempSea_Level_Pressure = float(interpolate(Sea_Level_Pressure, vtx, wts))
                grid_tempRH = float(interpolate(RH, vtx, wts))
                grid_tempSurface_Skin_Temperature = float(interpolate(Surface_Skin_Temperature, vtx, wts))
                grid_tempSurface_Pressure = float(interpolate(Surface_Pressure, vtx, wts))
                grid_tempU_10 = float(interpolate(U_10, vtx, wts))
                grid_tempV10 = float(interpolate(V10, vtx, wts))
                grid_tempWind_Speed = float(interpolate(Wind_Speed, vtx, wts))
                grid_tempSNOWNC = float(interpolate(SNOWNC, vtx, wts))
                grid_tempRAINNC = float(interpolate(RAINNC, vtx, wts))
                grid_tempSWDNB = float(interpolate(SWDNB, vtx, wts))
                grid_tempSWUPB = float(interpolate(SWUPB, vtx, wts))
                grid_tempLWDNB = float(interpolate(LWDNB, vtx, wts))
                grid_tempHFX = float(interpolate(HFX, vtx, wts))
                grid_tempLH = float(interpolate(LH, vtx, wts))
                grid_tempICEDEPTH = float(interpolate(ICEDEPTH, vtx, wts))

                array_date.append(date)
                array_gridded_Potential_Temperature_2_meters.append(float(grid_tempPotential_Temperature_2_meters))
                array_gridded_Temperature_2_meters.append(grid_tempTemperature_2_meters)
                array_gridded_Sea_Level_Pressure.append(grid_tempSea_Level_Pressure)
                array_gridded_RH.append(grid_tempRH)
                array_gridded_Surface_Skin_Temperature.append(grid_tempSurface_Skin_Temperature)
                array_gridded_Surface_Pressure.append(grid_tempSurface_Pressure)
                array_gridded_U_10.append(grid_tempU_10)
                array_gridded_V10.append(grid_tempV10)
                array_gridded_Wind_Speed.append(grid_tempWind_Speed)
                array_gridded_SNOWNC.append(grid_tempSNOWNC)
                array_gridded_RAINNC.append(grid_tempRAINNC)
                array_gridded_SWDNB.append(grid_tempSWDNB)
                array_gridded_SWUPB.append(grid_tempSWUPB)
                array_gridded_LWDNB.append(grid_tempLWDNB)
                array_gridded_HFX.append(grid_tempHFX)
                array_gridded_LH.append(grid_tempLH)
                array_gridded_ICEDEPTH.append(grid_tempICEDEPTH)
                
    return array_gridded_Potential_Temperature_2_meters, array_gridded_Temperature_2_meters, array_gridded_Sea_Level_Pressure, array_gridded_RH, array_gridded_Surface_Skin_Temperature, array_gridded_Surface_Pressure, array_gridded_U_10, array_gridded_V10, array_gridded_Wind_Speed, array_gridded_SNOWNC, array_gridded_RAINNC, array_gridded_SWDNB, array_gridded_SWUPB, array_gridded_LWDNB, array_gridded_HFX, array_gridded_LH, array_gridded_ICEDEPTH



                
fig, ax1 = plt.subplots(figsize=(9, 9))

#Reading the weather station data.
col_names = ["Station_name","Station_code","Height","Lat","Long","Date","Hour","Air_Temp","P", "RH", "Wind_Direction", "Wind_Speed"]
df3 = pd.read_csv(in_dir+ "/Vastmarkum_Weather_Station_NoICE_ALL.csv", delimiter=",", skiprows=1, names=col_names, engine='python')

############ SMHI STATION ###############################

Station_name_3 = np.array(df3["Station_name"])
Lat = np.array(df3["Lat"])
Long = np.array(df3["Long"])

Date_3 = np.array(df3["Date"][0:24])
Date_4 = np.array(df3["Date"][24:48])
Date_5 = np.array(df3["Date"][48:72])
Date_6 = np.array(df3["Date"][72:96])
Date_7 = np.array(df3["Date"][96:120])
Date_8 = np.array(df3["Date"][120:144])
Date_9 = np.array(df3["Date"][144:168])

Hour_3 = np.array(df3["Hour"][0:24])
Hour_4 = np.array(df3["Hour"][24:48])
Hour_5 = np.array(df3["Hour"][48:72])
Hour_6 = np.array(df3["Hour"][72:96])
Hour_7 = np.array(df3["Hour"][96:120])
Hour_8 = np.array(df3["Hour"][120:144])
Hour_9 = np.array(df3["Hour"][144:168])

Air_Temp_3 = np.array(df3["Air_Temp"][0:24])   #Celsius
Air_Temp_4 = np.array(df3["Air_Temp"][24:48])
Air_Temp_5 = np.array(df3["Air_Temp"][48:72])
Air_Temp_6 = np.array(df3["Air_Temp"][72:96])
Air_Temp_7 = np.array(df3["Air_Temp"][96:120])
Air_Temp_8 = np.array(df3["Air_Temp"][120:144])
Air_Temp_9 = np.array(df3["Air_Temp"][144:168])

P_3 = np.array(df3["P"][0:24]) #(hPa) 
P_4 = np.array(df3["P"][24:48])
P_5 = np.array(df3["P"][48:72])
P_6 = np.array(df3["P"][72:96])
P_7 = np.array(df3["P"][96:120])
P_8 = np.array(df3["P"][120:144])
P_9 = np.array(df3["P"][144:168])

RH_3 = np.array(df3["RH"][0:24])  # (%) 
RH_4 = np.array(df3["RH"][24:48])
RH_5 = np.array(df3["RH"][48:72])
RH_6 = np.array(df3["RH"][72:96])
RH_7 = np.array(df3["RH"][96:120])
RH_8 = np.array(df3["RH"][120:144])
RH_9 = np.array(df3["RH"][144:168])

Wind_Direction_3 = np.array(df3["Wind_Direction"][0:24])  # (DegC)
Wind_Direction_4 = np.array(df3["Wind_Direction"][24:48])
Wind_Direction_5 = np.array(df3["Wind_Direction"][48:72])
Wind_Direction_6 = np.array(df3["Wind_Direction"][72:96])
Wind_Direction_7 = np.array(df3["Wind_Direction"][96:120])
Wind_Direction_8 = np.array(df3["Wind_Direction"][120:144])
Wind_Direction_9 = np.array(df3["Wind_Direction"][144:168])

Wind_Speed_3 = np.array(df3["Wind_Speed"][0:24]) # (m/s) 
Wind_Speed_4 = np.array(df3["Wind_Speed"][24:48])
Wind_Speed_5 = np.array(df3["Wind_Speed"][48:72])
Wind_Speed_6 = np.array(df3["Wind_Speed"][72:96])
Wind_Speed_7 = np.array(df3["Wind_Speed"][96:120])
Wind_Speed_8 = np.array(df3["Wind_Speed"][120:144])
Wind_Speed_9 = np.array(df3["Wind_Speed"][144:168])

############ WRF DATA ###############################
wrf_list_temp=[]
wrf_list_RH=[]
wrf_list_P=[]
wrf_list_WS=[]

wrf_03_12_16_array_gridded_Potential_Temperature_2_meters, wrf_03_12_16_array_gridded_Temperature_2_meters, wrf_03_12_16_array_gridded_Sea_Level_Pressure, wrf_03_12_16_array_gridded_RH, wrf_03_12_16_array_gridded_Surface_Skin_Temperature, wrf_03_12_16_array_gridded_Surface_Pressure, wrf_03_12_16_array_gridded_U_10, wrf_03_12_16_array_gridded_V10, wrf_03_12_16_array_gridded_Wind_Speed, wrf_03_12_16_array_gridded_SNOWNC, wrf_03_12_16_array_gridded_RAINNC, wrf_03_12_16_array_gridded_SWDNB, wrf_03_12_16_array_gridded_SWUPB, wrf_03_12_16_array_gridded_LWDNB, wrf_03_12_16_array_gridded_HFX, wrf_03_12_16_array_gridded_LH, wrf_03_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "03")
wrf_04_12_16_array_gridded_Potential_Temperature_2_meters, wrf_04_12_16_array_gridded_Temperature_2_meters, wrf_04_12_16_array_gridded_Sea_Level_Pressure, wrf_04_12_16_array_gridded_RH, wrf_04_12_16_array_gridded_Surface_Skin_Temperature, wrf_04_12_16_array_gridded_Surface_Pressure, wrf_04_12_16_array_gridded_U_10, wrf_04_12_16_array_gridded_V10, wrf_04_12_16_array_gridded_Wind_Speed, wrf_04_12_16_array_gridded_SNOWNC, wrf_04_12_16_array_gridded_RAINNC, wrf_04_12_16_array_gridded_SWDNB, wrf_04_12_16_array_gridded_SWUPB, wrf_04_12_16_array_gridded_LWDNB, wrf_04_12_16_array_gridded_HFX, wrf_04_12_16_array_gridded_LH, wrf_04_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "04")
wrf_05_12_16_array_gridded_Potential_Temperature_2_meters, wrf_05_12_16_array_gridded_Temperature_2_meters, wrf_05_12_16_array_gridded_Sea_Level_Pressure, wrf_05_12_16_array_gridded_RH, wrf_05_12_16_array_gridded_Surface_Skin_Temperature, wrf_05_12_16_array_gridded_Surface_Pressure, wrf_05_12_16_array_gridded_U_10, wrf_05_12_16_array_gridded_V10, wrf_05_12_16_array_gridded_Wind_Speed, wrf_05_12_16_array_gridded_SNOWNC, wrf_05_12_16_array_gridded_RAINNC, wrf_05_12_16_array_gridded_SWDNB, wrf_05_12_16_array_gridded_SWUPB, wrf_05_12_16_array_gridded_LWDNB, wrf_05_12_16_array_gridded_HFX, wrf_05_12_16_array_gridded_LH, wrf_05_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "05")
wrf_06_12_16_array_gridded_Potential_Temperature_2_meters, wrf_06_12_16_array_gridded_Temperature_2_meters, wrf_06_12_16_array_gridded_Sea_Level_Pressure, wrf_06_12_16_array_gridded_RH, wrf_06_12_16_array_gridded_Surface_Skin_Temperature, wrf_06_12_16_array_gridded_Surface_Pressure, wrf_06_12_16_array_gridded_U_10, wrf_06_12_16_array_gridded_V10, wrf_06_12_16_array_gridded_Wind_Speed, wrf_06_12_16_array_gridded_SNOWNC, wrf_06_12_16_array_gridded_RAINNC, wrf_06_12_16_array_gridded_SWDNB, wrf_06_12_16_array_gridded_SWUPB, wrf_06_12_16_array_gridded_LWDNB, wrf_06_12_16_array_gridded_HFX, wrf_06_12_16_array_gridded_LH, wrf_06_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "06")
wrf_07_12_16_array_gridded_Potential_Temperature_2_meters, wrf_07_12_16_array_gridded_Temperature_2_meters, wrf_07_12_16_array_gridded_Sea_Level_Pressure, wrf_07_12_16_array_gridded_RH, wrf_07_12_16_array_gridded_Surface_Skin_Temperature, wrf_07_12_16_array_gridded_Surface_Pressure, wrf_07_12_16_array_gridded_U_10, wrf_07_12_16_array_gridded_V10, wrf_07_12_16_array_gridded_Wind_Speed, wrf_07_12_16_array_gridded_SNOWNC, wrf_07_12_16_array_gridded_RAINNC, wrf_07_12_16_array_gridded_SWDNB, wrf_07_12_16_array_gridded_SWUPB, wrf_07_12_16_array_gridded_LWDNB, wrf_07_12_16_array_gridded_HFX, wrf_07_12_16_array_gridded_LH, wrf_07_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "07")
wrf_08_12_16_array_gridded_Potential_Temperature_2_meters, wrf_08_12_16_array_gridded_Temperature_2_meters, wrf_08_12_16_array_gridded_Sea_Level_Pressure, wrf_08_12_16_array_gridded_RH, wrf_08_12_16_array_gridded_Surface_Skin_Temperature, wrf_08_12_16_array_gridded_Surface_Pressure, wrf_08_12_16_array_gridded_U_10, wrf_08_12_16_array_gridded_V10, wrf_08_12_16_array_gridded_Wind_Speed, wrf_08_12_16_array_gridded_SNOWNC, wrf_08_12_16_array_gridded_RAINNC, wrf_08_12_16_array_gridded_SWDNB, wrf_08_12_16_array_gridded_SWUPB, wrf_08_12_16_array_gridded_LWDNB, wrf_08_12_16_array_gridded_HFX, wrf_08_12_16_array_gridded_LH, wrf_08_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "08")
wrf_09_12_16_array_gridded_Potential_Temperature_2_meters, wrf_09_12_16_array_gridded_Temperature_2_meters, wrf_09_12_16_array_gridded_Sea_Level_Pressure, wrf_09_12_16_array_gridded_RH, wrf_09_12_16_array_gridded_Surface_Skin_Temperature, wrf_09_12_16_array_gridded_Surface_Pressure, wrf_09_12_16_array_gridded_U_10, wrf_09_12_16_array_gridded_V10, wrf_09_12_16_array_gridded_Wind_Speed, wrf_09_12_16_array_gridded_SNOWNC, wrf_09_12_16_array_gridded_RAINNC, wrf_09_12_16_array_gridded_SWDNB, wrf_09_12_16_array_gridded_SWUPB, wrf_09_12_16_array_gridded_LWDNB, wrf_09_12_16_array_gridded_HFX, wrf_09_12_16_array_gridded_LH, wrf_09_12_16_array_gridded_ICEDEPTH = grid_WRF_OBS(Lat[1], Long[1], "2016", "12", "09")



wrf_list_temp.extend(wrf_03_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_04_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_05_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_06_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_07_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_08_12_16_array_gridded_Temperature_2_meters)
wrf_list_temp.extend(wrf_09_12_16_array_gridded_Temperature_2_meters)

wrf_list_RH.extend(wrf_03_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_04_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_05_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_06_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_07_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_08_12_16_array_gridded_RH)
wrf_list_RH.extend(wrf_09_12_16_array_gridded_RH)

wrf_list_P.extend(wrf_03_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_04_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_05_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_06_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_07_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_08_12_16_array_gridded_Sea_Level_Pressure)
wrf_list_P.extend(wrf_09_12_16_array_gridded_Sea_Level_Pressure)

wrf_list_WS.extend(wrf_03_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_04_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_05_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_06_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_07_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_08_12_16_array_gridded_Wind_Speed)
wrf_list_WS.extend(wrf_09_12_16_array_gridded_Wind_Speed)

axis_x = (np.linspace(0,len(wrf_list_temp)-1,len(wrf_list_temp)))

#CORRELATIONS
# print(np.corrcoef(Air_Temp, wrf_03_12_16_array_gridded_Temperature_2_meters))
# print("")
# print(pearsonr(Air_Temp, wrf_03_12_16_array_gridded_Temperature_2_meters))
# print("")
print("")

print(wrf_04_12_16_array_gridded_Temperature_2_meters)
print("")

print(spearmanr(Air_Temp_3, wrf_03_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_4, wrf_04_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_5, wrf_05_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_6, wrf_06_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_7, wrf_07_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_8, wrf_08_12_16_array_gridded_Temperature_2_meters))
print(spearmanr(Air_Temp_9, wrf_09_12_16_array_gridded_Temperature_2_meters))
print("")
print(spearmanr(P_3, wrf_03_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_4, wrf_04_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_5, wrf_05_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_6, wrf_06_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_7, wrf_07_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_8, wrf_08_12_16_array_gridded_Sea_Level_Pressure))
print(spearmanr(P_9, wrf_09_12_16_array_gridded_Sea_Level_Pressure))
print("")
print(spearmanr(RH_3, wrf_03_12_16_array_gridded_RH))
print(spearmanr(RH_4, wrf_04_12_16_array_gridded_RH))
print(spearmanr(RH_5, wrf_05_12_16_array_gridded_RH))
print(spearmanr(RH_6, wrf_06_12_16_array_gridded_RH))
print(spearmanr(RH_7, wrf_07_12_16_array_gridded_RH))
print(spearmanr(RH_8, wrf_08_12_16_array_gridded_RH))
print(spearmanr(RH_9, wrf_09_12_16_array_gridded_RH))
print("")
print(spearmanr(Wind_Speed_3, wrf_03_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_4, wrf_04_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_5, wrf_05_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_6, wrf_06_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_7, wrf_07_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_8, wrf_08_12_16_array_gridded_Wind_Speed))
print(spearmanr(Wind_Speed_9, wrf_09_12_16_array_gridded_Wind_Speed))
print("")

# #### BIAS: as bias is the difference on average between the true parameter and the estimate and unless you have simulated the data you will not know this.

print(np.mean(Air_Temp_3)-np.mean(wrf_03_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_4)-np.mean(wrf_04_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_5)-np.mean(wrf_05_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_6)-np.mean(wrf_06_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_7)-np.mean(wrf_07_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_8)-np.mean(wrf_08_12_16_array_gridded_Temperature_2_meters))
print(np.mean(Air_Temp_9)-np.mean(wrf_09_12_16_array_gridded_Temperature_2_meters))
# print((np.mean(Air_Temp_3)+np.mean(Air_Temp_4)+np.mean(Air_Temp_5)+np.mean(Air_Temp_6))-(np.mean(wrf_03_12_16_array_gridded_Temperature_2_meters)+np.mean(wrf_04_12_16_array_gridded_Temperature_2_meters)+np.mean(wrf_05_12_16_array_gridded_Temperature_2_meters)+np.mean(wrf_06_12_16_array_gridded_Temperature_2_meters)))
print("")
print(np.mean(P_3)-np.mean(wrf_03_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_4)-np.mean(wrf_04_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_5)-np.mean(wrf_05_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_6)-np.mean(wrf_06_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_7)-np.mean(wrf_07_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_8)-np.mean(wrf_08_12_16_array_gridded_Sea_Level_Pressure))
print(np.mean(P_9)-np.mean(wrf_09_12_16_array_gridded_Sea_Level_Pressure))
# print((np.mean(P_3)+np.mean(P_4)+np.mean(P_5)+np.mean(P_6))-(np.mean(wrf_03_12_16_array_gridded_Sea_Level_Pressure)+np.mean(wrf_04_12_16_array_gridded_Sea_Level_Pressure)+np.mean(wrf_05_12_16_array_gridded_Sea_Level_Pressure)+np.mean(wrf_06_12_16_array_gridded_Sea_Level_Pressure)))

print("")

print(np.mean(RH_3)-np.mean(wrf_03_12_16_array_gridded_RH))
print(np.mean(RH_4)-np.mean(wrf_04_12_16_array_gridded_RH))
print(np.mean(RH_5)-np.mean(wrf_05_12_16_array_gridded_RH))
print(np.mean(RH_6)-np.mean(wrf_06_12_16_array_gridded_RH))
print(np.mean(RH_7)-np.mean(wrf_07_12_16_array_gridded_RH))
print(np.mean(RH_8)-np.mean(wrf_08_12_16_array_gridded_RH))
print(np.mean(RH_9)-np.mean(wrf_09_12_16_array_gridded_RH))
# print((np.mean(RH_3)+np.mean(RH_4)+np.mean(RH_5)+np.mean(RH_6))-(np.mean(wrf_03_12_16_array_gridded_RH)+np.mean(wrf_04_12_16_array_gridded_RH)+np.mean(wrf_05_12_16_array_gridded_RH)+np.mean(wrf_06_12_16_array_gridded_RH)))

print("")

print(np.mean(Wind_Speed_3)-np.mean(wrf_03_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_4)-np.mean(wrf_04_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_5)-np.mean(wrf_05_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_6)-np.mean(wrf_06_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_7)-np.mean(wrf_07_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_8)-np.mean(wrf_08_12_16_array_gridded_Wind_Speed))
print(np.mean(Wind_Speed_9)-np.mean(wrf_09_12_16_array_gridded_Wind_Speed))
# print((np.mean(Wind_Speed_3)+np.mean(Wind_Speed_4)+np.mean(Wind_Speed_5)+np.mean(Wind_Speed_6))-(np.mean(wrf_03_12_16_array_gridded_Wind_Speed)+np.mean(wrf_04_12_16_array_gridded_Wind_Speed)+np.mean(wrf_05_12_16_array_gridded_Wind_Speed)+np.mean(wrf_06_12_16_array_gridded_Wind_Speed)))

print("")

mpl.rcParams['axes.linewidth'] = 2.5 #set the value globally
font = {
        'weight' : 'bold',
        'size'   : 15}
mpl.rc('font', **font)
# myFmt = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
# fig.autofmt_xdate()
# ax1.xaxis.set_major_formatter(myFmt)
# ax1.plot(axis_x, Air_Temp_3, color ='green',marker = '^',  linestyle='-', label='Air_temp OBS')
# ax1.plot(axis_x, wrf_03_12_16_array_gridded_Temperature_2_meters, color ='red',marker = '^',  linestyle='-', label='Air_temp WRF')
# ax1.set_ylabel('Temp (K)',weight='bold')
# ax1.set_xlabel('Local Time (h)',weight='bold' )
# ax1.legend(loc = 'upper left')


# ax2 = ax1.twinx()
# ax2.plot(array_time, array_4, color ='blue',marker = '^',  linestyle='-', label='P-OBS')
# ax2.plot(array_time, array_4, color ='orange',marker = '^',  linestyle='-', label='P-WRF')
# ax2.legend(loc = 'upper right')
# ax2.set_ylabel('Pressure (hPa)',weight='bold')
# ax2.plot(array_time, array_5, color ='blue',marker = '^',  linestyle='-', label='RH-OBS')
# ax2.plot(array_time, array_6, color ='orange',marker = '^',  linestyle='-', label='RH-WRF')
# ax2.set_ylabel('RH (%)',weight='bold')
# ax2.legend(loc = 'upper right')
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, sharex=True)

MEDIUM_SIZE = 15
# # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# # plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
# #plt.ylim(top=297.5) # adjust the top leaving bottom unchanged
#plt.ylim(bottom=296.5)  # adjust the bottom leaving top unchanged
#plt.grid()


#plt.title('Umea Station Comparison',weight='bold')
fig.suptitle('Vastmarkum Station Comparison',weight='bold')

# plt.show()

ax1.plot(axis_x, df3["Air_Temp"], color ='green',marker = '^',  linestyle='-', label='Air_temp OBS')
ax1.plot(axis_x, wrf_list_temp, color ='red',marker = '^',  linestyle='-', label='Air_temp WRF')
ax1.set_ylabel('Temp (K)',weight='bold')
# plt.xlabel('Local Time (h)',weight='bold' )
ax1.legend(loc = 'upper left')

ax2.plot(axis_x, df3["P"], color ='green',marker = '^',  linestyle='-', label='P OBS')
ax2.plot(axis_x, wrf_list_P, color ='red',marker = '^',  linestyle='-', label='P WRF')
ax2.set_ylabel('P (hPa)',weight='bold')
# plt.xlabel('Local Time (h)',weight='bold' )
ax2.legend(loc = 'upper left')

ax3.plot(axis_x, df3["RH"], color ='green',marker = '^',  linestyle='-', label='RH OBS')
ax3.plot(axis_x, wrf_list_RH, color ='red',marker = '^',  linestyle='-', label='RH WRF')
ax3.set_ylabel('RH (%)',weight='bold')
# plt.xlabel('Local Time (h)',weight='bold' )
ax3.legend(loc = 'upper left')

ax4.plot(axis_x, df3["Wind_Speed"], color ='green',marker = '^',  linestyle='-', label='Wind_Speed OBS')
ax4.plot(axis_x, wrf_list_WS, color ='red',marker = '^',  linestyle='-', label='Wind_Speed WRF')
ax4.set_ylabel('Wind Speed (m/s^2)',weight='bold')
ax4.set_xlabel('Local Time (h)',weight='bold' )
ax4.legend(loc = 'upper left')
plt.show()

    