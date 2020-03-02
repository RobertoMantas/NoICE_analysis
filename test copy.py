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
                print(file)
                wrf_file = Dataset(in_dir+file)
                with open(out_dir+ '/OUTPUT_'+ str(file)+".dat" ,'w') as output_file:
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

def WRF_vs_OBS(lat_station, long_station, year, month, day):

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
                print(date)
                Latitude = np.array(df["Latitude"])
                Longitude = np.array(df["Longitude"])
                Potential_Temperature_2_meters = np.array(df["Potential_Temperature_2_meters"])
                Temperature_2_meters = np.array(df["Temperature_2_meters"])
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

                grid_tempPotential_Temperature_2_meters = griddata((Latitude, Longitude), Potential_Temperature_2_meters, (lat_station, long_station), method='cubic')
                grid_tempTemperature_2_meters = griddata((Latitude, Longitude), Temperature_2_meters, (lat_station, long_station), method='cubic')
                grid_tempSea_Level_Pressure = griddata((Latitude, Longitude), Sea_Level_Pressure, (lat_station, long_station), method='cubic')
                grid_tempRH = griddata((Latitude, Longitude), RH, (lat_station, long_station), method='cubic')
                grid_tempSurface_Skin_Temperature = griddata((Latitude, Longitude), Surface_Skin_Temperature, (lat_station, long_station), method='cubic')
                grid_tempSurface_Pressure = griddata((Latitude, Longitude), Surface_Pressure, (lat_station, long_station), method='cubic')
                grid_tempU_10 = griddata((Latitude, Longitude), U_10, (lat_station, long_station), method='cubic')
                grid_tempV10 = griddata((Latitude, Longitude), V10, (lat_station, long_station), method='cubic')
                grid_tempWind_Speed = griddata((Latitude, Longitude), Wind_Speed, (lat_station, long_station), method='cubic')
                grid_tempSNOWNC = griddata((Latitude, Longitude), SNOWNC, (lat_station, long_station), method='cubic')
                grid_tempRAINNC = griddata((Latitude, Longitude), RAINNC, (lat_station, long_station), method='cubic')
                grid_tempSWDNB = griddata((Latitude, Longitude), SWDNB, (lat_station, long_station), method='cubic')
                grid_tempSWUPB = griddata((Latitude, Longitude), SWUPB, (lat_station, long_station), method='cubic')
                grid_tempLWDNB = griddata((Latitude, Longitude), LWDNB, (lat_station, long_station), method='cubic')
                grid_tempHFX = griddata((Latitude, Longitude), HFX, (lat_station, long_station), method='cubic')
                grid_tempLH = griddata((Latitude, Longitude), LH, (lat_station, long_station), method='cubic')
                grid_tempICEDEPTH = griddata((Latitude, Longitude), ICEDEPTH, (lat_station, long_station), method='cubic')
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

WRF_vs_OBS(63.7947, 20.2918, "2016", "12", "03")

with open(in_dir + '/Umea_Weather_Station_NoICE2_CSV.csv','r') as input_file, open(in_dir + '/WRF_Umea-data.csv','r') as input_file2, open(in_dir + '/METAR_Umea-data.csv','r') as input_file3:
                    
    fig, ax1 = plt.subplots(figsize=(9, 9))
    
    for (i, j, z) in zip(input_file, input_file2, input_file3):
        col1=i.split(",")
        if i.startswith("Station_name"):
            continue
        col2=j.split(",")
        if j.startswith("Date"):
            continue
        col3=z.split(",")
        if z.startswith("station"):
            continue

        ############ SMHI STATION ###############################

        station_name = (col1[0])
        station_code = float(col1[1])
        station_height = float(col1[2])
        lat = float(col1[3])
        long = float(col1[4])
        date = (col1[5])
        hour = (col1[6])
        #Celsius
        air_Temp =float(col1[7])+kelvin
        #(hPa) 
        P =float(col1[8])
        # (%) 
        RH = float(col1[9])
        # (DegC)
        wind_direction =float(col1[10])
        # (m/s) 
        wind_speed = float(col1[11])
        ############ WRF DATA ###############################

        date2 = (col2[0])
        hour2 = (col2[1])
        lat2 = float(col2[2])
        long2 = float(col2[3])
        Potential_Temperature_2_meters = float(col2[4])
        #K
        air_Temp2 = float(col2[5])
        Sea_Level_Pressure = float(col2[6])
        #%
        RH2 =float(col2[7])
        #(K) 
        SFT2 =float(col2[8])
        # (%) 
        P2 = float(col2[9])/100 #TODO remove +2
        wind_speed2 = float(col2[10])
        wind_speed2_2 = float(col2[11])

        ############ METAR STATION ###############################

        station_name3 = (col3[0])
        date3 = (col3[1])
        lon3 = float(col3[2])
        lat3 = float(col3[3])
        air_temp3 = float(col3[4])+kelvin
        RH3 = float(col3[5])
        WD3 = float(col3[6])
        #Celsius
        WS3 = (float(col3[8])/2.237)
        #(hPa) 
        #SLP3 = float(col3[9])



        timestamp = date + ' '+ hour
        timestamp2 = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        
        #if date == "2016-12-03":
        array_time.append(timestamp2)
        array_1.append(air_Temp)
        array_2.append(air_Temp2)
        array_3.append(P)
        array_4.append(Sea_Level_Pressure)
        array_5.append(RH)
        array_6.append(RH2)
        array_7.append(SFT2)
        array_8.append(wind_speed)
        array_9.append(wind_speed2_2)
        array_10.append(air_temp3)
        array_11.append(WS3)
        array_12.append(RH3)
        
    res_airT = np.subtract(array_1, array_2)
    res_P = np.subtract(array_3, array_4)
    res_RH = np.subtract(array_5, array_6)
    
    # print(res_airT)
    # print(res_P)
    # print(res_RH)

    #CORRELATIONS
    print(np.corrcoef(array_1, array_2))
    print(np.corrcoef(array_3, array_4))
    print(np.corrcoef(array_5, array_6))
    print(np.corrcoef(array_8, array_9))
    print("")
    print(pearsonr(array_1, array_2))
    print(pearsonr(array_3, array_4))
    print(pearsonr(array_5, array_6))
    print(pearsonr(array_8, array_9))
    print("")
    print(spearmanr(array_1, array_2))
    print(spearmanr(array_3, array_4))
    print(spearmanr(array_5, array_6))
    print(spearmanr(array_8, array_9))
    MSE = np.square(np.subtract(array_8,array_9)).mean()
    print(MSE)
    #### BIAS: as bias is the difference on average between the true parameter and the estimate and unless you have simulated the data you will not know this.

    # plt.scatter(array_1, array_2)
    # plt.show()
    col_names = ["Station_name","Station_code","Height","Lat","Long","Date","Hour", "AirTemp", "P" ,"RH" ,"Wind_Direction", "Wind_Speed"]

    df = pd.read_csv(in_dir+ "/Output_WRF_list/"+file, delimiter=",", skiprows=1, names=col_names, engine='python')
    AirTemp_OBS = df["AirTemp"]

    mpl.rcParams['axes.linewidth'] = 2.5 #set the value globally
    font = {
            'weight' : 'bold',
            'size'   : 15}

    mpl.rc('font', **font)
    myFmt = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(myFmt)
    ax1.plot(array_time, AirTemp_OBS, color ='green',marker = '^',  linestyle='-', label='WS OBS')
    #ax1.plot(array_time, array_11, color ='red',marker = '^',  linestyle='-', label='WS METAR')
    #ax1.plot(array_time, array_9, color ='blue',marker = '^',  linestyle='-', label='WS WRF')
    ax1.set_ylabel('Wind Speed (m/s)',weight='bold')
    ax1.set_xlabel('Local Time',weight='bold' )
    ax1.legend(loc = 'upper left')
    # ax2 = ax1.twinx()
    # ax2.plot(array_time, array_4, color ='blue',marker = '^',  linestyle='-', label='P-OBS')
    # ax2.plot(array_time, array_4, color ='orange',marker = '^',  linestyle='-', label='P-WRF')
    # ax2.legend(loc = 'upper right')
    # ax2.set_ylabel('Pressure (hPa)',weight='bold')
    # ax2.plot(array_time, array_5, color ='blue',marker = '^',  linestyle='-', label='RH-OBS')
    # ax2.plot(array_time, array_6, color ='orange',marker = '^',  linestyle='-', label='RH-WRF')
    # ax2.set_ylabel('RH (%)',weight='bold')
    # ax2.legend(loc = 'upper right')
    
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
    plt.grid()
    

    plt.title('Umea Station Comparison',weight='bold')
    plt.show()

