from __future__ import print_function
import os
from netCDF4 import Dataset
from wrf import getvar
from wrf import getvar, to_np, get_cartopy, latlon_coords
import pandas as pd
import numpy as np

in_dir = "/Volumes/GAS_HDD/SABER2/026/"

def WRF_extract_list(in_dir, out_dir):
    for base, dirs, files in os.walk(in_dir):
        for file in files:
            if str(file).startswith('SABER_L2A_2002026_00730_02.07.nc'):
                wrf_file = Dataset(in_dir+file)
                with open(out_dir+ '/OUTPUT_'+ str(file)+".dat" ,'w') as output_file:
                    #If we are reading ncfiles:
                    #lons = wrf_file.variables['TH2'][:]
                    # Changing data from xarrays to numpy arrays and reading each variable from WRF file and 
                    # lat: latitude, o
                    lat = getvar(wrf_file, "pressure")
                    data_lat = np.asarray(lat)
                    print(data_lat)
                    # lon: longitude, o

WRF_extract_list(in_dir, in_dir)
