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
                #If we are reading ncfiles:
                #lons = wrf_file.variables['TH2'][:]
                # Changing data from xarrays to numpy arrays and reading each variable from WRF file and 
                # lat: latitude, o
                lat = getvar(wrf_file, "lat")
                data_lat = np.asarray(lat)
                MMAKO = getvar(wrf_file, "RMAKO")
                data_MMAKO = np.asarray(MMAKO)
                print((data_MMAKO[0]).shape)
                print((data_lat).shape)
WRF_extract_list("/Volumes/Samsung_T5/WRF/WRF_Code_updated/", os.path.dirname(__file__)+"/Output_WRF_list")
