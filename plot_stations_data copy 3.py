import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
in_dir = os.path.dirname(__file__) 
print(in_dir)
col_names = ["Station_name","Station_code","Height","Lat","Long","Date","Hour","Air_Temp","P", "RH", "Wind_Direction", "Wind_Speed"]
df3 = pd.read_csv(in_dir+ "/Umea_Weather_Station_NoICE_ALL.csv", delimiter=",", skiprows=1, names=col_names, engine='python')

fig = plt.figure(figsize=(12,9))

m = Basemap(projection='mill',
           llcrnrlat = -90,
           urcrnrlat = 90,
           llcrnrlon = -180,
           urcrnrlon = 180,
           resolution = 'c')

m.drawcoastlines()

m.drawparallels(np.arange(-90,90,10),labels=[True,False,False,False])
m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1])

sites_lat_y = df['latitude'].tolist()
sites_lon_x = df['longitude'].tolist()

colors = ['green', 'darkblue', 'yellow', 'red', 'blue', 'orange']

m.scatter(sites_lon_x,sites_lat_y,latlon=True, s=500, c=colors, marker='o', alpha=1, edgecolor='k', linewidth=1, zorder=2)
m.scatter(-135,60,latlon=True, s=5000, c='blue', marker='^', alpha=1, edgecolor='k', linewidth=1, zorder=1)

plt.title('Basemap tutorial', fontsize=20)

plt.show()