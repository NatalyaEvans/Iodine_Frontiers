# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 07:17:55 2022

@author: Natalya Evans

Purpose: plot a map of stations
"""
#%% Initialize packages


#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import netCDF4


#%% set up map boundaries and graphical properties

def format_string(lonlat): #to remove the +/= and degree marks on ticks when drawing parallels and meridians
    if(lonlat>180):
       return "{num}".format(num=lonlat-360)
    return "{num}".format(num=lonlat)

# Nedge=27
# Sedge=7
Nedge=30
Sedge=3.5
Wedge=-123
Eedge=-100

longlineE=Eedge
longlineW=Wedge
longlinespacing=5

latlineN=Nedge
latlineS=5
latlinespacing=5


#%% load in data to plot

df = pd.read_csv ('FK_ompa_04no3.csv')

#%% Process data to plot
df2 = df[['long', 'lat']].copy()
df2['lat']=round(df2['lat'],3) # coverts everything to 3 decimal places
lat=df2['lat']
long=df2['long']

#%% Read in the O2 data

file2read = netCDF4.Dataset('data-set4.nc','r')
print(file2read)
print(file2read.variables.keys())

x = np.array(file2read.variables["Longitude"][:]) # var can be 'defined by he file2read variable keys.
y = np.array(file2read.variables["Latitude"][:]) # var can be 'defined by he file2read variable keys.
Density = np.array(file2read.variables["Density"][:]) # var can be 'defined by he file2read variable keys.
fODZ = np.array(file2read.variables["fODZ"][:]) # var can be 'defined by he file2read variable keys.
fODZ_slice=fODZ[12,:,:]
X, Y = np.meshgrid(x, y) # mesh the x and y

#%% Plot the map with the previous settings

fig = plt.figure(figsize=(8, 8))
plt.rcParams["font.family"] = "sans"
m = Basemap(projection='cyl', resolution='i',
            llcrnrlat=Sedge, urcrnrlat=Nedge,
            llcrnrlon=Wedge, urcrnrlon=Eedge, )
m.shadedrelief(scale=0.2, alpha=0.5)
m.drawcoastlines()
m.drawparallels(np.arange(latlineS,latlineN,latlinespacing),labels=[1,0,0,0], fontsize=11, labelstyle= "+/-", fmt=format_string)
m.drawmeridians(np.arange(longlineW,longlineE,longlinespacing),labels=[1,1,0,1], fontsize=11, labelstyle= "+/-", fmt=format_string)

clev = np.arange(0,1,0.01) #Adjust the .001 to get finer gradient
# C=plt.contourf(X,Y,fODZ_slice,clev,cmap='Blues_r',)
# cbar = plt.colorbar(shrink = 0.825,)
C=plt.contourf(X,Y,fODZ_slice,clev,cmap='BuPu',alpha=0.75, antialiased=True)
cbar = plt.colorbar(shrink = 0.945,)
cbar.ax.set_ylabel('ODZ fraction at 26.5 kg $m^{-3}$')
plt.plot(long,lat,linestyle='None',marker="o",color="k")


plt.xlabel('Longitude/°E', labelpad=20, fontsize=11)
plt.ylabel('Latitude/°N', labelpad=35, fontsize=11)
plt.legend()
plt.savefig('FK180624_map.tif', dpi='figure')


