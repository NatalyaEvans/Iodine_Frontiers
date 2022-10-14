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
import scipy.io as sio

#%% set up map boundaries and graphical properties

def format_string(lonlat): #to remove the +/= and degree marks on ticks when drawing parallels and meridians
    if(lonlat>180):
       return "{num}".format(num=lonlat-360)
    return "{num}".format(num=lonlat)

Nedge=45.5
Sedge=43.5
# Wedge=-127.01
Wedge=-125.5
Eedge=-123.75

longlineE=Eedge
longlineW=Wedge
longlinespacing=1

latlineN=Nedge
latlineS=44
latlinespacing=0.5

markers=["o","^","s","d","P","p","x","o"]
#colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] # extracts the default color cycle
#colors[3], colors[4]=colors[4], colors[3]

#%% load in data to plot

df = pd.read_csv ('OC Master Sample Sheet - Data.csv')
bathy = sio.loadmat('OC2107A_ETOP01.mat')
bathyx=bathy['coast_x']
bathyy=bathy['coast_y']
bathyz=bathy['coastal_relief_3sec']
bathyx2=np.tile(bathyx,[len(bathyy),1])
bathyy2=np.tile(bathyy,[1,7200])

cval=1500
bathyz[bathyz<-cval]=-cval # force values to not mess up the colorbar


#%% Process data to plot
# df2 = df[['Cruise','Station','Cast type', 'Longitude', 'Latitude','Depth/m','Bottom depth/m','Fe(II)/nM','O2/umolkg-1']].copy()
df2=df[df['Station'].isin(['3','4','7','13','31','32','33','34.5','35','36'])]
df3=df[df['Station'].isin(['MT0','MT2','NT5'])]
df3=df3[df3['Type of  cast']=='BBG']

#%% Plot the map with the OC2107A data for Fe(II)

plotdf=df2[df2['Cruise']=='OC2107A']
plotdf.reset_index(drop=True, inplace=True)

fig = plt.figure(figsize=(8, 8))
plt.rcParams["font.family"] = "sans"
m = Basemap(projection='cyl', resolution='i',
            llcrnrlat=Sedge, urcrnrlat=Nedge,
            llcrnrlon=Wedge, urcrnrlon=Eedge, )
#m.shadedrelief(scale=0.2, alpha=0.5)
# plt.contourf(bathyx2,bathyy2,bathyz,levels=np.linspace(-4000,50,50),cmap='ocean',vmin=-6000,vmax=0)
plt.contourf(bathyx2,bathyy2,-bathyz,levels=np.linspace(0,cval,50),cmap='Blues',vmin=0,vmax=cval)
m.drawcoastlines()
m.drawparallels(np.arange(latlineS,latlineN,latlinespacing),labels=[1,0,0,0], fontsize=11, labelstyle= "+/-", fmt=format_string)
m.drawmeridians(np.arange(longlineW,longlineE,longlinespacing),labels=[1,1,0,1], fontsize=11, labelstyle= "+/-", fmt=format_string)

plt.clim(0,cval)
cbar = plt.colorbar(shrink = 0.95,ticks=np.linspace(0,cval,5))
cbar.ax.set_ylabel('Seafloor depth/m',fontsize=11,)

s=plt.scatter(plotdf['Longitude/(ºW)'],plotdf['Latitude/(ºN)'],s=50,edgecolors='m',facecolor='m')
for i in range(len(plotdf['Station'])):  
        plt.annotate(plotdf['Station'][i], (plotdf['Longitude/(ºW)'][i], plotdf['Latitude/(ºN)'][i]),textcoords="offset points",xytext=(0,7),ha='center')

plotdf=df3
plotdf.reset_index(drop=True, inplace=True)
s2=plt.scatter(plotdf['Longitude/(ºW)'],plotdf['Latitude/(ºN)'],s=50,edgecolors='k',facecolor='k',marker='d')
for i in range(len(plotdf['Station'])):  
        plt.annotate(plotdf['Station'][i], (plotdf['Longitude/(ºW)'][i], plotdf['Latitude/(ºN)'][i]),textcoords="offset points",xytext=(0,7),ha='center')

plt.legend([s,s2], ['OC2017A','OC2111A'],loc='upper left')
plt.ylabel('Latitude/°N', labelpad=35, fontsize=11)
plt.xlabel('Longitude/°E', labelpad=35, fontsize=11)
plt.savefig('Iodine_OC2107A_OC2111A_map.tif', dpi='figure')

