import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

detectors = pd.read_csv (r'detectors_params.csv')
detectors = detectors[0:5]

names = pd.DataFrame(detectors, columns= ['name'])['name']
latitudes = pd.DataFrame(detectors, columns= ['North lat'])['North lat']
longitudes = pd.DataFrame(detectors, columns= ['East lon'])['East lon']
arms1 = pd.DataFrame(detectors, columns= ['arm1 azimuth'])['arm1 azimuth']
arms2 = pd.DataFrame(detectors, columns= ['arm2 azimuth'])['arm2 azimuth']

latitudes = latitudes*180/np.pi
for k in range(len(detectors)):
    if longitudes[k]>np.pi:
        longitudes[k]=longitudes[k]-2*np.pi
longitudes = longitudes*180/np.pi

fig = plt.figure()

plt.title('Earth-Based Detector Locations')

plt.xlabel('Longitude')
plt.xlim(-180,180)
xticks = [str(n)+'°W' for n in range(180,0,-30)]
xticks += ['0°']
xticks += [str(n)+'°E' for n in range(30,210,30)]
lon_axis = np.arange(-180,210,30)
plt.xticks(lon_axis,xticks)

plt.ylabel('Latitude')
plt.ylim(-90,90)
yticks = [str(n)+'°S' for n in range(90,0,-30)]
yticks += ['0°']
yticks += [str(n)+'°N' for n in range(30,120,30)]
lat_axis = np.arange(-90,120,30)
plt.yticks(lat_axis,yticks)

colors = ['red','green','purple','blue','orange','black']
loc_txt = [[3,3],[3,3],[5,-5],[3,-8],[1,-5],[10,10]]
for k in range(len(detectors)):
    plt.scatter(longitudes[k],latitudes[k],marker='.',color=colors[k])
    plt.text(longitudes[k]+loc_txt[k][0], latitudes[k]+loc_txt[k][1], names[k], fontsize=9,color=colors[k])
    direction1 = [np.sin(arms1[k]),np.cos(arms1[k])]
    plt.quiver(longitudes[k],latitudes[k],direction1[0],direction1[1],color=colors[k],headlength=0,headwidth=1)
    direction2 = [np.sin(arms2[k]),np.cos(arms2[k])]
    plt.quiver(longitudes[k],latitudes[k],direction2[0],direction2[1],color=colors[k],headlength=0,headwidth=1)
    
plt.grid(True)

fig_name = 'detector_map.png';
plt.savefig(fig_name)
fig_name = 'detector_map.pdf';
plt.savefig(fig_name)

plt.show()