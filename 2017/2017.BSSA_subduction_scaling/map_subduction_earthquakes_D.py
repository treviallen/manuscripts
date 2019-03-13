import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import array, arange, log10, where, median, std, mean, concatenate, vstack, argsort
from mpl_toolkits.basemap import Basemap

plt.close('all')
# for bilinear sloping segment

# use historical interface data
usehistoric = True
plt_inter_only = True

import pickle
from shakemap_tools import *
from mapping_tools import distance, reckon
from numpy import array, sqrt, nan, isnan, arange, abs, unique, hstack, savetxt, \
                  logical_and, mean, median, std, log10, ones, logspace, exp, max, signbit, \
                  percentile, random, cov, dot, argsort
#from make_slab_fault import make_slab_matrix
from os import path, sep
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(10,10))
plt.tick_params(labelsize=16)

# set map
urcrnrlat = 55
llcrnrlat = 35
urcrnrlon = 190
llcrnrlon = 140

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
'''
m = Basemap(projection='merc', \
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
'''            
#ax = fig.add_subplot(211)

#m.drawcoastlines(color='0.5', linewidth=0.5)
m.shadedrelief()
#m.etopo()
meridians = arange(20,360,10)
m.drawmeridians(meridians,labels=[0,0,0,1], fontsize=14, linewidth=0.5, color='0.5')
parallels = arange(-90,90,5)
m.drawparallels(parallels,labels=[1,0,0,0], fontsize=14, linewidth=0.5, color='0.5')
m.drawcountries(linewidth=0.75, color='0.2')


# now load trimmed faults
ftxt = open('//Users//tallen//Documents//PAGER//Data_Prep//trimmed_faults_dip_type.csv').readlines()[1:]
date = []
lon = []
lat = []
typ = []
mag = []

for line in ftxt:
    dat = line.strip().split(',')
    if dat[-1] == 'i' or dat[-1] == 's' or dat[-1] == 't' \
       or dat[-1] == 'o' or dat[-1] == 'h':
        date.append(dat[0])
        mag.append(float(dat[3]))
        lon.append(float(dat[4]))
        lat.append(float(dat[5]))
        typ.append(dat[-1])

evtypes = unique(array(typ))
evidx = argsort(argsort(date)) 
evidx += 1 # inc up one

# get types
ilat = []
ilon = []
imag = []
slat = []
slon = []
smag = []
tlat = []
tlon = []
tmag = []
olat = []
olon = []
omag = []

for i, t in enumerate(typ):
    plttxt = False
    x, y = m(lon[i], lat[i]) 
    xt, yt = m(lon[i]-0.5, lat[i]+0.5)
    if t == 'i' or t == 'h':
        hi = plt.plot(x, y, '^', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1.5, markersize=12)
        plttxt = True
    elif t == 's':
        hs = plt.plot(x, y, 'H', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1.5, markersize=13)
        plttxt = True
    elif t == 't':
        ht = plt.plot(x, y, 'd', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1.5, markersize=12)
        plttxt = True
    elif t == 'o':
        ho = plt.plot(x, y, 's', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1.5, markersize=12)
        plttxt = True
        
    if plttxt == True:
           plt.text(xt, yt, str(evidx[i]), size=13, ha='right', va='center', weight='normal').set_clip_on(True)
    
plt.legend((hi[0],hs[0],ho[0]), ('Interface', 'Intraslab', 'Outer-Rise'), \
            fontsize=14, loc='lower right', numpoints=1)

plt.savefig('subduction_events_D.pdf', format='pdf', dpi=300, bbox_inches='tight')    
plt.savefig('subduction_events_D.png', format='png', dpi=300, bbox_inches='tight')
plt.show()