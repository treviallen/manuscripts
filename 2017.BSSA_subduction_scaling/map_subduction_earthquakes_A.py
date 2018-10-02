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

fig = plt.figure(figsize=(18,10.5))
plt.tick_params(labelsize=16)

# set map
m = Basemap(projection='robin', lat_0=0, lon_0=180, 
            resolution='l', area_thresh=10000.0)
#ax = fig.add_subplot(211)

#m.drawcoastlines(color='0.5', linewidth=0.5)
m.shadedrelief()
#m.etopo()
meridians = arange(20,360,40)
m.drawmeridians(meridians,labels=[0,0,0,1], fontsize=14, linewidth=0.5, color='0.5')
parallels = arange(-90,90,30)
m.drawparallels(parallels,labels=[1,0,0,0], fontsize=14, linewidth=0.5, color='0.5')


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

swp_urcrnrlat = 10
swp_llcrnrlat = -32
swp_urcrnrlon = 192
swp_llcrnrlon = 89 

sam_urcrnrlat = 3  
sam_llcrnrlat = -25.5 
sam_urcrnrlon = -65 
sam_llcrnrlon = -85 

kur_urcrnrlat = 58
kur_llcrnrlat = 35
kur_urcrnrlon = 190
kur_llcrnrlon = 138

for i, t in enumerate(typ):
    plttxt = False
    x, y = m(lon[i], lat[i]) 
    xt, yt = m(lon[i]+3.5, lat[i]+2.5)
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
    
    # check if inside SAM box
    if lon[i] > sam_llcrnrlon and lon[i] < sam_urcrnrlon \
       and lat[i] > sam_llcrnrlat and lat[i] < sam_urcrnrlat:
       plttxt = False
    
    # check if inside SWP box
    if lon[i] < 0:
        lon[i] += 360
    
    if lon[i] > swp_llcrnrlon and lon[i] < swp_urcrnrlon \
       and lat[i] > swp_llcrnrlat and lat[i] < swp_urcrnrlat:
       plttxt = False
       
    # check if inside KUR box
    if lon[i] < 0:
        lon[i] += 360
    
    if lon[i] > kur_llcrnrlon and lon[i] < kur_urcrnrlon \
       and lat[i] > kur_llcrnrlat and lat[i] < kur_urcrnrlat:
       plttxt = False
    
    if plttxt == True:
           plt.text(xt, yt, str(evidx[i]), size=13, horizontalalignment='left', verticalalignment='center', weight='normal').set_clip_on(True)

plt.legend((hi[0],hs[0],ht[0],ho[0]), ('Interface', 'Intraslab', 'Strike-Slip', 'Outer-Rise'), \
            fontsize=14, loc=7, numpoints=1)

##
# plt SWP box
swplon = [swp_llcrnrlon, swp_llcrnrlon, swp_urcrnrlon, swp_urcrnrlon, swp_llcrnrlon]
swplat = [swp_llcrnrlat, swp_urcrnrlat, swp_urcrnrlat, swp_llcrnrlat, swp_llcrnrlat]

x, y = m(swplon, swplat)
plt.plot(x, y, 'k--', lw=1.5)

# write text
ft2lat = min(swplat)-1.5
ft2lon = min(swplon)
x, y = m(ft2lon, ft2lat)
plt.text(x, y, 'Figure 2A', size=16, ha='left', va='top', weight='normal', style='italic')

# plt SAM box
samlon = [sam_llcrnrlon, sam_llcrnrlon, sam_urcrnrlon, sam_urcrnrlon, sam_llcrnrlon]
samlat = [sam_llcrnrlat, sam_urcrnrlat, sam_urcrnrlat, sam_llcrnrlat, sam_llcrnrlat]

x, y = m(samlon, samlat)
plt.plot(x, y, 'k--', lw=1.5)

# write text
ft2lat = min(samlat)-1.5
ft2lon = min(samlon)
x, y = m(ft2lon, ft2lat)
plt.text(x, y, 'Figure 2D', size=16, ha='left', va='top', weight='normal', style='italic')

# plt KUR box
kurlon = [kur_llcrnrlon, kur_llcrnrlon, kur_urcrnrlon, kur_urcrnrlon, kur_llcrnrlon]
kurlat = [kur_llcrnrlat, kur_urcrnrlat, kur_urcrnrlat, kur_llcrnrlat, kur_llcrnrlat]

x, y = m(kurlon, kurlat)
plt.plot(x, y, 'k--', lw=1.5)

# write text
ft2lat = min(kurlat)-1.5
ft2lon = min(kurlon)
x, y = m(ft2lon, ft2lat)
plt.text(x, y, 'Figure 2C', size=16, ha='left', va='top', weight='normal', style='italic')


plt.savefig('subduction_events_A.pdf', format='pdf', dpi=300, bbox_inches='tight')

plt.savefig('subduction_events_A.png', format='png', dpi=300, bbox_inches='tight')
plt.show()