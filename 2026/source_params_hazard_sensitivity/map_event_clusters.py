#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort, log10, unique, nanmean, nanstd
#from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap, remove_last_cmap_colour, remove_first_cmap_colour
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, get_mpl2_colourlist
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn
import shapefile

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -1.0
urcrnrlat = -6.0
llcrnrlat = -45.
urcrnrlon = 153.
llcrnrlon = 106
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
#ax = fig.add_subplot(111)
ax = plt.subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=500.) #i

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawmapboundary(fill_color='0.8')
m.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# set colours
##########################################################################################

# get stress range
logstressbins = arange(-.7, 1.7, 0.2)
logstressrng = logstressbins[-1] - logstressbins[0]

ncols = len(logstressbins) - 1

#logstressrng = float(round(maxstress - minstress))

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
#cptfile = '//Users//trev//Documents//DATA//GMT//cpt//plasma.cpt'

cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cols = (cmap(arange(ncols)))
'''
cmap = plt.get_cmap('viridis_r', len(logstressbins))
cs = (cmap(arange(len(logstressbins))))
'''
##########################################################################################
# add earthquakes
##########################################################################################

#csvfile = '../../2023/au_stress_drop/brune_stats.csv'
csvfile = 'brune_stats_cluster.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

print(csvfile)
lines = open(csvfile).readlines()[1:]
lats = []
lons = []
mags = []
qual = []
stressdrops = []
cluster = []

for line in lines:
    dat = line.strip().split(',')
    lons.append(float(dat[2]))
    lats.append(float(dat[3]))
    mags.append(float(dat[8]))
    qual.append(float(dat[-2]))
    stressdrops.append(float(dat[10]))
    cluster.append(int(float(dat[-1])))

logstress = log10(array(stressdrops))
cluster = array(cluster)
lons = array(lons)
lats = array(lats)
qual = array(qual)

unique_clusters = unique(cluster)

# map clusters
#cols = get_mpl2_colourlist()
clust_stats = 'Cluster,Region,log10 SD +- STD,N\n'

# get national average
meanlogstress = nanmean(logstress)
stdlogstress = nanstd(logstress)
clust_stats += '0,National,' + str('%0.2f' % meanlogstress) + ' +- ' \
                   + str('%0.2f' % stdlogstress) + ',' + str(len(logstress)) + '\n'


syms = ['o', 's', 'd', 'H', '^', 'p', 'o', 's', 'd', 'H', '^', 'p']
for i, uc in enumerate(unique_clusters):
    idx = where((cluster == uc) & (qual == 1))[0]
    x, y = m(lons[idx], lats[idx])
    
    if i == 1 or i == 7:
        ms = 8
    else:
        ms = 9
    m.plot(x, y, ls='none', marker=syms[i], ms=ms, c=cols[i], label='Cluster '+str(i+1))
    
    meanlogstress = nanmean(logstress[idx])
    stdlogstress = nanstd(logstress[idx])
    print(uc+1, meanlogstress, stdlogstress, len(logstress[idx]))
    
    clust_stats += str(uc+1) + ',,' + str('%0.3f' % meanlogstress) + ' +- ' \
                   + str('%0.3f' % stdlogstress) + ',' + str(len(logstress[idx])) + '\n'

plt.legend(loc=3, numpoints=1, fontsize=14, ncol=2)

#########################################################################################
# finish
'''
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Voronoi.html
https://www.geeksforgeeks.org/machine-learning/plotting-boundaries-of-cluster-zone-with-scikit-learn/
'''
	
##########################################################################################

f = open('cluster_stats.csv','w')
f.write(clust_stats)
f.close()

plt.savefig('figures/mapped_event_clusters.png',format='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()