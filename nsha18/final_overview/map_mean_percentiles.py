from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil, linspace
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, discrete_cmap
from mapping_tools import drawshapepoly, labelpolygon

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse hazard percentiles
##########################################################################################

# parse mean percentile
fracFolder = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/fractile_tables'
fracFile = 'fractile_table_PGA-0.1.csv'
#fracFile = 'fractile_table_PGA-0.02.csv'

fracPath = path.join(fracFolder, fracFile)

lines = open(fracPath).readlines()[1:]

lons = []
lats = []
pctl = []

for line in lines:
    dat = line.split(',')
    lons.append(float(dat[1]))
    lats.append(float(dat[2]))
    pctl.append(float(dat[-1]))

pctl = array(pctl)
prob = fracFile.split('-')[-1].strip('.csv')

##########################################################################################
#108/152/-44/-8
llcrnrlat = -45
urcrnrlat = -5.5
llcrnrlon = 106.
urcrnrlon = 153.
                
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(13,101))
ax = fig.add_subplot(111)
plt.tick_params(labelsize=8)

m = Basemap(projection='lcc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            lat_1=lat_1,lat_2=lat_2,lon_0=lon_0, \
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
m.drawstates()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt percentiles
##########################################################################################
if prob == '0.1':
    vmin = 30
    vmax = 80
    ncols = 10
else:
    vmin = 30
    vmax = 90
    ncols = 12
ticks = range(vmin, vmax+1, 5)

cmap = discrete_cmap(ncols, base_cmap=plt.cm.RdBu_r)

x, y = m(array(lons), array(lats))
m.scatter(x, y, c=pctl, s=80, cmap=cmap, vmin=vmin, vmax=vmax, alpha=1.)
#plt.title('Risk Coefficient at Sa '+pertxt+' s', fontsize=14)
'''
cax = fig.add_axes([0.915,0.25,0.02,0.5])
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cb = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm)
#cb.set_ticks(ticks)
cb.ax.tick_params(labelsize=14)
#cb = plt.colorbar(cax)
cb.set_label('Percentile at Mean Hazard', rotation=270, fontsize=18, labelpad=22)
'''
pngFile = 'fractile_map_' + fracFile.split('_')[-1].strip('.csv')+'.png'

plt.savefig(pngFile, format='png', bbox_inches='tight') #, dpi=200)
plt.show()

