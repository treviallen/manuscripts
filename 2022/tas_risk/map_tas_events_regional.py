#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, remove_first_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn
import shapefile
from io_catalogues import parse_ga_event_query

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -39.3
llcrnrlat = -44.2
urcrnrlon = 150
llcrnrlon = 142.7
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='f',area_thresh=100.) #i

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.fillcontinents(color='w',lake_color='w')
#m.drawmapboundary(fill_color='w')
m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), labels=[0,0,0,1], fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print( 'Reading netCDF file...')
try:
    nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//au_indo_gebco_2020.nc')
except:
    nc = NetCDFFile('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =30. #gray
#zscale =30. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print( 'Getting colormap...')
# get colormap
try:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//lightgrey.cpt'
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//grey_fade_2.cpt'
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'

#cmap, zvals = cpt2colormap(cptfile, 30) # wiki
cmap, zvals = cpt2colormap(cptfile, 31) # grey
#cmap = remove_last_cmap_colour(cmap)

#cmap = ListedColormap
        
# make shading
print( 'Making map...')
ls = LightSource(azdeg = 135, altdeg = 10)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
#norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
norm = mpl.colors.Normalize(vmin=-2000/zscale, vmax=3500/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)

##########################################################################################
# add earthquakes
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

##########################################################################################
# set colours
##########################################################################################

ncols = 18

# get year range
minyear = 1850
maxyear = 2030
yearrng = float(round(maxyear - minyear))

#cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cs = (cmap(arange(ncols)))

##########################################################################################
# add earthquakes
##########################################################################################

evdict = parse_ga_event_query('earthquakes_export.csv')

lons = dictlist2array(evdict, 'lon')
lats = dictlist2array(evdict, 'lat')
mags = dictlist2array(evdict, 'mag')
years = dictlist2array(evdict, 'year')

# get zorder for plotting
sortidx = argsort(argsort(mags))
for i in range(0, len(mags)): #[0:100])):
    #get colour idx
    colidx = int(floor((ncols) * (years[i]-minyear+0.1) / yearrng))
    x, y = m(lons[i], lats[i])
    zo = sortidx[i] + 20
    plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
             markersize=(5. * mags[i] - 4), zorder=zo, alpha=0.8)
  
# make legend
legmag = [2.5, 4.5, 6.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', mfc='k', mec='w', mew=0.25, markersize=(5 * lm - 4), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

legtxt = ('$\mathregular{M_W}$ 2.5', '$\mathregular{M_W}$ 4.5', '$\mathregular{M_W}$ 6.5')
l = plt.legend(legh, legtxt, loc=3, numpoints=1, fontsize=10, title="Magnitude", labelspacing=0.75)
l.set_zorder(len(mags)+5)

##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,

ticks = arange(0, (ncols+1)/2)
years = float(minyear) + ticks*20
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.315,0.035,0.4,0.025]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.set_label('Year of Earthquake', rotation=0, fontsize=15, labelpad=5) # 

##########################################################################################
# finish
##########################################################################################

plt.savefig('tas_eq_map.png',fmt='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()