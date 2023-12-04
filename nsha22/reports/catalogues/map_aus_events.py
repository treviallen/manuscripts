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
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn, make_street_map
import shapefile
#from io_catalogues import parse_ga_event_query
from tools.mfd_tools import parse_hmtk_cat
from obspy import UTCDateTime

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -6.
llcrnrlat = -42
urcrnrlon = 156.5
llcrnrlon = 101
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)
'''
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
m.drawparallels(arange(-90.,90.,1.),fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)
'''

##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
ll_buffer = (urcrnrlon - llcrnrlon) / 2.
'''
plt, m, ax = make_street_map(clat, clon, service='World_Imagery', \
             xpixels = 150, plt_inset = False, plt_marker = False, inset_loc=4, inset_multiplier=0.1)
'''
from mpl_toolkits.basemap import Basemap

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, epsg=3112)
#http://server.arcgisonline.com/arcgis/rest/services

m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)

'''
    Map Services:
        ESRI_Imagery_World_2D
        ESRI_StreetMap_World_2D
        NatGeo_World_Map
        NGS_Topo_US_2D
        Ocean_Basemap
        USA_Topo_Maps
        World_Imagery
        World_Physical_Map
        World_Shaded_Relief
        World_Street_Map
        World_Terrain_Base
        World_Topo_Map
        
ESRI_Imagery_World_2D - good sat image
NatGeo_World_Map - nice road map - hard to read
'''

##########################################################################################
# plot gebco
##########################################################################################
"""
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
"""
##########################################################################################
# add earthquakes
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

##########################################################################################
# set colours
##########################################################################################

ncols = 15

# get year range
minyear = 1880
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

hmtk_csv = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/merged_NSHA23-ISCGEM_hmtk.csv'

evdict = []
lines = open(hmtk_csv).readlines()[1:]
	
for line in lines:
    dat = line.strip().split(',')
    tmp = {'lon':float(dat[7]), 'lat':float(dat[8]), 'mag':float(dat[10]), 'year':float(dat[1])}
    	
    evdict.append(tmp) 

lons = dictlist2array(evdict, 'lon')
lats = dictlist2array(evdict, 'lat')
mags = dictlist2array(evdict, 'mag')
years = dictlist2array(evdict, 'year')

# get zorder for plotting
sortidx = argsort(argsort(mags))
for i in range(0, len(mags)): #[0:100])):
    #get colour idx
    if mags[i] >= 3.:
        colidx = int(floor((ncols) * (years[i]-minyear+0.1) / yearrng))
        x, y = m(lons[i], lats[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=(2. * mags[i] - 1), zorder=zo, alpha=0.8)
 
# make legend
legmag = [3., 5., 7.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', mfc='w', mec='k', mew=0.25, markersize=(2 * lm - 1), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

legtxt = ('$\mathregular{M}$ 3.0', '$\mathregular{M}$ 5.0', '$\mathregular{M}$ 7.0')
l = plt.legend(legh, legtxt, loc=4, numpoints=1, fontsize=10, title="Magnitude", labelspacing=0.75, frameon=False)
for text in l.get_texts():
    text.set_color("w")
l.set_zorder(len(mags)+5)
'''
import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(weight='bold',
                                   style='normal')
'''
plt.setp(l.get_title(), color='white')

##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,

ticks = arange(0, (ncols+1)/2)
years = float(minyear) + ticks*20
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.245,0.16,0.3,0.025]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.ax.tick_params(color='w', labelcolor='w')
for spine in cb.ax.spines.values():
    spine.set_edgecolor('w')
cb.set_label('Year of Earthquake', rotation=0, fontsize=15, labelpad=5, color='w')
cb.outline.set_edgecolor('k') # 

##########################################################################################
# finish
##########################################################################################

plt.savefig('nat_eq_map_satellite.png',fmt='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()