#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort, nan
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, remove_first_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn, make_street_map
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
import shapefile
#from io_catalogues import parse_ga_event_query
from tools.mfd_tools import parse_hmtk_cat
from obspy import UTCDateTime

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -6.
llcrnrlat = -41.5
urcrnrlon = 153
llcrnrlon = 107
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

#m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)

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
# add earthquakes
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

##########################################################################################
# set colours
##########################################################################################

bounds = array([1, 3, 5, 10, 15, 20, 30, 1000])
cmap = plt.get_cmap('viridis_r', len(bounds))
cs = (cmap(arange(len(bounds))))

##########################################################################################
# add earthquakes
##########################################################################################

recs = pickle.load(open('../../../2023/au_stress_drop/fft_data.pkl', 'rb' ))

# convert mags to MW
for i, rec in enumerate(recs):
    if rec['magType'].lower().startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].lower().startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(rec['ev']) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.07
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag']) 
        
    # get chan details
    if len(rec['channels']) > 0:
        chan = rec['channels'][0]
        snr = rec[chan]['sn_ratio'][48] # 1 Hz  
        
        recs[i]['snr'] = snr  
        
    else:
        recs[i]['snr'] = nan

lons = dictlist2array(recs, 'eqlo')
lats = dictlist2array(recs, 'eqla')
mags = dictlist2array(recs, 'mag')
dts = dictlist2array(recs, 'ev')
snrs = dictlist2array(recs, 'snr')

# get unique event data
udts = unique(dts)

ulats = []
ulons = []
umags = []
ucnt = []
for u in udts:
    idx = where((dts == u) & (snrs >= 10))[0] # make sure exceeds snr
    if len(
    ucnt.append(len(idx))
    ulats.append(recs[idx[0]]['eqla'])
    ulons.append(recs[idx[0]]['eqlo'])
    umags.append(recs[idx[0]]['mag'])
  
# get zorder for plotting
sortidx = argsort(argsort(umags))
for i in range(0, len(umags)): #[0:100])):
    #get colour idx
    cidx = where(bounds >= ucnt[i])[0][0] # take first idx
    #col= tuple(cs[cidx][:-1])
    
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

ticks = arange(0, (len(bounds)+1))
years = float(minyear) + ticks*20
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.17,0.16,0.3,0.025]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.ax.tick_params(color='w', labelcolor='w')
cb.set_ticks(bounds)
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