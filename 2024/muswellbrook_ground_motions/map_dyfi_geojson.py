# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:32:05 2019

@author: u56903
"""
from mapping_tools import make_street_map, get_map_polygons, mask_outside_polygons, annotate_cities
from numpy import mean, percentile, arange
from mpl_toolkits.basemap import Basemap
from matplotlib import colors, colorbar
from gmt_tools import cpt2colormap, remove_last_cmap_colour
import matplotlib.pyplot as plt
from os import getcwd
from sys import argv
import matplotlib as mpl
mpl.style.use('classic')

import json
from numpy import array

#jsonFilePath = 'aggregation/canberra_1_1.geojson'
#jsonFilePath = 'events/20190714.Broome/felt_reports_5km.geojson'

jsonFilePath = argv[1]

with open(jsonFilePath) as f:
    data = json.load(f)

dyfi_dict = []
for feature in data['features']:
    tmp = {'geomerty':feature['geometry']['coordinates'][0],
           'centroid':feature['properties']['center']['coordinates'],
           'intensity':feature['properties']['intensityFine'],
           'nresp':feature['properties']['nresp']}
    
    # append to list
    dyfi_dict.append(tmp)
        
###############################################################################
# set eq details
###############################################################################

# Banda Sea
eqla = -32.36
eqlo = 150.85
degrng = 2.5
latoff = 0.
lonoff = 0.
mstr = '4.7'
place = 'Muswellbrook, NSW'
evid = '2024-08-23'

##########################################################################################
# set up map
##########################################################################################
'''
plt, m = make_street_map(eqla, eqlo, service='ESRI_Imagery_World_2D', ll_buffer = degrng, \
             xpixels = 1500, plt_inset = True, inset_state = 'nsw', inset_loc = 3, \
             plt_marker = True, marker='*', ms = 14, mew = 1., mfc = 'none', mec='r')
'''

##########################################################################################
# set up map
##########################################################################################

# set map on fly
if degrng < 0.5:
    parSpace = 0.2
    merSpace = 0.2
    scaleLength = 15
    res = 'f'
elif degrng < 1:
    parSpace = 0.4
    merSpace = 0.4
    scaleLength = 30
    res = 'f'
elif degrng >= 1 and degrng <= 2:
    parSpace = 1.
    merSpace = 1.
    res = 'h'
    scaleLength = 50
elif degrng > 2 and degrng < 5:
    parSpace = 1.
    merSpace = 2.
    res = 'h'
    scaleLength = 100
elif degrng >= 5.:
    parSpace = 2.
    merSpace = 2.
    res = 'i'
    scaleLength = 400
    
# make bounds on the fly - adjust to make square-ish
urcrnrlat = eqla + degrng*0.9 + latoff
llcrnrlat = eqla - degrng*0.9 + latoff
urcrnrlon = eqlo + degrng + lonoff
llcrnrlon = eqlo - degrng + lonoff

# set up figure
fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

# set projection
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

from mpl_toolkits.basemap import Basemap

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, epsg=3857) #3112)
m.drawparallels(arange(-90.,90.,0.5), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
m.drawmeridians(arange(0.,360.,0.5), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
#http://server.arcgisonline.com/arcgis/rest/services

m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)
'''
m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution=res,area_thresh=5000)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()

m.fillcontinents(color='0.9',lake_color='lightskyblue')
m.drawparallels(arange(-90.,90.,parSpace), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,merSpace), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
'''
# get scale position and length
slon = m.llcrnrlon + 0.12*(m.urcrnrlon - m.llcrnrlon)
slat = m.llcrnrlat + 0.06*(m.urcrnrlat - m.llcrnrlat)
slon0 = mean([m.urcrnrlon, m.llcrnrlon])
slat0 = mean([m.urcrnrlat, m.llcrnrlat])

m.drawmapscale(slon, slat, slon0, slat0, scaleLength, fontsize = 12, barstyle='fancy', zorder=10000)

##########################################################################################
# get land & lake polygons for masking
##########################################################################################
'''
polys = get_map_polygons(m)

mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons and mask
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
    polygons.append(poly)
'''
##########################################################################################
# get cmap
##########################################################################################
# get colormap
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/mi_pop.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//mi_pop.cpt'
ncols = 10
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

##########################################################################################
# plt dyfi
##########################################################################################

# loop thru grid cells values and plt MMI
nresp = 0
min_resp = 0
for dyfi in dyfi_dict:        
    # add to list greater than minObs
    if dyfi['nresp'] > min_resp:
        
        # now plot
        pltx = array(dyfi['geomerty'])[:,0]
        plty = array(dyfi['geomerty'])[:,1]
        
        x, y = m(pltx, plty)
        colidx = int(round(dyfi['intensity']))-1
        c= tuple(cs[colidx][:-1])
        plt.fill(x, y, fc=c, ec='none', lw=0.25, zorder=100, alpha=0.8)
        
    nresp += dyfi['nresp']
        
##########################################################################################
# annotate
##########################################################################################

# plt earthquake epicentre
x, y = m(eqlo, eqla)
plt.plot(x, y, '*', markersize=14, markerfacecolor='None', mec='r', mew=1.5, zorder=100000)

# annotate cities
numCities = 20
annotate_cities(numCities, plt, m, blacklist=['Narrabri'], fs=11)

# add n responses txt
lonrng = m.urcrnrlon - m.llcrnrlon
latrng = m.urcrnrlat - m.llcrnrlat
plttxt = 'Number of Responses = '+str(nresp)
x, y = m(m.llcrnrlon+0.02*lonrng, m.llcrnrlat+0.98*latrng)

props = dict(boxstyle='round', facecolor='w', alpha=1)
plt.text(x, y, plttxt, size=16, ha='left', va='top', bbox=props, zorder=10000)

##########################################################################################
# make colorbar
##########################################################################################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.12)
cax = fig.add_axes([0.34,0.06,0.33,0.03]) # setup colorbar axes.

norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', alpha=0.8)

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=16)

plt.savefig('_'.join((evid.replace('-',''), place.split(',')[0].replace(' ', '_'), \
            'gridded_mmi_data_gridded.png')), format='png', bbox_inches='tight', dpi=300)
plt.show()