def return_dyfi_data(jsonfile):
    import json


    with open(jsonfile) as f:
        data = json.load(f)

    dyfi_dict = []
    for feature in data['features']:
        tmp = {'geomerty':feature['geometry']['coordinates'][0],
               'centroid':feature['properties']['center']['coordinates'],
               'intensity':feature['properties']['intensityFine'],
               'nresp':feature['properties']['nresp']}
        
        # append to list
        dyfi_dict.append(tmp)
        
    return dyfi_dict

###############################################################################
# start main code here
###############################################################################

from io_catalogues import parse_ga_event_query
from misc_tools import listdir_extension
from numpy import array, isnan, nan, unique, mean, percentile, arange, where
from os import path, getcwd
from mapping_tools import make_street_map, get_map_polygons, mask_outside_polygons, annotate_cities
from datetime import datetime

###############################################################################
# parse catalogue
###############################################################################

gacsv = 'earthquakes_export_edit.csv'
events = parse_ga_event_query(gacsv)

###############################################################################
# grab 10-km aggregated
###############################################################################
json_folder = 'all_felt_reports_aggregated'

# getting aggregated files
all_geojson = listdir_extension(json_folder, 'geojson')

# now just keep 10-km agg
geojson_10 = []
for jsonfile in all_geojson:
    #if jsonfile.endswith('10km_filtered.geojson'):
    if jsonfile.endswith('5km.geojson'):
        geojson_10.append(jsonfile)
        
###############################################################################
# loop through events and match json files
###############################################################################

first_date = datetime(2030, 12, 12)
for i, event in enumerate(events):
    events[i]['dyfi'] = nan # initialise
    
    for jsonfile in geojson_10:
        if jsonfile.startswith(event['event_id']):
            dyfi_dict = return_dyfi_data(path.join(json_folder, jsonfile))
            
            events[i]['dyfi'] = dyfi_dict
            
            if event['datetime'] < first_date:
                first_date = event['datetime']
                
print('First Date', first_date)

###############################################################################
# loop through events get unique centroids
###############################################################################

cent_lonlist = []
cent_latlist = []
cent_mmi = []
cent_nresp = []
total_nresp = 0
ev_count = 0

centroids = []
for i, event in enumerate(events):
    if isinstance(event['dyfi'],list):
       # loop through dyfi recs
       for rec in event['dyfi']:
           cent_lonlist.append(rec['centroid'][0])
           cent_latlist.append(rec['centroid'][1])
           cent_mmi.append(rec['intensity'])
           cent_nresp.append(rec['nresp'])
           centroids.append(array(rec['centroid']))
           
           # add nresp
           total_nresp += rec['nresp']
           len_recs = len(event['dyfi'])
       
       if len_recs > 2:
           ev_count += 1

crash = crash
cent_lonlist = array(cent_lonlist)
cent_latlist = array(cent_latlist)
cent_mmi = array(cent_mmi)
cent_nresp = array(cent_nresp)

'''
# get unique values
unique_lons = unique(array(cent_lonlist))
unique_lats = unique(array(cent_latlist))
unique_cent = []

i = 0
for clo, cla in zip(cent_lonlist, cent_latlist):
    cent_match = False
    if i == 0:
        unique_cent.append([clo, cla])
    else:
        for ucent in unique_cent:
            if ucent[0] == clo and ucent[1] == cla:
                cent_match = True
                
        if cent_match == False:
            if clo > 105 and clo < 165 and cla > -50 and cla < -5:
                unique_cent.append([clo, cla])
        
    i += 1
'''
'''
# get max mmi vals
max_mmi_10k = []
for ucent in unique_cent:
    max_mmi = 0.
    for cent in centroids:
        if cent == ucent:
           if cm > max_mmi:
               max_mmi = cm
    max_mmi_10k.append(max_mmi)
'''

###############################################################################
# start mapping
###############################################################################

from mpl_toolkits.basemap import Basemap
from matplotlib import colors, colorbar
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from mapping_tools import distance, reckon
import matplotlib.pyplot as plt
import matplotlib as mpl     
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mpl.style.use('classic')

##########################################################################################
# get cmap
##########################################################################################
# get colormap
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/qual-dark-06.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//qual-dark-06.cpt'
ncols = 6
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

##########################################################################################
# make map
##########################################################################################

bbox = '110.5/156.5/-44.5/-9.5'

bbox = bbox.split('/')
minlon = float(bbox[0])
maxlon = float(bbox[1])
minlat = float(bbox[2])
maxlat = float(bbox[3])
mbuff = 1.

# make bounds on the fly - adjust to make square-ish
urcrnrlat = maxlat
llcrnrlat = minlat
urcrnrlon = maxlon
llcrnrlon = minlon

# set up figure
fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

# set projection
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=1000)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()

m.fillcontinents(color='0.9',lake_color='lightskyblue')
m.drawparallels(arange(-90.,90.,6), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,10), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)

##########################################################################################
# get land & lake polygons for masking
##########################################################################################

polys = get_map_polygons(m)

mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons and mask
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
    polygons.append(poly)

###############################################################################
# make 50-km grid the ugly way
###############################################################################

cell_size = 35.
min_nresp = 2

grd_lons = []
grd_lats = []
grd_geom = []
grd_mmi = []
grd_count = []
lon_grd = minlon

while lon_grd < maxlon:
    lat_grd = minlat
    nxt_lon, tmp_lat = reckon(lat_grd, lon_grd, cell_size*0.95, 90)
    
    while lat_grd < maxlat:
       tmp_lon, nxt_lat = reckon(lat_grd, lon_grd, cell_size, 0)
       
       # fill grids
       geom = [[lon_grd, lat_grd], [lon_grd, nxt_lat], [nxt_lon, nxt_lat], [nxt_lon, lat_grd], [lon_grd, lat_grd]]
       grd_geom.append(geom)
       
       # check if mmi data exists
       idx = where((cent_lonlist >= lon_grd) & (cent_lonlist < nxt_lon) \
                   & (cent_latlist >= lat_grd) & (cent_latlist < nxt_lat) & (cent_nresp >= min_nresp))[0]
                   
       if len(idx) > 0:
           grd_mmi.append(max(cent_mmi[idx]))
           grd_count.append(len(idx))
           #!!! need to add datetimes and hstack then get unique!!!
       else:
           grd_mmi.append(nan)
           grd_count.append(nan)
           
       lat_grd = nxt_lat
       
    lon_grd = nxt_lon

# make coulour idx for counts
bounds = array([1, 3, 5, 10, 20, 100])
for geom, gcnt in zip(grd_geom, grd_count):
    if not isnan(gcnt):
    
        # map grid
        pltx = array(geom)[:,0]
        plty = array(geom)[:,1]
        
        x, y = m(pltx, plty) 
        
        cidx = where(bounds <= gcnt)[0][-1] # take last idx
        c= tuple(cs[cidx][:-1])
        plt.fill(x, y, fc=c, ec='0.45', lw=0.25, zorder=100)
        
       # plt.plot(x, y, '-', c='0.5', lw=0.5)   
'''
# map centroids
ulos = array(unique_cent)[:,0]
ulas = array(unique_cent)[:,1]
x, y = m(ulos, ulas)
plt.plot(x, y, 'b+')
'''
##########################################################################################
# annotate
##########################################################################################

# add n responses txt
lonrng = m.urcrnrlon - m.llcrnrlon
latrng = m.urcrnrlat - m.llcrnrlat
plttxt = 'Number of Responses = '+str(total_nresp)
x, y = m(m.llcrnrlon+0.02*lonrng, m.llcrnrlat+0.98*latrng)

props = dict(boxstyle='round', facecolor='w', alpha=1)
plt.text(x, y, plttxt, size=16, ha='left', va='top', bbox=props, zorder=10000)

##########################################################################################
# add colorbar
##########################################################################################

axins = inset_axes(ax,
                   width="50%",  # width = 50% of parent_bbox width
                   height="4%",  # height : 5%   
                   loc=3,
                   bbox_to_anchor=(0.01,0.062,1,1), bbox_transform=ax.transAxes)
                    
norm = mpl.colors.Normalize(vmin=0, vmax=6)#myb
cb = colorbar.ColorbarBase(axins, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = arange(0.5,6)
cnt_rng = ['1', '2-3', '4-5', '6-10', '11-20', '21+']
cb.set_ticks(ticks)
cb.set_ticklabels(cnt_rng)

titlestr = 'Number of Felt Events'
cb.set_label(titlestr, fontsize=16)


plt.savefig('map_nat_agg_mmi_count.png', fmt='png', bbox_inches='tight')
plt.show()