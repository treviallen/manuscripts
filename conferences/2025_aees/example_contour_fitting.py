from sys import argv
from scipy.interpolate import griddata
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, median, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor, round
from mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap, annotate_cities, get_field_data # drawshapepoly, labelpolygon, 
import shapefile
from shapely.geometry import Point, Polygon
import warnings
warnings.filterwarnings("ignore")
from misc_tools import get_mpl2_colourlist
#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

# set map resolution

region = argv[1]  # nat, vic, wa, etc
region = 'swwa'

bbox = '107.0/150.0/-43.0/-10.0'


'''
# get map bounds
llcrnrlat = minlat
urcrnrlat = maxlat
llcrnrlon = minlon
urcrnrlon = maxlon
lon_0 = mean([llcrnrlon, urcrnrlon])
lon_0 = 134.
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
'''

if region == 'nat':
    bbox = '106.0/152.5/-43.0/-7.0'
    ccbyX = 0.66
    ylabel = 6
    xlabel = 6
    legloc = 3
elif region == 'vic':
    bbox = '140.6/150.4/-39.5/-33.75'
    bbox = '143.2/150./-39.25/-34.25'
    insetLoc = 4
    insetSize = 0.034
    ccbyX = 0.185
    ylabel = 2
    xlabel = 2
    legloc = 2
elif region == 'wa':
    bbox = '112.5/129.5/-35.5/-13.5'
    insetLoc = 2
    insetSize = 0.1
    ccbyX = 0.585
elif region == 'swwa':
    bbox = '115.4/119/-32.5/-29.6'
    insetLoc = 2
    insetSize = 0.022
    ccbyX = 0.68
    ylabel = 1
    xlabel = 1
    legloc = 4
elif region == 'sa':
    bbox = '137.9/139.2/-35.5/-34.3'
    insetLoc = 2
    insetSize = 0.008
    ccbyX = 0.68
    ylabel = 0.5
    xlabel = 0.5
    legloc = 3
elif region == 'nt':
    bbox = '127/140/-28./-11'
    insetLoc = 3
    insetSize = 0.08
    ccbyX = 0.68
    ylabel = 3
    xlabel = 3
    legloc = 'center left'
elif region == 'nwwa':
    bbox = '118/129.25/-25./-13'
    insetLoc = 2
    insetSize = 0.065
    ccbyX = 0.68
    ylabel = 3
    xlabel = 3
    legloc = 4


bbox = bbox.split('/')
urcrnrlon = float(bbox[1])
llcrnrlon = float(bbox[0])
urcrnrlat = float(bbox[3])
llcrnrlat = float(bbox[2])

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
lat_0 = mean([lat_1, lat_2])

# set map
# Projection used for National Mapping
if region == 'nat':
    '''
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution=res,area_thresh=1000.)
    '''
    res = 'h' 
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, \
                resolution=res, area_thresh=1000., epsg=3112)
else:
    if region == 'sa' or region == 'swwa':
        res = 'f' 
        
        m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='merc',rsphere=6371200.,
                    resolution=res,area_thresh=0.) 
                    
        '''
        m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, epsg=3857)
        
        #http://server.arcgisonline.com/arcgis/rest/services
        
        #m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)
        '''
        
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

##############################################################################
# parse shapefile
##############################################################################

    else:
        res = 'h' 
        
        m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='merc',rsphere=6371200.,
                    resolution=res,area_thresh=100.) 
                    
        
        
irregular_shp = '2025_contours/1170.4_2025_PGA-0.033_0667_07_contours_08.shp' # original
sf = shapefile.Reader(irregular_shp)
irregular_shapes = sf.shapes()

#shpfile = '2025_contours/AS1170_4_fitted_contours.shp'
fitted_shp = '2025_contours/AS1170_4_fitted_modified_contours.shp'
sf = shapefile.Reader(fitted_shp)

fitted_shapes = sf.shapes()
levels = get_field_data(sf, 'LEVELS', 'float')

cols = get_mpl2_colourlist()

##############################################################################
# set some default values here
##############################################################################

fig = plt.figure(figsize=(19,12))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)

m.drawcoastlines(linewidth=0.5,color='k')
m.drawcountries(color='0.2')
m.drawstates(color='0.2')
m.drawmapboundary(fill_color='0.9')
m.fillcontinents(color='w', lake_color='0.9')

m.drawparallels(arange(-90.,90.,ylabel), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.5)
m.drawmeridians(arange(0.,360.,xlabel), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.5)
#m.arcgisimage(service='ESRI_Imagery_World_2D', verbose= True)
'''
# get mapsscale pos
widthkm = (urcrnrlon - llcrnrlon) * 110 * 0.25
ordmag = floor(log10(widthkm))
scalewith = round(widthkm / 10**ordmag) * 10**ordmag
slat = llcrnrlat + (urcrnrlat - llcrnrlat) * 0.06
if region == 'swwa' or region == 'nwwa':
    slon = llcrnrlon + (urcrnrlon - llcrnrlon) * 0.13
    m.drawmapscale(slon, slat, lon_0, lat_0, scalewith, fontsize = 14, barstyle='fancy')
elif region == 'nat':
    slon = llcrnrlon + (urcrnrlon - llcrnrlon) * 0.5
    m.drawmapscale(slon, slat, lon_0, lat_0, scalewith, fontsize = 14, barstyle='fancy')
elif region == 'vic':
    slon = llcrnrlon + (urcrnrlon - llcrnrlon) * 0.15
    m.drawmapscale(slon, slat, lon_0, lat_0, scalewith, fontsize = 14, barstyle='fancy')
else:
    slon = llcrnrlon + (urcrnrlon - llcrnrlon) * 0.83
    m.drawmapscale(slon, slat, lon_0, lat_0, scalewith, fontsize = 14, barstyle='fancy', zorder=100000)
'''    
##############################################################################
# now plot levels
##############################################################################

plt_levs = [0.08, 0.1]

for i, shapes in enumerate([irregular_shapes, fitted_shapes]):
    labeLegend = True
    for shape in shapes:
        x = []
        y = []
        
        p = 0
        parts = shape.parts
        parts.append(len(shape.points))
        
        for prt in parts[1:]:
            #print('part',prt
            while p < prt:
                x.append(shape.points[p][0])
                y.append(shape.points[p][1])
                p += 1
            
            # close polygon if not polyline
            if x[0] != x[-1] or y[0] != y[-1]:
                x.append(x[0])
                y.append(y[0])
        
            # plot each polygon                    
            xx, yy = m(x,y)
            #print(edgecolor)
            if labeLegend == True:
                if i == 0:
                    plt.plot(xx,yy, '--', c=cols[0], linewidth=2, label='Raw Contours')
                else:
                    plt.plot(xx,yy, '-', c=cols[1], linewidth=2, label='Fitted Contours')
                
                labeLegend = False
            else:
                if i == 0:
                    plt.plot(xx,yy, '--', c=cols[0], linewidth=2)
                else:
                    plt.plot(xx,yy, '-', c=cols[1], linewidth=2)

###########################################################################################
# annotate cities
###########################################################################################
capfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/shared/capitals_names.csv'

# annotate cities
if not region == 'nat':
    if region == 'swwa' or region == 'sa':
        numCities = 15
    elif region == 'nt' or region == 'nwwa':
        numCities = 15
    else:
        numCities = 40
    blacklist = ['Cunderdin', 'Dalwallinu', 'Mount Helena', 'Gabbadah']    
    annotate_cities(numCities, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, blacklist=blacklist)
    
    if region == 'vic':
        blacklist = ['Melbourne', 'Berwick', 'Traralgon', 'Sale', 'Bairnsdale', 'Mansfield', 'Glen Waverley']
        filterloc=[145, 150, -40, -37]
        annotate_cities(15, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        blacklist=blacklist, filterloc=filterloc)
                        
        filterloc=[145.5, 147, -38, -37.5]
        annotate_cities(4, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        blacklist=blacklist, filterloc=filterloc)
                        
    
    if region == 'swwa':
        blacklist = ['Northam', 'Wongan Hills', 'Merredin', 'Kellerberrin', 'Cunderdin', 'York', 'Mount Helena'\
                     'Mukinbudin', 'Castletown', 'Dalwallinu', 'Buntine', 'Bunji', 'Three Springs', 'South Kumminin']
        filterloc=[116.2, 119, -32.5, -29.5]
        annotate_cities(15, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        blacklist=blacklist, filterloc=filterloc)
                        
        '''
        filterloc=[118, 123, -35.5, -32]
        annotate_cities(6, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        blacklist=blacklist, filterloc=filterloc)
        filterloc=[115.5, 120, -30, -26.5]
        annotate_cities(6, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        blacklist=blacklist, filterloc=filterloc)
        '''
    if region == 'nwwa':
        filterloc=[120.5, 129, -25, -19.5]
        annotate_cities(4, plt, m, markersize=5, marker='o', markeredgewidth=0., markerfacecolor='k', fs=12.5, \
                        filterloc=filterloc)
    
else:
    import matplotlib.patheffects as PathEffects
    pe = [PathEffects.withStroke(linewidth=3, foreground="w")]
              
    llat = []
    llon = []
    locs = []
    textoffset = []
    
    # read data
    #capfile = path.join('..', '..', 'shared', 'capitals_names.csv')
    lines = open(capfile).readlines()
    for line in lines:
        llon.append(float(line.strip().split(',')[0]))
        llat.append(float(line.strip().split(',')[1]))
        locs.append(line.strip().split(',')[2])
        textoffset.append(float(line.strip().split(',')[3]))
    
    # plot locs on map
    x, y = m(array(llon), array(llat))
    plt.plot(x, y, 's', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1.5, markersize=8)
    
    # label cities
    for i, loc in enumerate(locs):
        if textoffset[i] == 0.:
            x, y = m(llon[i]-0.35, llat[i]+0.12)
            plt.text(x, y, loc, size=15, ha='right', va='bottom', weight='light', path_effects=pe)
        elif textoffset[i] == 1.:
            x, y = m(llon[i]+0.35, llat[i]+0.12)
            #plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
            plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
        elif textoffset[i] == 2.:
            x, y = m(llon[i]+0.3, llat[i]-0.3)
            plt.text(x, y, loc, size=15, ha='left', va='top', weight='light', path_effects=pe)
        elif textoffset[i] == 3.:
            x, y = m(llon[i]-0.3, llat[i]-0.2)
            plt.text(x, y, loc, size=15, ha='right', va='top', weight='light', path_effects=pe)
    

leg1 = plt.legend(loc=legloc, fontsize=16)
#plt.rcParams['legend.title_fontsize'] = 18
#leg1.get_title().set_ha("center")
leg1._legend_box.align = "center"
leg1.set_zorder(40000)

###############################################################################
# make map insetf
###############################################################################

if not region == 'nat':
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    if region == 'nt':
        x, y = m(136, -25)
        axins = zoomed_inset_axes(ax, insetSize, loc='center right') #, bbox_to_anchor=(x, y, 0.5, 0.5))
    else:
        axins = zoomed_inset_axes(ax, insetSize, loc=insetLoc)

    m2 = Basemap(projection='merc',\
        llcrnrlon=111,llcrnrlat=-45, \
        urcrnrlon=155,urcrnrlat=-9,\
        rsphere=6371200.,resolution='l',area_thresh=10000)
        
    m2.drawmapboundary(fill_color='0.9')
    m2.fillcontinents(color='w', lake_color='0.9') #, zorder=0)
    m2.drawcoastlines()
    m2.drawcountries()
    m2.drawstates()
    
    # fill main area
    xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
    yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
    x, y = m2(xv, yv)
    plt.plot(x, y, 'k-', lw=2)

plt.savefig(region+'_fitting_example_sat.png', dpi=300, format='png', bbox_inches='tight')
#plt.savefig(region+'_as1170_4_map.eps', dpi=300, format='eps', bbox_inches='tight')
plt.show()
