'''
script to map shakemap grid.xml data

Usage:
    python map_shakemap.py <path/to/grid.xml>
    
    e.g.:
        
    run map_shakemap.py ../../testing/grid.xml
'''

from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, mgrid, ogrid
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from mapping_tools import get_map_polygons, mask_outside_polygons, annotate_cities, return_city_list
from shakemap_tools import parse_dataxml
from sys import argv
from os import getcwd, mkdir, path
import shapefile
import matplotlib.gridspec as gridspec

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')
cwd = getcwd()

##########################################################################################
# read topo
##########################################################################################

mdiv = 500.

print 'Reading topo file...'
# use GMTED2010
#netcdffile = '../../../DATA/GEBCO/au_gebco.nc'
netcdffile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GMTED2010/50S120E_20101117_gmted_mea075.grd'
nc = NetCDFFile(netcdffile)

zscale = 30. #colour


'''
# if using GEBCO 30 arcsec
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]
'''

##########################################################################################
# parse shakemap grid
##########################################################################################

xmlPaths = ['M6.9_Adelaide_NF/current/products/grid.xml', \
            '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHA_18/Earthquake_Scenarios/Scenario_list1_contours/M6.9_Adelaide/current/products/grid.xml']
faultPath = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHA_18/Earthquake_Scenarios/Scenario_list1_contours/M6.9_Adelaide/current/Adelaide_fault.txt'

# loop through grids
fig = plt.figure(1, figsize=(18, 10))
gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.03, hspace=0.03)
pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']

for i, xmlPath in enumerate(xmlPaths):
    print(xmlPath)
    ax = fig.add_subplot(gs1[i])
    event, gridspecs, fields, xmldata = parse_dataxml(xmlPath)
    
    # set data fields
    xmllons = xmldata[:,0]
    xmllats = xmldata[:,1]
    
    # find MMI col
    keys = fields.keys()
    for key in keys:
        if fields[key] == 'mmi':
            mmicol = int(key[-1]) - 1
    
    mmi = xmldata[:,mmicol]
    
    ##########################################################################################
    # set up map
    ##########################################################################################
    # get map extents
    buff = 0.0
    urcrnrlat = gridspecs['lat_max'] - buff
    llcrnrlat = gridspecs['lat_min'] + buff
    urcrnrlon = gridspecs['lon_max'] - buff
    llcrnrlon = gridspecs['lon_min'] + buff
    lonspan = urcrnrlon - llcrnrlon
    
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    plt.tick_params(labelsize=16)
    
    m = Basemap(projection='lcc',\
                llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                rsphere=6371200.,resolution='h',area_thresh=200)
    '''
    epsg = 3107
    
    # note, you need change the epsg for different region, 
    #US is 4269, and you can google the region you want
    plt.figure(1, figsize = (10, 10))
    ax = plt.subplot(111)
    
    m = Basemap(projection='mill',llcrnrlon=llcrnrlon ,llcrnrlat=llcrnrlat,
        urcrnrlon=urcrnrlon ,urcrnrlat=urcrnrlat, resolution = 'h', epsg = epsg, area_thresh=200)
    '''    
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    #m.drawmapboundary(fill_color='blue', zorder=50)
    #m.fillcontinents(color='coral',lake_color='aqua')
    if lonspan <= 3.0:
        tickspace = 0.5
        scale_len = 50
        latoff = 0.13
        lonoff = 0.3
    elif lonspan <= 4.0:
        tickspace = 0.5
        scale_len = 100
        latoff = 0.2
        lonoff = 0.64
    elif lonspan > 4.0:
        tickspace = 1.0
        scale_len = 100
        latoff = 0.23
        lonoff = 0.64
    #m.drawparallels(arange(-90.,90.,tickspace), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
    #m.drawmeridians(arange(0.,360.,tickspace), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
    
    m.drawmapscale(llcrnrlon+lonoff, llcrnrlat+latoff, llcrnrlon+lonoff, llcrnrlat+latoff, scale_len, \
                   fontsize = 14, barstyle='fancy', zorder=100)
    
    ##########################################################################################
    # plot intensity grid from shakemap
    ##########################################################################################
    
    # set topo data
    data = nc.variables['z'][:] / zscale
    lons = nc.variables['x'][:]
    lats = nc.variables['y'][:]

    # transform topo to metres
    nx = int((m.xmax-m.xmin)/mdiv)+1
    ny = int((m.ymax-m.ymin)/mdiv)+1
    topodat = m.transform_scalar(data,lons,lats,nx,ny)
    
    print 'Resampling data...'
    N = 500j
    extent = (llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)
    
    xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    resampled = griddata(xmllons, xmllats, mmi, xs, ys, interp='linear')
    
    # get 1D lats and lons for map transform
    lons = ogrid[extent[0]:extent[1]:N]
    lats = ogrid[extent[2]:extent[3]:N]
    
    # nas interprets grids differently
    if cwd.startswith('/nas'):
        mmidat = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    else:
        mmidat = m.transform_scalar(resampled,lons,lats,nx,ny)
    
    print 'Getting colormap...'
    # get colormap
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/mi_pop.cpt'
    cmap, zvals = cpt2colormap(cptfile, 256)
    #cmap=plt.get_cmap('Spectral',256)
    
    # make shading
    print 'Making map...'
    ls = LightSource(azdeg = 180, altdeg = 0)
    norm = mpl.colors.Normalize(vmin=0, vmax=11)#myb
    rgb = ls.shade_rgb(cmap(norm(mmidat)), topodat, blend_mode='hsv', vert_exag=1) 
    im = m.imshow(rgb)
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    polys = get_map_polygons(m)
    
    #mask_outside_polygon(polys[1][::-1], ax=None)
    mask_outside_polygons(polys, 'lightskyblue', plt)
    
    # get lake ploygons
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
        polygons.append(poly)
        
    
    ##########################################################################################
    # add cities
    ##########################################################################################
    numCities=15
    clatList, clonList, cityList = return_city_list(numCities, m)
    txtoff = 0.03
    
    import matplotlib.patheffects as PathEffects
    path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
    
    for cla, clo, cn in zip(clatList, clonList, cityList):
        if cla+(txtoff+0.05) < m.urcrnrlat and cla+(txtoff+0.05) > m.llcrnrlat:
        
            x, y = m(clo, cla)
            plt.plot(x, y, 'o', markerfacecolor='k', \
                     markeredgecolor='k', markeredgewidth=1, \
                     markersize=6, zorder=11000)
            
            x, y = m(clo-txtoff, cla+txtoff) # changed for camden fig
            plt.text(x, y, cn, size=11, ha='right', weight='normal', path_effects=path_effects, zorder=11000)
        
    ##########################################################################################
    # annotate earthquake
    ##########################################################################################
    if i == 0:
        eqlat = event['lat']
        eqlon = event['lon']
        
        x, y = m(eqlon, eqlat)
        plt.plot(x, y, '*', markerfacecolor='None', markeredgecolor='k', markeredgewidth=2, markersize=18)
    
    ##########################################################################################
    # add fault
    ##########################################################################################
    if i == 1:
        lines = open(faultPath).readlines()
        flon = []
        flat = []
        for line in lines:
            dat = line.strip().split('\t')
            flon.append(float(dat[0]))
            flat.append(float(dat[1]))
        
        flon = array(flon)
        flat = array(flat)
        x, y = m(flon, flat)
        plt.plot(x, y, 'k-', lw=2)
     
    ##########################################################################################
    # add title
    ##########################################################################################
    
    '''
    titlestr = ' '.join(('MW',str(event['magnitude']), event['event_timestamp'])) + '\n' \
                         + event['event_description']
    plt.title(titlestr, fontsize=16)
    '''
    # plt letter
    xlim = ax.get_xlim()
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.96
    plt.text(xtxt, ytxt, pltlett[i], fontsize=20, va='top', ha='left',path_effects=path_effects)
                         
##########################################################################################
# make colorbar
##########################################################################################

# set colourbar
#plt.gcf().subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.34,0.11,0.33,0.03]) # setup colorbar axes.

cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=16)

##########################################################################################
# save file
##########################################################################################
    
pngFile = 'ade_fault_nofault_cmp.png'
plt.savefig(pngFile, format='png', bbox_inches='tight', dpi=150)
"""
##########################################################################################
# make shapefile of contour lines
##########################################################################################

# check to see if shapefile contours exists
if path.isdir('contours') == False:
    mkdir('contours')
    
# make list of levels - old levels array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.12, 0.18, 0.24])
allLevels = [arange(0.5, 11., 1.)]
 
levelNames = ['mmi_contours'] #, 'lev_0_01', 'lev_0_02', 'lev_0_05'


# resomple to smooth contours
N = 80j
xs2,ys2 = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
smoothed = griddata(xmllons, xmllats, mmi, xs2, ys2, interp='linear')

resampled = smoothed
xs = xs2
ys = ys2          

# loop thru levels
for levels, levelName in zip(allLevels, levelNames):
    
    # setup shapefile
    '''
    outshp = path.join('contours', '_'.join((levelNames[0].replace(' ','_'), key, \
                       levelName, 'contours.shp')))
    '''
    outshp = 'mmi_contours.shp'

    # set shapefile to write to
    w = shapefile.Writer(shapefile.POLYLINE)
    w.field('MMI','F', 5, 2)
        
    # have to re-contour using un-transformed lat/lons
    cs = plt.contour(xs, ys, resampled, levels, colors='k')
    
    plt.close(fig)
    
    # loop through contour levels
    for l, lev in enumerate(cs.levels):
        contours = cs.collections[l].get_paths()
        
        # now loop through multiple paths within level
        if len(contours) > 0:
            for cnt in contours:
                
                # add polyline to shapefile
                w.line(parts=[cnt.vertices], shapeType=shapefile.POLYLINE)
                
                # add level attribute
                w.record(lev)

    # now save area shapefile
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'wb')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()
"""
plt.show()
