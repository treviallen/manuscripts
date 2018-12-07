#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system, getcwd
#from obspy.imaging.beachball import Beach
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour
from mapping_tools import get_field_data
import shapefile
from sys import argv

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')


def get_basemap(llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat, res, fig, useGEBCO):
    '''
    useGEBCO = True: use GEBCO data
    useGEBCO = False: use shadedrelief
    '''
    from mpl_toolkits.basemap import Basemap
    import matplotlib as mpl
    from matplotlib.colors import LightSource
    from netCDF4 import Dataset as NetCDFFile
    from os import getcwd

    cwd = getcwd()
    
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    plt.tick_params(labelsize=16)
    ax = fig.add_subplot(111)
    
    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
                llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                rsphere=6371200.,resolution=res,area_thresh=100.)
    
    # draw coastlines, state and country boundaries, edge of map.
    if useGEBCO == False:
        m.shadedrelief()
        
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,2.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
    #m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)
    
    ##########################################################################################
    # plot gebco
    ##########################################################################################
    if useGEBCO == True:
        print 'Reading netCDF file...'
        if cwd.startswith('/nas'):
            nc = NetCDFFile('//nas//gemd//ehp//georisk_earthquake//hazard//DATA//GEBCO//au_gebco.nc')
            cptfile = '//nas//gemd//ehp//georisk_earthquake//hazard//DATA//cpt//mby_topo-bath.cpt'
            #cptfile = '//nas//gemd//ehp//georisk_earthquake//hazard//DATA//cpt//gray.cpt'
            #cptfile = '//nas//gemd//ehp//georisk_earthquake//hazard//DATA//cpt//mby_topo-bath_mod.cpt'
        else:
            nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
            cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mby_topo-bath.cpt'
            cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
        
        zscale =20. #gray
        zscale =30. #colour
        data = nc.variables['elevation'][:] / zscale
        lons = nc.variables['lon'][:]
        lats = nc.variables['lat'][:]
        
        # transform to metres
        nx = int((m.xmax-m.xmin)/500.)+1
        ny = int((m.ymax-m.ymin)/500.)+1
        
        topodat = m.transform_scalar(data,lons,lats,nx,ny)
        
        print 'Getting colormap...'
        # get colormap
        
        #cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gray.cpt'
        cmap, zvals = cpt2colormap(cptfile, 30)
        cmap = remove_last_cmap_colour(cmap)
        #cmap = remove_last_cmap_colour(cmap)
        
        # make shading
        print 'Making map...'
        ls = LightSource(azdeg = 180, altdeg = 5)
        #norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
        norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
        rgb = ls.shade(topodat, cmap=cmap, norm=norm)
        im = m.imshow(rgb)

    
    return m, ax


reg = argv[1]

plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
if reg == 'sea':
    #108/152/-44/-8
    urcrnrlat = -34.
    llcrnrlat = -39.2
    urcrnrlon = 151.
    llcrnrlon = 143.8
    legloc = 2
    maploc = 4
    
elif reg == 'ade':
    #108/152/-44/-8
    urcrnrlat = -31.5
    llcrnrlat = -36.25
    urcrnrlon = 140.5
    llcrnrlon = 135.5
    legloc = 3
    maploc = 2

fig = plt.figure(figsize=(14,10))
useGEBCO = True

m, ax = get_basemap(llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat, 'h', fig, useGEBCO)

##########################################################################################
# parse dams
##########################################################################################
cwd = getcwd()
'''
if cwd.startswith('/nas'):
    damshp = '/nas/gemd/ehp/georisk_earthquake/hazard/DATA/GIS/Water_Storage/aus_dam_gt_100k.shp'
else:
    damshp = '/Users/tallen/Documents/DATA/GIS/ArcGIS/Water_Storage/aus_dam_gt_100k.shp'
    
sf = shapefile.Reader(damshp)
shapes = sf.shapes()

# map dams
dms = 12
for shape in shapes:
    x, y = m(shape.points[0][0], shape.points[0][1])
    plt.plot(x, y, 'h', mfc='b', mec='w', mew=0.75, markersize=dms)

# replot last for label
plt.plot(x, y, 'h', mfc='b', mec='w', mew=0.75, markersize=dms, label='Large Dams')
'''
##########################################################################################
# add simple faults
##########################################################################################
if cwd.startswith('/nas'):
    nfsmshp = '//nas//gemd//ehp//georisk_earthquake//neotectonic//Seismicity_Scenario_models//Hazard Map working 2018//ARCGIS//FSM lines//FSD_simple_faults.shp'
else:
    nfsmshp = '/Users/tallen/Documents/Geoscience_Australia/Neotectonics/FSD_simple_faults.shp'
    
sf = shapefile.Reader(nfsmshp)
shapes = sf.shapes()
records = sf.records()
lt_rates = get_field_data(sf, 'SL_RT_LT', 'float') # long term slip rate
lt_rates = array(lt_rates)
st_rates = get_field_data(sf, 'SL_RT_ST', 'float') # short term slip rate
st_rates = array(st_rates)

if cwd.startswith('/nas'):
    cptfile = '//nas//gemd//ehp//georisk_earthquake//hazard//DATA//cpt//temperature.cpt'
else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//temperature.cpt'
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//humidity.cpt'
    
ncols = 17
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False) # a
cmap = remove_last_cmap_colour(cmap)
cmap = remove_last_cmap_colour(cmap)
#cmap = remove_last_cmap_colour(cmap)

#cmap =  mpl.cm.jet
cs = (cmap(arange(cmap.N)))

minslip = 0.
maxslip = 160.
slip_diff = maxslip - minslip

sortidx = argsort(argsort(lt_rates))

i = 0
for shape, ltr in zip(shapes, lt_rates):
    lons = []
    lats = []
    for xy in shape.points:
        lons.append(xy[0])
        lats.append(xy[1])
    
    x, y = m(lons, lats)
    
    if ltr > maxslip:
        ltr = maxslip
    
    # set colour by 
    colidx = int(round((cmap.N-1) * (ltr-minslip) / slip_diff))
    plt.plot(x, y, '-', c=cs[colidx], lw=2.5, zorder=sortidx[i])
    i += 1


##########################################################################################
# add cities
##########################################################################################

capfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/capitals_names.csv'

# read data
lines = open(capfile).readlines()
for line in lines:
    clon = float(line.strip().split(',')[0])
    clat = float(line.strip().split(',')[1])
    loc = line.strip().split(',')[2]
    textoffset = float(line.strip().split(',')[3])
    
    if clon > llcrnrlon and clon < urcrnrlon \
       and clat > llcrnrlat and clat < urcrnrlat:
       
        x, y = m(clon, clat)
        plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', \
                 markeredgewidth=1.5, markersize=10, label='Capital Cities')
                 
        if loc == 'Canberra' or loc == 'Adelaide':
            x, y = m(clon-0.1, clat+0.07)
            plt.text(x, y, loc, size=20, ha='right', va='bottom', weight='normal')
        elif loc == 'Melbourne':
            x, y = m(clon+0.1, clat+0.07)
            #plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
            plt.text(x, y, loc, size=20, ha='left', va='bottom', weight='normal')
        elif textoffset == 2.:
            x, y = m(clon+0.1, clat-0.3)
            plt.text(x, y, loc, size=20, ha='left', va='top', weight='normal')
        elif textoffset == 3.:
            x, y = m(clon-0.1, clat-0.2)
            plt.text(x, y, loc, size=20, ha='right', va='top', weight='normal')
            

# plt letter
xlim = ax.get_xlim()
xtxt = xlim[1] * 0.02
ylim = ax.get_ylim()
if reg == 'sea':
    plt.text(xtxt, ytxt, '(b)', fontsize=30, va='top', ha='left') 
    ytxt = ylim[1] * 0.96

else:
    plt.text(xtxt, ytxt, '(a)', fontsize=30, va='bottom', ha='left')
    ytxt = ylim[1] * 0.04

#plt.legend(numpoints=1, loc=legloc, fontsize=20)

##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.03, loc=maploc)

m2 = Basemap(projection='merc',\
            llcrnrlon=111,llcrnrlat=-45, \
            urcrnrlon=156,urcrnrlat=-9,\
            rsphere=6371200.,resolution='c',area_thresh=1000)
            
#m2 = Basemap(llcrnrlon=-20,llcrnrlat=3,urcrnrlon=0,urcrnrlat=18, ax=axins)
m2.fillcontinents(color='w', lake_color='0.8', zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()
#m2.drawmapboundary(fill_color='0.8')
m2.drawlsmask(land_color='w', ocean_color='0.8')

# fill main area
xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'k-', lw=1.5)

##########################################################################################
# label states
##########################################################################################

state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -29.5, -23.0, -32.5, -37.1, -42.]
slon = [122, 133.5, 135.0, 144.5, 146.5, 143.6, 147.0]
for i, st in enumerate(state):
    x, y = m2(slon[i], slat[i])
    plt.text(x, y, st, size=11, horizontalalignment='center', verticalalignment='center', weight='normal')

##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,
#ticks = arange(0, 19, 2)
#years = float(minyear) + ticks*10
#labels = [str('%0.0f' % x) for x in years]
if reg == 'sea':
    # normalise
    norm = mpl.colors.Normalize(vmin=minslip, vmax=maxslip)
    
    cax = fig.add_axes([0.85,0.3,0.02,0.4]) # setup colorbar axes.
    cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=1., norm=norm) # \
    cb.ax.tick_params(labelsize=15)
    #cb.set_ticks(years)
    #cb.set_ticklabels(labels)
    cb.set_label('Long-Term Slip Rate (m/Ma)', rotation=270, labelpad=25, fontsize=17)


plt.savefig(reg+'_faults_map.png', format='png', bbox_inches='tight', dpi=100)
plt.show()
