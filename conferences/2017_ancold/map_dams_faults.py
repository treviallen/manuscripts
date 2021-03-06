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
                rsphere=6371200.,resolution=res,area_thresh=10.)
    
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
        cmap, zvals = cpt2colormap(cptfile, 256)
        cmap = remove_last_cmap_colour(cmap)
        #cmap = remove_last_cmap_colour(cmap)
        
        # make shading
        print 'Making map...'
        ls = LightSource(azdeg = 300, altdeg = 45)
        norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
        rgb = ls.shade(topodat, cmap=cmap, norm=norm) #norm=norm)
        im = m.imshow(rgb, alpha=0.4)
    
    return m, ax


plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
#108/152/-44/-8
urcrnrlat = -34.
llcrnrlat = -39.2
urcrnrlon = 151.
llcrnrlon = 143.8

fig = plt.figure(figsize=(14,10))
useGEBCO = True

m, ax = get_basemap(llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat, 'h', fig, useGEBCO)

##########################################################################################
# parse dams
##########################################################################################
cwd = getcwd()
if cwd.startswith('/nas'):
    damshp = '/nas/gemd/ehp/georisk_earthquake/hazard/DATA/GIS/Water_Storage/aus_dam_gt_100k.shp'
    
sf = shapefile.Reader(damshp)
shapes = sf.shapes()

# map dams
for shape in shapes:
    x, y = m(shape.points[0][0], shape.points[0][1])
    plt.plot(x, y, 'h', mfc='b', mec='w', mew=0.75, markersize=9)

# replot last for label
plt.plot(x, y, 'h', mfc='b', mec='w', mew=0.75, markersize=9, label='Large Dams')
##########################################################################################
# add simple faults
##########################################################################################

nfsmshp = '//nas//gemd//ehp//georisk_earthquake//neotectonic//Seismicity_Scenario_models//Hazard Map working 2018//ARCGIS//FSM lines//FSD_simple_faults.shp'
sf = shapefile.Reader(nfsmshp)
shapes = sf.shapes()
records = sf.records()
lt_rates = get_field_data(sf, 'SL_RT_LT', 'float') # long term slip rate
lt_rates = array(lt_rates)
st_rates = get_field_data(sf, 'SL_RT_ST', 'float') # short term slip rate
st_rates = array(st_rates)


if cwd.startswith('/nas'):
    cptfile = '//nas//gemd//ehp//georisk_earthquake//hazard//DATA//cpt//temperature.cpt'
ncols = 18
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)

#cmap =  mpl.cm.jet
cs = (cmap(arange(cmap.N)))

minslip = 0.
maxslip = 160.
slip_diff = maxslip - minslip

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
    plt.plot(x, y, '-', c=cs[colidx], lw=2.)


"""
ncols = 18

# get year range
minyear = 10*floor(cat.data['year'][0]/10.)
maxyear = cat.data['year'][-1]
yearrng = float(round(maxyear - minyear))
#cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/grayscale08.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:1], cs[2:], [1., 1., 1., 1.]))

# get zorder for plotting
sortidx = argsort(argsort(cat.data['magnitude']))
for i in range(0, len(cat.data['magnitude'])): #[0:100])):
    if cat.data['magnitude'][i] >= 3.0:
        #get colour idx
        year = cat.data['year'][i]
        colidx = int(round((ncols-1) * (year-minyear) / yearrng))
        x, y = m(cat.data['longitude'][i], cat.data['latitude'][i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', markerfacecolor=cs[colidx], markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + cat.data['magnitude'][i]*2.8, zorder=zo, alpha=0.8)
    
# make legend
legmag = [3., 5., 7.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-4 + lm*2.8), alpha=1., zorder=len(cat.data['magnitude'])+1, lw=2)
    legh.append(h[0])

l = plt.legend(legh, ('MW 3.0', 'MW 5.0', 'MW 7.0'), loc=1, numpoints=1)
l.set_zorder(len(cat.data['magnitude'])+5)
"""


##########################################################################################
# add cities
##########################################################################################
'''
clat = [-37.814, -42.882, -35.282]
clon = [144.964, 147.324, 149.129]
cloc = ['Melbourne', 'Hobart', 'Canberra']
'''
clat = [-37.814,  -35.282]
clon = [144.964,  149.129]
cloc = ['Melbourne', 'Canberra']

x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', \
         markeredgewidth=1.5, markersize=10, label='Capital Cities')

plt.legend(loc=2)
# label cities
'''
off = 0.25
for i, loc in enumerate(cloc):
    if i == 1:
        x, y = m(clon[i]+0.12, clat[i]+0.1)
        plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal')
    else:
        x, y = m(clon[i]-0.12, clat[i]+0.1)
        plt.text(x, y, loc, size=18, horizontalalignment='right', weight='normal')
'''        
'''
##########################################################################################
# add stations for 4.4
##########################################################################################

folder = 'Moe_4.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/tallen/Documents/Code/process_waves/stationlist.dat').readlines()
for stn in stns:
    for line in lines:
        if line.startswith(stn):
            dat = line.strip().split('\t')
            lat = float(dat[5])
            lon = float(dat[4])
            net = dat[6]
            
    stlat.append(lat)
    stlon.append(lon)
    stnet.append(net)
    print stn, net, lat
stlat = array(stlat)
stlon = array(stlon)
#unet = unique(array(stnet))
#print stns, stnet


# loop thru networks and plot
unet = ['AU', 'MEL', 'S', 'UOM']
sym = ['^', 'H', 'd', 's'] 
ms = [12, 13, 12, 12]
hnd = []
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = 'w', markeredgecolor='k', markeredgewidth=1, markersize=ms[i])
    hnd.append(h[0])

##########################################################################################
# add stations for 5.4
##########################################################################################

folder = 'Moe_5.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/tallen/Documents/Code/process_waves/stationlist.dat').readlines()
for stn in stns:
    for line in lines:
        if line.startswith(stn+'\t'):
            dat = line.strip().split('\t')
            lat = float(dat[5])
            lon = float(dat[4])
            net = dat[6]
    stlat.append(lat)
    stlon.append(lon)
    stnet.append(net)
stlat = array(stlat)
stlon = array(stlon)

# loop thru networks and plot
hnd = []
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = '0.1', markeredgecolor='w', markeredgewidth=1.5, markersize=ms[i])
    hnd.append(h[0])
    
plt.legend(hnd, list(unet), fontsize=14, loc=4, numpoints=1)

##########################################################################################
# annotate earthquake
##########################################################################################

eqlat = -38.252
eqlon = 146.234
blats = [-39.75]
blons = [146.5]

# line connecting epi to beachball
x, y = m([blons[0], eqlon], [blats[0], eqlat])
plt.plot(x, y, 'k-', lw=1.5)

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25)

# Add beachballs for two events
x, y = m(blons, blats)

# Two focal mechanisms for beachball routine, specified as [strike, dip, rake]
focmecs = [[122, 22, 167]]
for i in range(len(focmecs)):
    b = Beach(focmecs[i], xy=(x[i], y[i]), width=100000, linewidth=1, facecolor='k')
    b.set_zorder(10)
    ax.add_collection(b)

##########################################################################################
# draw box for fig 2
##########################################################################################

# draw box
f2lat = [-38.55, -37.75, -37.75, -38.55, -38.55]
f2lon = [145.7, 145.7, 146.75, 146.75, 145.7]
x, y = m(f2lon, f2lat)
plt.plot(x, y, 'k--', lw=1.5)

# write text
ft2lat = min(f2lat)+0.
ft2lon = max(f2lon)+0.1
x, y = m(ft2lon, ft2lat)
plt.text(x, y, 'Figure 2', size=14, horizontalalignment='left', weight='bold', style='italic')
'''
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.03, loc=4)

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

# normalise
norm = mpl.colors.Normalize(vmin=minslip, vmax=maxslip)

cax = fig.add_axes([0.84,0.3,0.02,0.4]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=1., norm=norm) # 
#cb.set_ticks(years)
#cb.set_ticklabels(labels)
cb.set_label('Long-Term Slip Rate (m/Ma)', rotation=270, labelpad=20, fontsize=15)


plt.savefig('seaus_dams_faults_map.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
