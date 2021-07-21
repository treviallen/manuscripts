#from scipy.interpolate import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension
from mapping_tools import drawshapepoly, get_field_data

mpl.style.use('classic')
#plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
# parse eq epicentres
##########################################################################################

def parse_usgs_events(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
    lines = open(usgscsv).readlines()[1:]
    
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'time': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict


##########################################################################################
# parse SA files
##########################################################################################
# reads sa data files and returns period (T) and acceleration (SA) vectors
print('Loading pkl file...')
recs = pickle.load(open("stdict_ampfact.pkl", "rb" ))
#recs = pickle.load(open("stdict.pkl", "rb" ))

##########################################################################################
#108/152/-44/-8
urcrnrlat = -4.0
llcrnrlat = -26.
urcrnrlon = 142.
llcrnrlon = 117
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=500.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.fillcontinents(color='w',lake_color='0.9')
#m.drawmapboundary(fill_color='0.9')
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
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
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'

cmap, zvals = cpt2colormap(cptfile, 30)
cmap = remove_last_cmap_colour(cmap)

# make shading
print('Making map...')
ls = LightSource(azdeg = 180, altdeg = 5)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)

##########################################################################################
# add epicentres
##########################################################################################
# get colormap
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_10.cpt'
ncols = 10
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

# load shape
import shapefile
shpfile = 'shapefiles/2021_nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
zone_code = get_field_data(sf, 'CODE', 'str')
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
drawshapepoly(m, plt, sf, edgecolor='k', lw=2)

# get net IDs
netid = []
stas = []

unets = ['AU', 'GE', 'S1', 'II', 'IU', 'OA', '2O']
netsym = ['^', 'H', 'd', 's', 'v', 'p', 'D']
msize = [10, 10, 10, 8, 10, 10, 7]

# loop thru zones
i = 0
for j, poly in enumerate(polygons):
    
    if zone_code[j] == 'BS':
        eqlo = []
        eqla = []
        emag = []
        edep = []
        stlo = []
        stla = []
        stas = []
        netid = []
                
        # plt paths
        for rec in recs: #[0:100])):
            #if rec['dep'] >= 50:
            if rec['rhyp'] < 2000:
                pt = Point(rec['eqlo'], rec['eqla'])
                if pt.within(poly):
                    #plt eq
                    x1, y1 = m(rec['eqlo'], rec['eqla'])
                    x2, y2 = m(rec['stlo'], rec['stla'])
                    
                    plt.plot([x1, x2], [y1, y2], '-', c='0.5', lw=0.3, alpha=1)
                    
                    eqlo.append(rec['eqlo'])
                    eqla.append(rec['eqla'])
                    stlo.append(rec['stlo'])
                    stla.append(rec['stla'])
                    emag.append(rec['mag'])
                    edep.append(rec['dep'])
                    stas.append(rec['sta'])
                    netid.append(rec['network'])
                
            
        '''
        # plt stns and events
        for i in range(0, len(eqlo)):
        
            x, y = m(eqlo[i], eqla[i])
            plt.plot(x, y, 'o', mfc=cs[j*2+1], mec='k', markeredgewidth=0.5, markersize=(-15 + emag[i]*4.), alpha=1., zorder=10000+i)
            
            x, y = m(stlo[i], stla[i])
            plt.plot(x, y, '^', markerfacecolor='0.1', markeredgecolor='w', markeredgewidth=1., \
                     markersize=9, alpha=1)
        '''
###############################################################################
# map eqs
###############################################################################
# set colours
netid = array(netid)
deprngs = array([0, 100, 200, 300, 400, 500, 600, 700])
ncols = len(deprngs)-1
cmap = plt.get_cmap('viridis_r', ncols)
cs = (cmap(arange(ncols)))

# get zorder for plotting
emag = array(emag)
sortidx = argsort(argsort(emag))
for i in range(0, len(emag)):
    #get colour idx
    
    colidx = int(round((ncols-1) * edep[i] / max(deprngs)))
    x, y = m(eqlo[i], eqla[i])
    zo = sortidx[i] + 20
    plt.plot(x, y, 'o', mfc=list(cs[colidx]), markeredgecolor='k', markeredgewidth=0.5, \
             markersize=(-15 + emag[i]*4.), zorder=zo, alpha=0.8)

###############################################################################
# add stations
###############################################################################
unets = ['AU', 'GE', 'S1', 'II', 'IU', 'OA', '2O']
netsym = ['^', 'H', 'd', 's', 'v', 'p', 'D']

stas = array(stas)
ustas = unique(stas)

for us in ustas:
    pltsta = True
    for i in range(0, len(netid)):
       if stas[i] == us and pltsta == True:
           if netid[i] == 'AU' or netid[i] == 'DPH' or netid[i] == 'DRS':
               x, y = m(stlo[i], stla[i])
               print(stlo[i], stla[i])
               sym = netsym[0]
               ms = msize[0]
           elif netid[i] == 'GE':
               x, y = m(stlo[i], stla[i])
               sym = netsym[1]
               ms = msize[1]
           elif netid[i] == 'S1' or netid[i] == 'S':
               x, y = m(stlo[i], stla[i])
               sym = netsym[2]
               ms = msize[2]
           elif netid[i] == 'II':
               x, y = m(stlo[i], stla[i])
               sym = netsym[3]
               ms = msize[3]
           elif netid[i] == 'IU':
               x, y = m(stlo[i], stla[i])
               sym = netsym[4] 
               ms = msize[4]
           elif netid[i] == 'OA':
               x, y = m(stlo[i], stla[i])
               sym = netsym[5] 
               ms = msize[5]
           elif netid[i] == '2O':
               x, y = m(stlo[i], stla[i])
               sym = netsym[6] 
               ms = msize[6]
              
           plt.plot(x, y, sym, markerfacecolor='0.1', markeredgecolor='w', markeredgewidth=0.5, \
                    markersize=ms, alpha=1, zorder=1000)
           pltsta = False


###############################################################################
# add legends
###############################################################################

# make legend
legmag = [6., 7., 8.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mec='k', markersize=(-15 + lm*4.), alpha=1., lw=0.5)
    legh.append(h[0])

labels = ['$\mathregular{M_W}$ 6.0', '$\mathregular{M_W}$ 7.0', '$\mathregular{M_W}$ 8.0'] 
l = plt.legend(legh, labels, loc=3, numpoints=1)
#l.set_zorder(99)

# add 2nd legend here
legh = []
for un, ns in zip(unets, netsym):
    x, y = m(0, 0)
    h = plt.plot(x, y, ns, markerfacecolor='0.1', markeredgecolor='w', markeredgewidth=0.5, \
                 markersize=10, alpha=1)
    legh.append(h[0])
plt.legend(legh, unets, loc=2, numpoints=1)

# re-add leg 1
plt.gca().add_artist(l)

##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.11, loc=1)

m2 = Basemap(projection='merc',\
            llcrnrlon=104,llcrnrlat=-45, \
            urcrnrlon=160,urcrnrlat=5,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
            
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()

# fill main area
xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'r-', lw=1.5)

##########################################################################################
# add colourbar
##########################################################################################

#plt.gcf().subplots_adjust(bottom=0.07)
norm = mpl.colors.Normalize(vmin=0, vmax=max(deprngs))#myb
cax = fig.add_axes([0.77,0.3,0.015,0.4]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=0.8, norm=norm)
cb.set_ticks(deprngs)
cb.set_label('Hypocentral Depth (km)', rotation=270, labelpad=20, fontsize=18)


# set cb labels
#ticks = arange(0.5,len(bounds))
'''
cnt_rng = ['1', '2-3', '4-5', '6-10', '11-20', '21-30', '31-50', '51-70', '71-100', '101-150', '151+']
cb.set_ticks(ticks)
cb.set_ticklabels(cnt_rng)
'''
##########################################################################################
# label states
##########################################################################################
'''
state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -29.5, -23.0, -32.5, -37.1, -42.]
slon = [122, 133.5, 135.0, 144.5, 146.5, 143.6, 147.0]
for i, st in enumerate(state):
    x, y = m(slon[i], slat[i])
    #plt.text(x, y, st, size=11, horizontalalignment='center', verticalalignment='center', weight='normal')
'''
plt.savefig('figures/ncc_gm_paths.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
