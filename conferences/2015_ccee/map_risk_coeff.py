def rgb2tuple(tup):
    return (tup[0]/255.,tup[1]/255.,tup[2]/255.)

'''
Start main
'''
from mpl_toolkits.basemap import Basemap
from numpy import array, arange, mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from sys import argv
import shapefile
from mapping_tools import drawshapepoly

mpl.rcParams['pdf.fonttype'] = 42
figure = plt.figure(1,figsize=(11,10))

# Projection used for National Mapping
m = Basemap(width=5500000,height=5000000,
            rsphere=(6378137.00,6356752.3142), \
            resolution='l',area_thresh=1000.,projection='lcc', \
            lat_1=49.,lat_2=77,lat_0=63,lon_0=-90.)

# get inputs 
rtgmfile = argv[1]

# make output file
outpdf = rtgmfile.split('.')[0]

# get period
pertxt = outpdf.split('_')[-1]
pertxt = pertxt[0]+'.'+pertxt[1]

m.shadedrelief()            
m.drawcoastlines(linewidth=0.4)
m.drawparallels(arange(42.,91.,6.), linewidth=0.3, labels=[1,0,0,0], color='0.25')
m.drawmeridians(arange(-180.,181.,10.), linewidth=0.3, labels=[0,0,0,1], color='0.25')
m.drawcoastlines(linewidth=0.5,color='k')
#m.drawmapboundary(fill_color='lightgray')
#m.fillcontinents(color='white',lake_color='lightgray',zorder=0)

# load lake shpfiel
lakeshp = '//Users//tallen//Documents//DATA/GIS//ArcGIS//Lakes//AC_1M_Waterbodies_large_WGS84.shp'
sf = shapefile.Reader(lakeshp)
drawshapepoly(m, plt, sf, col='lightblue', fillshape=True, lw=0.15)

m.drawcountries(color='k')

# read risk ceof
lat = []
lon = []
riskcoeff = []
lines = open(rtgmfile).readlines()[1:]
for line in lines:
    dat = line.strip().split(',')
    lon.append(float(dat[1]))
    lat.append(float(dat[2]))
    riskcoeff.append(float(dat[4]))

print 'riskcoeff',min(riskcoeff),'-', max(riskcoeff)

x, y = m(array(lon), array(lat))
m.scatter(x, y, c=riskcoeff, s=40, cmap=plt.cm.Spectral_r, vmin=0.84, vmax=1.0, alpha=1.)
#plt.title('Risk Coefficient at Sa '+pertxt+' s', fontsize=14)

cax = figure.add_axes([0.915,0.25,0.02,0.5])
norm = mpl.colors.Normalize(vmin=0.84, vmax=1.0)
cb = mpl.colorbar.ColorbarBase(cax,cmap=plt.cm.Spectral_r,norm=norm)
#cb = plt.colorbar(cax)
cb.set_label('Risk Coefficient', rotation=270, fontsize=13, labelpad=18)

#cax = fig.add_axes([0.34,0.05,0.34,0.05]) # setup colorbar axes.
#
#cb = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')

plt.savefig(outpdf+'.jpg' ,format='jpg', dpi=300)

plt.show()
