def rgb2tuple(tup):
    return (tup[0]/255.,tup[1]/255.,tup[2]/255.)

'''
Start main

to run:

run map_risk_coeff_westcoast.py NBCC2015Loc_RTGM_10.csv NBCC2015Loc_mean_hazcurves_10.csv

'''
from mpl_toolkits.basemap import Basemap
from numpy import array, arange, mean, percentile
import matplotlib.pyplot as plt
import matplotlib as mpl
from sys import argv

mpl.rcParams['pdf.fonttype'] = 42
figure = plt.figure(1,figsize=(11,10))

# set map bounds
llcrnrlat = 47.5
urcrnrlat = 54.5
llcrnrlon = -132.5
urcrnrlon = -120.5

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
            resolution ='i',area_thresh=1000.)
            
# get inputs 
rtgmfile = argv[1]
locnfile = argv[2]

# make output file
outpdf = rtgmfile.split('.')[0]

# get period
pertxt = outpdf.split('_')[-1]
pertxt = pertxt[0]+'.'+pertxt[1]

#m.bluemarble()            
m.drawcoastlines(linewidth=0.4)
m.drawparallels(arange(42.,91.,2.), linewidth=0.3, labels=[1,0,0,0], color='0.25')
m.drawmeridians(arange(-180.,181.,3.), linewidth=0.3, labels=[0,0,0,1], color='0.25')
m.drawcoastlines(linewidth=0.5,color='k')
m.drawmapboundary(fill_color='lightgray')
m.fillcontinents(color='white',lake_color='lightgray',zorder=0)
m.drawcountries(color='k')

'''
# load lake shpfiel
lakeshp = '//Users//tallen//Documents//DATA/GIS//ArcGIS//Lakes//AC_1M_Waterbodies_large_WGS84.shp'
sf = shapefile.Reader(lakeshp)
drawshapepoly(m, plt, sf, col='lightblue', fillshape=True, lw=0.15)
'''

# read risk ceof
lat = []
lon = []
riskcoeff = []
rtref = []
lines = open(rtgmfile).readlines()[1:]
for line in lines:
    dat = line.strip().split(',')
    rtref.append(int(dat[0]))
    lon.append(float(dat[1]))
    lat.append(float(dat[2]))
    riskcoeff.append(float(dat[4]))
    
# read location file
locref = []
locs = []

lines = open(locnfile).readlines()[4:]
for line in lines:
    dat = line.strip().split(',')
    locref.append(int(dat[-3]))
    locs.append(dat[-2].strip())
    
# list key locations to plot
keylocs = ['Victoria', 'Vancouver (city hall)', 'Sidney', 'Prince George', 'Tofino', 'Terrace', \
           'Kitimat Townsite', 'Queen Charlotte City', 'Whistler', 'Kelowna', 'Penticton', 'Bella Coola',
           'Masset', 'Ucluelet', 'Campbell River', 'Nanaimo', 'Prince Rupert', 'Port Hardy',
           'Williams Lake', 'Quesnel', 'Squamish', 'Lillooet', 'Hope', 'Burns Lake', 'Gold River', '100 Mile House']

x, y = m(array(lon), array(lat))
m.scatter(x, y, c=riskcoeff, s=120, cmap=plt.cm.Spectral_r, vmin=0.84, vmax=1.0, alpha=1.)
#plt.title('Risk Coefficient at SA '+pertxt+' s', fontsize=14)

# now annotate with site locations
for kl in keylocs:
    for i, l in enumerate(locs):
        if l == kl:
            r1 = locref[i]
            # now loop thru RT file
            for j, r2 in enumerate(rtref):
                if r1 == r2:
                    # now we've gotten this far, annotate!
                    x, y = m(array(lon[j]+.06), array(lat[j]+.06))
                    plt.annotate(kl , (0,0), xycoords='data', xytext=(x, y), \
                                 color='k', fontsize=13)
                    
cax = figure.add_axes([0.91,0.25,0.02,0.5])
norm = mpl.colors.Normalize(vmin=0.84, vmax=1.0)
cb = mpl.colorbar.ColorbarBase(cax,cmap=plt.cm.Spectral_r,norm=norm)
#cb = plt.colorbar(cax)
cb.set_label('Risk Coefficient', rotation=270, fontsize=13, labelpad=18)

#cax = fig.add_axes([0.34,0.05,0.34,0.05]) # setup colorbar axes.
#
#cb = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')

plt.savefig(outpdf+'_westcoast.pdf' ,format='pdf', dpi=300)

plt.show()
