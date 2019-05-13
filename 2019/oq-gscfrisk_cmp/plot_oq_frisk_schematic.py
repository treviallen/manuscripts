from mpl_toolkits.basemap import Basemap
from numpy import array, arange, mean, percentile, array, unique, where, ones_like, zeros_like, sin, radians
import matplotlib.pyplot as plt
from mapping_tools import drawshapepoly, distance, reckon
import shapefile
from mapping_tools import get_field_data
from fault_tools import *
from shapely.geometry import Point, Polygon

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
import matplotlib.gridspec as gridspec

figure = plt.figure(1, figsize=(19, 10))

gs1 = gridspec.GridSpec(1, 2)
maprows=1
hspace = 0
gs1.update(wspace=0.05, hspace=hspace) # negative looks bad in "show", but ok in pngs

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']
##############################################################################
# loop thru files and parse hazard grids
########################################

cnrs = [-81, -57, 38.6, 56]
mapres = 'l'
grdsize = 4.
loc_lon = -76.49
loc_lat = 44.23
dist_cutoff = 800.

lagrd = []
logrd = []
for lo in arange(cnrs[0]-5, cnrs[1]+3):
    for la in arange (cnrs[2], cnrs[3]+3):
        logrd.append(lo)
        lagrd.append(la)
        
lagrd = array(lagrd)
logrd = array(logrd)

# raed shp
shpfile = 'random_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()

llcrnrlon = cnrs[0]
urcrnrlon = cnrs[1]
llcrnrlat = cnrs[2]
urcrnrlat = cnrs[3]

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

#plt.tick_params(labelsize=12)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution=mapres,area_thresh=600)

###############################################################
# OQ schematic
###############################################################
ax = figure.add_subplot(gs1[1])
#plt.subplot(122)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color='0.8', zorder=0)
m.fillcontinents(color='w',lake_color='0.8')
m.drawparallels(arange(-90.,90.,grdsize), labels=[0,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,grdsize*2), labels=[0,0,0,0], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
  
  
# drow zones
drawshapepoly(m, plt, sf, lw=2.)

# find grd points inside zones
shla = []
shlo = []
for shape in shapes:
    poly = Polygon(shape.points)
    
    for slo, sla in zip(logrd, lagrd):
        pt = Point(slo, sla)
        if pt.within(poly):
            shlo.append(slo)
            shla.append(sla)
    
shlo = array(shlo)
shla = array(shla)

# plt OQ grd points
x, y = m(shlo, shla)
m.plot(x, y, '+',c='0.5', ms=10, mew=1.5)

# get points LT dist cutoff
oq_dist = []
for la, lo in zip(shla, shlo):
    oq_dist.append(distance(la, lo, loc_lat, loc_lon)[0])
    
# plt pts within cutoff
oq_dist = array(oq_dist)
idx = where(oq_dist <= dist_cutoff)[0]
x, y = m(shlo[idx], shla[idx])
m.plot(x, y, 'r+', ms=10, mew=1.5)

# draw cut-off boundary
cutaz = arange(0,360, 0.5)
cutla = []
cutlo = []
for az in cutaz:
    cutlo.append(reckon(loc_lat, loc_lon, dist_cutoff, az)[0])
    cutla.append(reckon(loc_lat, loc_lon, dist_cutoff, az)[1])
cutlo.append(cutlo[0])
cutla.append(cutla[0])

# plt ring
x, y = m(cutlo, cutla)
m.plot(x, y, '-', c='seagreen', lw=3)

# plt site
x, y = m(loc_lon, loc_lat)
m.plot(x, y, 'r*', ms=30)

#plt.title('OpenQuake-engine', fontsize=22)
#plt.savefig('oq_schematic.png', fmt='png', bbox_inches='tight')

# plt letter
xlim = ax.get_xlim()
xtxt = xlim[1] * 0.98
ylim = ax.get_ylim()
ytxt = ylim[1] * 0.02
plt.text(xtxt, ytxt, pltlett[1], fontsize=30, va='bottom', ha='right')

###############################################################
# FRISK schematic
###############################################################
m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution=mapres,area_thresh=600)

#fig = plt.figure(2, figsize=(12,12))
#plt.subplot(121)
ax = figure.add_subplot(gs1[0])
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color='0.8', zorder=0)
m.fillcontinents(color='w',lake_color='0.8')
m.drawparallels(arange(-90.,90.,grdsize), labels=[0,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,grdsize*2), labels=[0,0,0,0], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
  
  
# drow zones
drawshapepoly(m, plt, sf, lw=2.)
cutaz = arange(0,360, 0.5)

# loop thru shapes and get distance
for shape in shapes:
    inCOD = False
    shla = []
    shlo = []
    d = []
    poly = Polygon(shape.points)
    
    for p in shape.points:
        # get distance
        d.append(distance(loc_lat, loc_lon, p[1], p[0])[0])
        shlo.append(p[0])
        shla.append(p[1])
        if min(d) < dist_cutoff:
            inCOD = True
            
    # get min/max d for slices
    mind = min(d)
    maxd = max(d)
    
    # get slice spacing
    slcint = (maxd-mind) / 10.
    slcdists = arange(mind-0*slcint, maxd+slcint, slcint)
    
    # loop thru slices
    for sd in slcdists:
        slcla = []
        slclo = []
    
        for az in cutaz:
            slo, sla = reckon(loc_lat, loc_lon, sd, az)
            # check if point in polygon
            pt = Point(slo, sla)
            if pt.within(poly):
                slclo.append(slo)
                slcla.append(sla)
            else: #try plotting and reset
                # plt ring
                if inCOD == True and len(slclo) > 0:
                    x, y = m(slclo, slcla)
                    m.plot(x, y, 'r--', lw=1.5)
                elif len(slclo) > 0:
                    x, y = m(slclo, slcla)
                    m.plot(x, y, 'k--', lw=1.5)
                slclo = []
                slcla = []
            
        # plt ring
        if inCOD == True:
            x, y = m(slclo, slcla)
            m.plot(x, y, 'r--', lw=1.5)
        else:
            x, y = m(slclo, slcla)
            m.plot(x, y, 'k--', lw=1.5)
    
    # now, make slices            
    if inCOD == True:
        x, y = m(shlo, shla)
        plt.fill(x, y, 'r', alpha=0.3)
        m.plot(x, y, 'r-', lw=2)
        

# draw cut-off boundary
cutla = []
cutlo = []
for az in cutaz:
    cutlo.append(reckon(loc_lat, loc_lon, dist_cutoff, az)[0])
    cutla.append(reckon(loc_lat, loc_lon, dist_cutoff, az)[1])
cutlo.append(cutlo[0])
cutla.append(cutla[0])

# plt ring
x, y = m(cutlo, cutla)
m.plot(x, y, '-', c='seagreen', lw=3)

# plt site
x, y = m(loc_lon, loc_lat)
m.plot(x, y, 'r*', ms=30)

#plt.title('GSCFRISK', fontsize=22)

# plt letter
'''
xlim = ax.get_xlim()
xtxt = xlim[1] * 0.02
ylim = ax.get_ylim()
ytxt = ylim[1] * 0.96
'''
plt.text(xtxt, ytxt, pltlett[0], fontsize=30, va='bottom', ha='right')


plt.savefig('oq_frisk_cartoon.png', fmt='png', bbox_inches='tight')
plt.show()


#points = sf.shapes()[i].points

'''
la = []
lo = []
for p in points:
    lo.append(p[0])
    la.append(p[1])

x, y = m(lo, la)

m.plot(x, y, 'r-', lw=2)

# get fault stats
ftype = 'ss'
area = flen[i] * seis_thickness / sin(radians(dip))
mchar = return_area2mag_for_ss(area)[1] - 0.25

# add title
plt.title('; '.join((fcode[i], 'L = '+str('%0.0f' % flen[i])+' km', 'Mc = '+str('%0.2f' % mchar))), fontsize=14)
'''

