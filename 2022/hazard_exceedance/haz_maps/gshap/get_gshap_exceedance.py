# -*- coding: utf-8 -*-
"""
Created on Mon Jul 03 15:48:17 2017

@author: u56903
"""
"""
def split_shape_parts(sf):
    from shapely.geometry import Polygon
    from numpy import array
    
    shapes = sf.shapes()
    polygons = []

    for i, shape in enumerate(shapes):
        parts = shape.parts
        print parts
        parts.append(len(shape.points)-1) # adding last point
        
        for i in range (0, len(parts)-1):
            '''
            while p <= parts[prt+1]:
                x.append(shape.points[p][0])
                y.append(shape.points[p][1])
                p += 1
            
                if x[0] != x[-1] or y[0] != y[-1]:
                    x.append(x[0])
                    y.append(y[0])
            '''        
            polygons.append(Polygon(shape.points[parts[i]:parts[i+1]]))
            print array(polygons[-1].bounds)
        #polygons.append(Polygon(shape.points[parts[-1]:]))
    return polygons
"""

from sys import argv
from shapely.geometry import Point, Polygon
import shapefile
from numpy import nan, where, delete, isnan, array, sort
from mapping_tools import split_shape_parts
from scipy.constants import g

grdlist = argv[1]
exout = argv[2]

# parse file
lines = open(grdlist).readlines()

lon = []
lat = []
val = []
for line in lines:
    dat = line.strip().split('\t')
    lon.append(float(dat[0]))
    lat.append(float(dat[1]))
    val.append(float(dat[2]) / g) # convert from m/s to g
    
# now load shapefile
shpfile = '/nas/gemd/ehp/georisk_earthquake/hazard/DATA/GIS/AUS_adm/AUS_adm_lores.shp'

sf = shapefile.Reader(shpfile)
#sf = sf.shapes()

polys = split_shape_parts(sf)

inidx = []
for poly in polys:
    print array(poly.bounds)
    for i in range (0, len(lon)):
        pt = Point(lon[i], lat[i])
        if pt.within(poly):
            inidx.append(i)

inidx = sort(array(inidx))

# remove nan values
val = array(val)[inidx]

# now get exceedances            
exceedance_vals = [0.02, 0.04, 0.06, 0.08, 0.1, 0.14, 0.2, 0.3, 0.4, 0.8]
outtxt = ''
for ev in exceedance_vals:
    idx = where(val > ev)[0]
    
    exrat = 1.0*len(idx) / len(val)
    print ev, exrat
    outtxt += ','.join((str(ev), str(exrat))) + '\n'

#exout = 'gshap_exceedance_floor.csv'  
f = open(exout, 'wb')
f.write(outtxt)
f.close()







