# -*- coding: utf-8 -*-
"""
Created on Fri Nov 02 16:11:24 2018

@author: u56903
"""
from numpy import array, zeros_like, where, delete
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import split_shape_parts

# parse maps file
pgagrid = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_maps_PGA/hazard_map-mean_1.csv'

lon = []
lat = []
pga = []

lines = open(pgagrid).readlines()[2:]

for line in lines:
    dat = line.strip().split(',')
    lon.append(float(dat[0]))
    lat.append(float(dat[1]))
    pga.append(float(dat[-1]))
    
lon = array(lon)
lat = array(lat)
pga = array(pga)
idx = zeros_like(pga)

# parse shapefile
inshape = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GIS/AUS_adm/AUS_adm_lores.shp'
    
sf = shapefile.Reader(inshape)

polygons = split_shape_parts(sf)

for i, poly in enumerate(polygons):
    print 'Poly', i
    for j in range(0,len(pga)):
        point = Point(lon[j], lat[j])
        if point.within(poly):
            idx[j] = 1.
            
delidx = where(idx == 0.)[0]

# now delete recs
lon = delete(lon, delidx)
lat = delete(lat, delidx)
pga = delete(pga, delidx)

'''
import matplotlib.pyplot as plt
plt.plot(lon, lat, '.')
'''
# check exceedances
exrng = array([0.02, 0.04, 0.06, 0.1, 0.14, 0.2, 0.3, 0.4, 0.8])

exstr = 'PGA_LEVELS,%_EXCEEDANCE\n'
for er in exrng:
    eidx = where(pga >= er)[0]
    
    percExeedance = 100. * (1. * len(eidx)) / len(pga)
    
    print er, percExeedance
    
    exstr += ','.join((str(er), str(str(percExeedance)))) + '\n'
    
f = open('nsha_2475_exceedance.csv', 'wb')
f.write(exstr)
f.close()    
    
    
    
#plt.show()
    
    
    
    
    
    
