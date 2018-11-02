# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 16:11:33 2018

@author: u56903
"""
from numpy import array
from interpolate_hazard_grid import get_nsha18_haz_curves
###############################################################################
# parse haz values
###############################################################################

pgahaz = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHM_18/reporting/NSHA18-grids_and_maps/data_pkg/3_locality_data/3.2_hazard_curves/nsha18_haz_curves_PGA.csv'

places = []
lats = []
lons = []
rp475 = []
rp2475 = []

lines = open(pgahaz).readlines()[2:]

for line in lines:
    dat = line.strip().split(',')
    places.append(dat[0])
    lons.append(float(dat[1]))
    lats.append(float(dat[2]))
    rp475.append(float(dat[5]))
    rp2475.append(float(dat[11]))
    
lons = array(lons)
lats = array(lats)
rp475 = array(rp475)
rp2475 = array(rp2475)

###############################################################################
# parse haz values
###############################################################################

for lon, lat, place in zip(lons, lats, places):
    hazCurveDict = get_nsha18_haz_curves()