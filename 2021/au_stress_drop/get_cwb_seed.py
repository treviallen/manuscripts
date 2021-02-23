# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 09:14:29 2018

@author: u56903
"""
import datetime as dt 
from numpy import floor
from misc_tools import checkfloat
from mapping_tools import distance
from io_catalogues import parse_ga_event_query
from data_fmt_tools import get_stn_dataless_seed, get_station_distance, \
                           get_sta_cwb_data, parse_iris_stationlist
##############################################################################
# parse eq list
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')
gadat = parse_ga_event_query('au_ge_4.4_earthquakes_export_edit.csv')

##############################################################################
# read ga sta list
##############################################################################

iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-noarray.txt')
network = 'AU'
#iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/iu-gmap-stations-autrim.txt')
#network = 'IU'

##############################################################################
# read ga sta list
##############################################################################

# get singel station data - use 4 mins padding
td_start = -240
td_end = 1560

##############################################################################
# loop through events
##############################################################################

mindist = 200
maxdist = 2200

# loop thru events
for ev in gadat[40:]:
    dt = ev['datetime']
    print(dt)
    if dt.year < 2007:
    
        for isl in iris_sta_list:
            
            # check if station is open
            if isl['stoptime'] >= dt:
                
                # check if in distance range
                repi = distance(ev['lat'], ev['lon'], isl['lat'], isl['lon'])[0]
                
                if repi >= mindist and repi <= maxdist:
                    st, msfile = get_sta_cwb_data(dt.year,dt.month,dt.day,dt.hour,dt.minute,td_start, td_end, isl['sta'])
        
