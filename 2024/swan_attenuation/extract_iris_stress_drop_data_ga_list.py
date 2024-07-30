# -*- coding: utf-8 -*-
"""
Created on Mon Nov 05 15:49:10 2018

@author: u56903
"""

from io_catalogues import parse_ga_event_query
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
from mapping_tools import distance
from data_fmt_tools import return_sta_data, parse_iris_stationlist, get_iris_data
import datetime as dt
from numpy import arange, array, where, zeros_like, histogram
#from gmt_tools import cpt2colormap 
#from misc_tools import remove_last_cmap_colour
from os import getcwd, path
from obspy import UTCDateTime

##############################################################################
# parse eq list
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')
#recent_csv = 'au_ge_4.4_earthquakes_export_edit.csv'
yilgarncsv = 'swsz_M3.2_4.0_2020-2024_earthquakes_export.csv'
#recent_csv = 'au_ge_4.4_earthquakes_export_bboo.csv'
gadat = parse_ga_event_query(yilgarncsv)

##############################################################################
# read ga sta list
##############################################################################

if getcwd().startswith('/nas'):
    
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-nwao.txt')
    network = 'AU'
    '''
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/iu-gmap-stations-autrim.txt')
    network = 'IU'
    
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/II/ii-gmap-stations-autrim.txt')
    network = 'II'
    '''
else:
    
    #iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/gmap-stations-noarray.txt')
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/au-gmap-stations-nwao.txt')
    network = 'AU'
    '''
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/S1/s1-gmap-stations.txt')
    network = 'S1'
    
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/IU/iu-gmap-stations-autrim.txt')
    network = 'IU'
    
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/II/ii-gmap-stations-autrim.txt')
    network = 'II'		
    
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/G/g-gmap-stations-autrim.txt')
    network = 'G'
    
    iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/2o-gmap-stations.txt')
    network = '2O'
    '''
##############################################################################
# loop through events
##############################################################################



# loop thru events
for ev in gadat: #[40:]:
    mindist = 0
    maxdist = 500
    '''
    if network == 'S1':
        maxdist = 1000
        #maxdist = 750
    else:
        if ev['mag'] >= 4.5:
            maxdist = 2200
        else:
            maxdist = 1500
    #maxdist = 1500
    #maxdist = 200 # already got 200 - 2200 km
    '''
    
    dt = ev['datetime']
    print(dt)
    # allow 2 mins pre-event - subs realised "get_iris_data" already pads by 2 mins, so have 4 mins
    dateTuple = (dt.year, dt.month, dt.day, dt.hour, dt.minute)
    
    # set string for dummy file
    evdate = str(UTCDateTime(dt.year, dt.month, dt.day, dt.hour, dt.minute))[0:16].replace(':','.')
    
    # loop thru stations
    #iris_sta_list = [{'CAAN']
    for isl in iris_sta_list:
        
        # check if station is open
        if isl['starttime'] <= dt and isl['stoptime'] >= dt and dt.year >= 2018 and dt.month > 9 and dt.day > 14:
            
            # check if in distance range
            repi = distance(ev['lat'], ev['lon'], isl['lat'], isl['lon'])[0]
            
            if repi >= mindist and repi <= maxdist: # and isl['sta'].startswith('KIM'):
                
                # make dummy file name and see if exists
                mseedfile = path.join('iris_dump', \
                           '.'.join((evdate, network, isl['sta'], 'mseed')))
                
                if not path.isfile(mseedfile):
                    st = get_iris_data(dateTuple, isl['sta'], network, durn=600)
                    
                else:
                    print('Skipping: ' + mseedfile)
