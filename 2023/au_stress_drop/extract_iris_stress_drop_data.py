# -*- coding: utf-8 -*-
"""
Created on Mon Nov 05 15:49:10 2018

@author: u56903
"""

from io_catalogues import parse_ga_event_query
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
from mapping_tools import distance
from data_fmt_tools import return_sta_data, parse_iris_stationlist, get_iris_data, get_auspass_data, get_swan_data
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
recent_csv = 'au_ge_4.4_earthquakes_export_recent.csv'
#recent_csv = 'au_ge_4.4_earthquakes_export_bboo.csv'	
gadat = parse_ga_event_query(recent_csv)

##############################################################################
# read ga sta list
##############################################################################
"""
if getcwd().startswith('/nas'):
    
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-noarray.txt')
    network = 'AU'
    '''
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/iu-gmap-stations-autrim.txt')
    network = 'IU'
    
    iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/II/ii-gmap-stations-autrim.txt')
    network = 'II'
    '''
else:
    '''
    
    #iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/gmap-stations-gswa.txt')
    network = 'AU'
    network = 'S1'
    network = 'IU'
    network = 'II'		
    network = 'G'
    network = '2O'
"""    
##############################################################################
# loop through networks
##############################################################################

networks = ['AU', 'S1', 'IU', 'II', 'G', '2O', 'M8', '3B', 'YW', 'WG', 'OZ']
#networks = ['OZ']
#networks = ['S1', 'IU', 'II', 'G', '2O', 'M8']

for network in networks:
    if network == 'AU':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/gmap-stations-noarray.txt')
        #iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/as17-gmap-stations.txt')
    elif network == 'S1':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/S1/s1-gmap-stations.txt')
    elif network == 'IU':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/IU/iu-gmap-stations-autrim.txt')
    elif network == 'II':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/II/ii-gmap-stations-autrim.txt')
    elif network == 'G':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/G/g-gmap-stations-autrim.txt')
    elif network == '2O':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/2o-gmap-stations.txt')
    elif network == 'M8':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AUSPASS/m8-gmap-stations.txt')
    elif network == '3B':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/3b-gmap-stations.txt') 
    elif network == 'YW':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AU/yw-gmap-stations.txt') 
    elif network == 'WG':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/GSWA/wg-gmap-stations.txt')
    elif network == '5G':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AUSPASS/5g-gmap-stations.txt')
    elif network == '4N':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AUSPASS/4n-gmap-stations.txt')
    elif network == '4N':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/AUSPASS/3o-gmap-stations.txt')
    elif network == 'OZ':
        iris_sta_list = parse_iris_stationlist('/Users/trev/Documents/Networks/SRC/oz-gmap-stations.txt')
       
##############################################################################
# loop through events
##############################################################################
    for ev in gadat: #[40:]:
        mindist = 0
        if network == 'S1' or network == 'M8' or network == '5G' or network == 'OZ':
            maxdist = 800
            #maxdist = 750
        else:
            if ev['mag'] >= 4.5:
                maxdist = 2200
            elif ev['mag'] >= 4.5:
                maxdist = 1000
            else:
                maxdist = 800
        #maxdist = 1500
        #maxdist = 200 # already got 200 - 2200 km
        
        dt = ev['datetime']
        print(dt)
        # allow 2 mins pre-event - subs realised "get_iris_data" already pads by 2 mins, so have 4 mins
        dt = UTCDateTime(ev['datetime']) - 120
        dateTuple = (dt.year, dt.month, dt.day, dt.hour, dt.minute)
        
        if dt >= UTCDateTime(2019, 7, 14, 22, 0):
        
            # set string for dummy file
            evdate = str(UTCDateTime(dt.year, dt.month, dt.day, dt.hour, dt.minute)-240)[0:16].replace(':','.')
            
            # loop thru stations
            #iris_sta_list = [{'CAAN']
            getauspass = True
            for isl in iris_sta_list:
                
                # check if station is open
                if isl['starttime'] <= dt and isl['stoptime'] >= dt: # and dt.year >= 2020:
                    #if isl['starttime'] <= dt and isl['stoptime'] < dt: # and dt.year >= 2020:
                    
                    # check if in distance range
                    repi = distance(ev['lat'], ev['lon'], isl['lat'], isl['lon'])[0]
                    #print(repi)
                    if repi >= mindist and repi <= maxdist: # and isl['sta'].startswith('KIM'):
                        
                        # make dummy file name and see if exists
                        mseedfile = path.join('iris_dump', \
                                   '.'.join((evdate, network, isl['sta'], 'mseed')))
                        
                        if not path.isfile(mseedfile):
                            if network == 'M8' and getauspass == True:
                                get_auspass_data(dateTuple, durn=1800, network='M8')
                                getauspass = False
                            elif network == 'WG' and getauspass == True:
                                get_auspass_data(dateTuple, durn=1800, network='WG')
                                getauspass = False
                            elif network == '5G' and getauspass == True:
                                #get_auspass_data(dateTuple, durn=1800, network='5G')
                                get_swan_data(dateTuple, durn=1800, network='5G')
                                getauspass = False
                            elif network == '4N' and getauspass == True:
                                #get_auspass_data(dateTuple, durn=1800, network='4N')
                                get_swan_data(dateTuple, durn=1800, network='4N')
                                getauspass = False
                            elif network == '3O' and getauspass == True:
                                get_auspass_data(dateTuple, durn=1800, network='3O')
                                getauspass = False
                            elif network == 'OZ' and getauspass == True:
                                get_swan_data(dateTuple, durn=1800, network='OZ')
                                getauspass = False 
                            else:
                                st = get_iris_data(dateTuple, isl['sta'], network, durn=1800)
                            
                        else:
                            print('Skipping: ' + mseedfile)
            