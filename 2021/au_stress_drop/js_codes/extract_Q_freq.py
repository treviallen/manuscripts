# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:02:23 2019

@author: u93322

Module looks up the Q value for a given frequecy band, provided the band
is in 'parameters.txt' attenuation files.  
"""

# script containing function to extract a Q value given an input lon/lat. 
# also contains the fucntion to produce a list of the frequecies to be used
# in the attenuation correction. 

import numpy as np

from obspy.geodetics import degrees2kilometers
from obspy.geodetics import gps2dist_azimuth
from math import cos, asin, sqrt

#import geographiclib.geodesic as geo


# this will be activated by a loop through the 50 frequency bands
# Freqs lower than 0.1 not included.  
def distance(lat1, lon1, lat2, lon2):
    """ distance on the globe using the Haversine formula
    """
    p = 0.017453292519943295
    a = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p)*cos(lat2*p) * (1-cos((lon2-lon1)*p)) / 2
    return 12742 * asin(sqrt(a))


def closest(data, v):
    return min(data, key=lambda p: distance(v['lat'],v['lon'],p['lat'],p['lon']))


def read_Q_layer_info(freq):
    """
    Reads the file parameters.txt and outputs the layer index 
    number given a queried frequency value.  
    """
    with open('/Users/trev/Documents/Manuscripts/manuscripts/2021/au_stress_drop/wei_etal_2017_lgq/parameters.txt') as f:
        next(f)
        lines = f.readlines()
        for line in lines:
            line = line.split()
            if np.float(line[1]) == freq:
                return line[0]
        
def get_Q_file(freq):
    """ returns a the file name for a given freq band
    """
    freq_index = read_Q_layer_info(freq)
    if len(freq_index) == 1:
        freq_index = "0" + str(freq_index)
    freq_file_string = "/Users/trev/Documents/Manuscripts/manuscripts/2021/au_stress_drop/wei_etal_2017_lgq/Q." + str(freq_index) +  "..txt"

    return freq_file_string

"""
def get_Q_value(freq, lon, lat):
    ''' returns Q value given a lon/lat pair at one point
    '''
    from numpy import loadtxt, where
    
    lon, lat = np.float(lon), np.float(lat)
    file_name = get_Q_file(freq)
    #file_name = '/Users/trev/Documents/Manuscripts/manuscripts/2021/au_stress_drop/wei_etal_2017_lgq/'+file_name
    #data_list = []
    '''
    with open(file_name) as f:
        lines = f.readlines()[1::]
        for line in lines:
            d = {}
            line = line.split()
            d['lon'] = np.float(line[0])
            d['lat'] = np.float(line[1])
            d['Q'] = np.float(line[2])
            data_list.append(d)
    '''
    data = loadtxt(file_name)
    
    idx = where((data[:,0] 
    

    v = {'lat': lat, 'lon': lon}
    return closest(data_list, v)
"""
def get_Q_data(freq):
    """ returns Q value given a lon/lat pair at one point
    """
    from numpy import loadtxt
    file_name = get_Q_file(freq)
    data = loadtxt(file_name)
    
    return closest(data_list, v)

def get_freq_list():
    """ function returns the list of frequecies to be used for attenuation
    correction in the stress drop study
    """
    freqs = []
    with open('/Users/trev/Documents/Manuscripts/manuscripts/2021/au_stress_drop/wei_etal_2017_lgq/parameters.txt') as f:
        lines = f.readlines()
        for line in lines[8:]:
            line = line.split()
            freqs.append(line[1])
    return freqs