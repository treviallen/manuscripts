# -*- coding: utf-8 -*-
"""
Created on Thur Nov 21 13:02:23 2019

@author: u93322

Take an input source and receiver lat/lon and plots the Q values
against every frequency band.  

"""

from js_codes import extract_Q_freq
from math import log
import numpy as np
#import geographiclib.geodesic as geo
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from mapping_tools import distance, reckon

def convert_list_to_log(input):
    """converts any list into log space
    """
    vals = []
    for val in input:
        val = float(val)
        val_log = np.log10(val)
        vals.append(val_log)
    return vals


def calculate_ray(elat, elon, slat, slon, evdp):
    """function calculates the path along a seismic ray
    and outputs a list of lonlat values to be used in calculations
    """
    #print(elat, elon, slat, slon, evdp)
    model = TauPyModel(model="ak135")
    arrivals = model.get_ray_paths_geo(evdp, elat, elon, slat, slon, \
                                       phase_list=('S'))
    arrival = arrivals[0]
    path_array = arrival.path
    latlons = []
    print(path_array)
    for point in path_array:
        latlon = [point[4], point[5]]
        latlons.append(latlon)

    latlons = np.array(latlons)

    return latlons


def get_all_Qs_ray(elat, elon, slat, slon, evdp):
    """ 
    OPTION 1: average Q along the ray.  
    function calculates the Q value for every frequency in
    the Q list and outputs a list of Q values. 
    It also averages the Q value along the ray so requires ray
    geometry as an input.  
    """
    freqs = extract_Q_freq.get_freq_list()
    Q_list = []
    for freq in freqs:
        latlons = calculate_ray(elat, elon, slat, slon, evdp)
        Q_cumulative = 0
        for lat, lon in latlons:
            Q_dict = extract_Q_freq.get_Q_value(float(freq), lon, lat)
            Q_cumulative += Q_dict['Q']
        average_Q = Q_cumulative / len(latlons)
        Q_list.append(average_Q)

    return Q_list


def get_centre_ray(elat, elon, slat, slon, evdp):
    """ locates the centre point of a ray.  An alternitive to this
    would be to locate the deepest point of the ray, but this is
    unessecery
    """
    '''
    ray = calculate_ray(elat, elon, slat, slon, evdp)
    centre_idx = int(len(ray)/2)
    '''
    
    rngkm, az, baz = distance(elat, elon, slat, slon)
    lonlat = reckon(elat, elon, rngkm/2., az)
    
    return lonlat


def random_points_within(poly, num_points):
    """ gets a user defined number of random points
    within any polygon
    """
    min_x, min_y, max_x, max_y = poly.bounds

    points = []

    while len(points) < num_points:
        random_point = Point([np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)])
        if (random_point.within(poly)):
            points.append(random_point)

    return points


def get_all_Qs_centre(elat, elon, slat, slon, evdp, half_bin_size=3.0, npoints=8):
    """
    OPTION 2:
    function draws a default 6 degree polygon around the point at the centre 
    of a ray and randomly samples points within the polygon.  
    Sized based on the spacing in the Q model.  Also returns a standard deviation
    Q value calculated in log space.  

    This option chosen over the ray as there is less variation in the
    values, and the code is much quicker to run (smaller nested loop).  
    """
    centre = get_centre_ray(elat, elon, slat, slon, evdp)
    lon, lat = centre[1], centre[0]
    #print(lon, lat)
    # create polygon around central point
    # arbitarily use 3 degrees... ? 
    minlon = lon - half_bin_size
    maxlon = lon + half_bin_size
    minlat = lat - half_bin_size
    maxlat = lat + half_bin_size
    
    # define polygon
    #poly = Polygon([left_lower, left_upper, right_lower, right_upper])
    #points = random_points_within(poly, npoints)
    #print(points)

    freqs = extract_Q_freq.get_freq_list()
    Q_list = []
    stdevQ_list = []
    #print('looping freqs')
    for freq in freqs:
        print(freq)
        qdat = extract_Q_freq.get_Q_data(freq)
        
        # get q vals within region
        idx = where((qdat[:,0] >= minlon) & (qdat[:,0] <= maxlon) \
        	          (qdat[:,1] >= minlat) & (qdat[:,1] <= maxlat))[0] 
        
        Q_vals = qdat[:,2][idx]
        	
        # get mean & stdev
        Q_list.append(np.mean(Q_vals))
        stdevQ_list.append(np.stdev(Q_vals))
        
        
        '''
        Q_cumulative = 0
        Q_each_point_list = []
        #print('looping points')
        for i,point in enumerate(points):
            point = points[i]
            lonp, latp = point.coords.xy[0][0], point.coords.xy[1][0]
            Q_dict = extract_Q_freq.get_Q_value(float(freq), lonp, latp)
            #print(Q_dict)
            Q_cumulative += Q_dict['Q']
            Q_each_point_list.append(Q_dict['Q'])
        average_Q = Q_cumulative / npoints
        # calculate stdev in log10 space
        diff_cumulative = 0
        for Q in Q_each_point_list:
            diff_sqrd = (np.log10(Q) - np.log10(average_Q))**2
            diff_cumulative += diff_sqrd
        #standard error calculation for sample
        stdev_Q = (np.sqrt(diff_cumulative / (npoints-1))) / np.sqrt(npoints)

        Q_list.append(average_Q)
        stdevQ_list.append(stdev_Q) # stdev in log10 space
        '''

    return Q_list, stdevQ_list