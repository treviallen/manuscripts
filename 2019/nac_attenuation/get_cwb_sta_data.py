# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 17:04:43 2017

@author: u56903
"""

from obspy import read
from data_fmt_tools import get_sta_cwb_data
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

usgscsv = '20190625_merged_events.csv'

def parse_usgs_events(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
    lines = open(usgscsv).readlines()[1:]
    
    #2017-05-29T14:35:21.510Z
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'time': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict

def check_cwb_data(sta, ev):
    from data_fmt_tools import return_sta_data
    from mapping_tools import distance
    
    sta_data = return_sta_data(sta)
    
    rngkm, az, baz = distance(ev['lat'], ev['lon'], sta_data['stla'], sta_data['stlo'])
    
    if rngkm <= 2000. and az > 130. and az < 230.:
        get_sta_cwb_data(ev['time'].year, ev['time'].month, ev['time'].day, ev['time'].hour, \
                         ev['time'].minute, 0, 2100, sta)

# start main code here - loop thru events and get data
evdict = parse_usgs_events(usgscsv)

# read shapefile 
shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')
       
for ev in evdict:
    if ev['time'].year <= 2006:
        pt = Point(ev['lon'], ev['lat'])
        for poly, zcode in zip(polygons, zone_code):
            # set Mmin
            if zcode == 'BS' or zcode == 'PNGT':
               mmin = 5.25
            else:
               mmin = 5.75
           
            if pt.within(poly) and ev['magType'].upper().startswith('MW') and ev['mag'] >= mmin:
                #stas = ['DRS', 'DPH', 'KAKA', 'MTN', 'KNA']
                stas = ['BATI']
                
                for sta in stas:
                    check_cwb_data(sta, ev)
        