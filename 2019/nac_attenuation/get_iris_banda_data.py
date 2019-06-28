# -*- coding: utf-8 -*-
"""
Created on Tue Jul 04 16:21:19 2017

@author: u56903
"""

def get_iris_event_data(bulk, folder, timestr, dataless, event):     
    from obspy import UTCDateTime
    from obspy.clients.fdsn.client import Client
    #from obspy.fdsn import Client
    from os import path
    from numpy import nan
    from mapping_tools import distance
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
            
    fdsn_client = Client("IRIS")
    sta = []
    #st = client.get_waveforms_bulk(bulk)
    for b in bulk:
        try:
            fname = '.'.join((timestr, b[0], b[1], 'mseed'))
            fpath = path.join(folder, fname.replace(':','.'))
            
            staloc = nan
            #first, check it site is in distance and azimuthal range
            for channel in ['SHZ', 'EHZ', 'BHZ', 'HHZ', 'BNZ', 'HNZ']:
                seedid = '.'.join((b[0], b[1], '00', channel)) #'AU.DPH.00.BNZ'
                try:
                    staloc = dataless.get_coordinates(seedid,b[4])
                except:
                    a=1 # dummy call
            
            # now get distance and azimuth
            rngkm, az, baz = distance(event['lat'], event['lon'], staloc['latitude'], staloc['longitude'])
            print(rngkm, az, baz)
            
            getRecord = False
            if rngkm <= 2000. and az > 130. and az < 230.:
                getRecord = True
                               
            # second, check if file exists
            if not path.isfile(fpath) and getRecord == True:
                 print(fpath)
                 st = fdsn_client.get_waveforms(network=b[0], station=b[1], location=b[2],
                                                channel=b[3], starttime=b[4], endtime=b[5],
                                                attach_response=True)
                 
                 print(st[0].stats.location)
                 st = st.merge(method=0, fill_value='interpolate')
                 sta += st 
                 
                 print('Writing file:', fpath)
                 st.write(fpath, format="MSEED")
            else:
                print('File exists:', fpath)
            #return st
        except:
            print('No data for', b[0], b[1])
            
    return sta
           
def get_arclink_event_data(bulk, fname, dataless, event):     
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.arclink.client import Client
    #import obspy.clients.arclink
    from os import path
    from numpy import nan
    from mapping_tools import distance
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
    try:
        staloc = nan
        #first, check it site is in distance and azimuthal range
        for channel in ['SHZ', 'EHZ', 'BHZ', 'HHZ', 'BNZ', 'HNZ']:
            seedid = '.'.join((b[0], b[1], '00', channel)) #'AU.DPH.00.BNZ'
            try:
                staloc = dataless.get_coordinates(seedid,b[4])
            except:
                a=1 # dummy call
        
        # now get distance and azimuth
        print(staloc)
        rngkm, az, baz = distance(event['lat'], event['lon'], staloc['latitude'], staloc['longitude'])
        print(rngkm, az, baz)
        
        getRecord = False
        if rngkm <= 2000. and az > 130. and az < 230.:
            getRecord = True
            
        # check if file already exists
        if not path.isfile(fname) and getRecord == True:
            client = Client(user='trevor.allen@ga.gov.au')
            st = client.get_waveforms(bulk[0], bulk[1], bulk[2], bulk[3], bulk[4], bulk[5])
            st = st.merge(method=0, fill_value='interpolate')
            
        print('Writing file:', fname)
        st.write(fname, format="MSEED")
    except:
        print('No data for:', fname)
    
    return st

################################################################################
# start main code here
################################################################################
from obspy.core.utcdatetime import UTCDateTime
import datetime as dt

################################################################################
# load usgs event list
################################################################################

usgscsv = '20190625_merged_events.csv'

def parse_usgs_events(usgscsv):
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

################################################################################
# load station lists
################################################################################

#from obspy.io.xseed import Parser
#from obspy.arclink.client import Client
#from obspy.core import read
from os import path, getcwd
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

# read shapefile 
shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')

from obspy.io.xseed import Parser

# read dataless seed volumes
print('Reading dataless seed volumes...')
if getcwd().startswith('/nas'):
    au_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU.IRIS.dataless')
    s1_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/S1/S1.IRIS.dataless')
    ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE.IRIS.dataless')
    iu_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/IU.IRIS.dataless')
    
else:
    au_parser = Parser('/Users/trev/Documents/Earthquake_Data/AU.dataless')
    
folder = 'mseed_dump'

cnt = 0
evdict = parse_usgs_events(usgscsv)
for ev in evdict:
    # check if event in polygons
    pt = Point(ev['lon'], ev['lat'])
    for poly, zcode in zip(polygons, zone_code):
        # set Mmin
        if zcode == 'BS':
           mmin = 5.25
        else:
           mmin = 5.75
           
        if pt.within(poly):
            print('Getting event data:', ev['time'], ev['lon'], ev['lat'])
    
            # only use MW
            if ev['magType'].upper().startswith('MW') and ev['mag'] >= mmin:
                cnt += 1
                
                # get AU network    
                t1 = ev['starttime'] + 60
                t2 = t1 + 2100
                
                bulk = [("AU", "PSAD2", "*", "*", t1, t2),
                        ("AU", "FITZ", "*", "*", t1, t2),
                        ("AU", "KNA", "*", "*", t1, t2),
                        ("AU", "KNRA", "*", "*", t1, t2),
                        ("AU", "MTN", "*", "*", t1, t2),
                        ("AU", "DPH", "*", "*", t1, t2),
                        ("AU", "DRS", "*", "*", t1, t2),
                        ("AU", "KDU", "*", "*", t1, t2),
                        ("AU", "KAKA", "*", "*", t1, t2),
                        ("AU", "COEN", "*", "*", t1, t2),
                        ("AU", "QIS", "*", "*", t1, t2),
                        ("AU", "CN1H", "*", "*", t1, t2),
                        ("AU", "CN2S", "*", "*", t1, t2),
                        ("AU", "MBWA", "*", "*", t1, t2),
                        ("AU", "MTSU", "*", "*", t1, t2),
                        ("AU", "WB2", "*", "*", t1, t2),
                        ("AU", "WRKA", "*", "*", t1, t2),
                        ("AU", "AS31", "*", "*", t1, t2),
                        ("AU", "RABL", "*", "*", t1, t2),
                        ("AU", "NTNOZ", "*", "*", t1, t2),
                        ("AU", "BOXOZ", "*", "*", t1, t2),
                        ("AU", "LAROZ", "*", "*", t1, t2),
                        ("AU", "GIRL", "*", "*", t1, t2),
                        ("AU", "CVQOZ", "*", "*", t1, t2),
                        ("AU", "MUKOZ", "*", "*", t1, t2),
                        ("AU", "CTA", "*", "*", t1, t2)]
                        
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], au_parser, ev)
                
                ###########################################################################
                
                # get S1 network
                bulk = [("S1", "AUDHS", "*", "*", t1, t2),
                        ("S1", "AUNHS", "*", "*", t1, t2),
                        ("S1", "AUKAT", "*", "*", t1, t2)]  
                
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], s1_parser, ev)
                
                ###########################################################################
                
                # get IU network
                bulk = [("IU", "MBWA", "*", "*", t1, t2),
                        ("IU", "CTAO", "*", "*", t1, t2),
                        ("IU", "PMG", "*", "*", t1, t2)]  
                
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], iu_parser, ev)
                
                ###########################################################################
                '''
                # get IA network
                t1 = ev['starttime'] - 60
                t2 = t1 + 1500 
                
                bulk = [("IA", "MMPI", "*", "*", t1, t2),
                        ("IA", "WAMI", "*", "*", t1, t2),
                        ("IA", "SRPI", "*", "*", t1, t2),
                        ("IA", "NBPI", "*", "*", t1, t2),
                        ("IA", "KMPI", "*", "*", t1, t2),
                        ("IA", "TLE2", "*", "*", t1, t2),
                        ("IA", "TLE", "*", "*", t1, t2),
                        ("IA", "ALKI", "*", "*", t1, t2),
                        ("IA", "ATNI", "*", "*", t1, t2),
                        ("IA", "BATI", "*", "*", t1, t2),
                        ("IA", "LRTI", "*", "*", t1, t2),
                        ("IA", "EDFI", "*", "*", t1, t2),
                        ("IA", "BASI", "*", "*", t1, t2),
                        ("IA", "WSI", "*", "*", t1, t2),
                        ("IA", "LBFI", "*", "*", t1, t2),
                        ("IA", "WBSI", "*", "*", t1, t2),
                        ("IA", "BMNI", "*", "*", t1, t2),
                        ("IA", "DBNI", "*", "*", t1, t2),
                        ("IA", "TWSI", "*", "*", t1, t2),
                        ("IA", "BASI", "*", "*", t1, t2),
                        ("IA", "WSI", "*", "*", t1, t2),
                        ("IA", "LBFI", "*", "*", t1, t2),
                        ("IA", "WBSI", "*", "*", t1, t2),
                        ("IA", "MBPI", "*", "*", t1, t2)]  
                
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16])
                '''
                
                ###########################################################################
                
                # get GE network
                t1 = ev['starttime'] - 60
                t2 = t1 + 1500
                bulk = [("GE", "SAUI", "*", "BH*", t1, t2),
                        ("GE", "SOEI", "*", "BH*", t1, t2),
                        ("GE", "GENI", "*", "BH*", t1, t2),
                        ("GE", "PMG", "*", "BH*", t1, t2),
                        ("GE", "MMRI", "*", "BH*", t1, t2),
                        ("GE", "BNDI", "*", "BH*", t1, t2),
                        ("GE", "PLAI", "*", "BH*", t1, t2),
                        ("GE", "FAKI", "*", "BH*", t1, t2)]  
                
                #print fname
                #st = get_iris_event_data(bulk, folder, ev['timestr'][:16], ge_parser, ev)
                
                #try:
                for b in bulk:
                    fname = '.'.join((ev['timestr'][:16],'GE',b[1],'mseed')).replace(':','.')
                    fpath = path.join('mseed_dump', fname)
                    try:
                        st = get_arclink_event_data(b, fpath, ge_parser, ev)
                    except:
                        b = 1
                #except:
                #    print('No GEOFON data...')
                