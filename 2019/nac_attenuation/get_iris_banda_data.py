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
    from numpy import nan, isnan
    from mapping_tools import distance
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
            
    fdsn_client = Client("IRIS")
    #client = Client("IRIS")
    sta = []
    #st = client.get_waveforms_bulk(bulk)
    for b in bulk:
        try:
            fname = '.'.join((timestr, b[0], b[1], 'mseed'))
            fpath = path.join(folder, fname.replace(':','.'))
            
            staloc = nan
            #first, check it site is in distance and azimuthal range
            for channel in ['SHZ', 'EHZ', 'BHZ', 'HHZ', 'BNZ', 'HNZ']:
                if b[0] == 'WRAB':
                    locCode = '10'
                else:
                    locCode = '00'
                seedid = '.'.join((b[0], b[1], locCode, channel)) # e.g., 'AU.DPH.00.BNZ'
                try:
                    staloc = dataless.get_coordinates(seedid,b[4])
                except:
                    a=1 # dummy call
                seedid = '.'.join((b[0], b[1], '', channel)) # e.g., 'AU.DPH..BNZ'
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
            elif rngkm <= 2000. and az > 120. and az < 240. and b[1] == 'RABL':
                getRecord = True
            elif rngkm <= 2000. and az > 120. and az < 240. and b[1] == 'PMG':
                getRecord = True
                               
            # second, check if file exists
            #print(path.isfile(fpath), getRecord)
            if not path.isfile(fpath) and getRecord == True:
                 bulk2 = [(b[0], b[1], b[2], "*", b[4], b[5])] #,
                 print('B2',bulk2)
#                         ("AU", "AFI", "1?", "BHE",  b[4], b[5])]
                 client = Client("IRIS")
                 #st = client.get_waveforms_bulk(bulk2)
                 st = client.get_waveforms(b[0], b[1], b[2], "*", b[4], b[5])
                 
                 '''
                 st = fdsn_client.get_waveforms(network=b[0], station=b[1], location=b[2],
                                                channel=b[3], starttime=b[4], endtime=b[5],
                                                attach_response=True)
                 '''
                 #print(st[0].stats.location)
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
    try:
        from obspy.arclink.client import Client
    except:
        from obspy.clients.arclink.client import Client
        #from obspy.clients.fdsn.client import Client
    from os import path
    from numpy import nan, isnan
    from mapping_tools import distance
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
    try:
        #first, check it site is in distance and azimuthal range
        for channel in ['SHZ', 'EHZ', 'BHZ', 'HHZ', 'BNZ', 'HNZ']:
            seedid = '.'.join((b[0], b[1], '00', channel)) #'AU.DPH.00.BNZ'
            
            try:    
                staloc = dataless.get_coordinates(seedid,b[4])
            except:
                a=1 # dummy call
            # try another seed id fmt
            seedid = '.'.join((b[0], b[1], '', channel)) #'AU.DPH.00.BNZ'
            try:
                staloc = dataless.get_coordinates(seedid,b[4])
            except:
                a=1 # dummy call
        
        # now get distance and azimuth
        rngkm, az, baz = distance(event['lat'], event['lon'], staloc['latitude'], staloc['longitude'])
        #print(rngkm, az, baz)
        print('arclink',seedid)
        getRecord = False
        if rngkm <= 2000. and az > 110. and az < 250.:
            getRecord = True
        elif rngkm <= 50.:
            getRecord = True
        
        # check if file already exists
        if not path.isfile(fname) and getRecord == True:
            print('Getting:', fname)
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
#from data_fmt_tools import get_iris_data

################################################################################
# load usgs event list
################################################################################

usgscsv = '20190625_merged_events.csv'
usgscsv = '2019-20-26_event.csv'
usgscsv = '20200511_events.csv'
usgscsv = '20200824_events.csv'
usgscsv = '20201008_events.csv'

def parse_usgs_events(usgscsv):
    lines = open(usgscsv).readlines()[1:]
    print(lines[0])
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

# parse catalogue
evdict = parse_usgs_events(usgscsv)

#a=b # kill

from obspy.io.xseed import Parser

# read dataless seed volumes
print('Reading dataless seed volumes...')
if getcwd().startswith('/nas'):
    
    au_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU.IRIS.dataless')
    s1_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/S1/S1.IRIS.dataless')
    ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE1.IRIS.dataless')
    ge2_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE1.IRIS.dataless')
    iu_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/IU.IRIS.dataless')
    ii_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/II/II.IRIS.dataless')
    
else:
    
    au_parser = Parser('/Users/trev/Documents/Networks/AU/AU.IRIS.dataless')
    s1_parser = Parser('/Users/trev/Documents/Networks/S1/S1.IRIS.dataless')
    iu_parser = Parser('/Users/trev/Documents/Networks/IU/IU.IRIS.dataless')
    ge1_parser = Parser('/Users/trev/Documents/Networks/GE/GE1.IRIS.dataless')
    ge2_parser = Parser('/Users/trev/Documents/Networks/GE/GE2.IRIS.dataless')
    ii_parser = Parser('/Users/trev/Documents/Networks/II/II.IRIS.dataless')
    
folder = 'mseed_dump'

cnt = 0

for ev in evdict[0:]:
    #print(ev.keys())
    # check if event in polygons
    pt = Point(ev['lon'], ev['lat'])
    for poly, zcode in zip(polygons, zone_code):
        # set Mmin
        if zcode == 'BS' or zcode == 'PNGT':
           mmin = 5.25
        else:
           mmin = 5.75
        
        if pt.within(poly):
            print('Getting event data:', ev['time'], ev['lon'], ev['lat'])
    
            # only use MW
            if ev['magType'].upper().startswith('MW') and ev['mag'] >= mmin:
                cnt += 1
                
                """
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
                        ("AU", "TV1H", "*", "*", t1, t2),
                        ("AU", "TV2S", "*", "*", t1, t2),
                        ("AU", "BW1H", "*", "*", t1, t2),
                        ("AU", "BW2S", "*", "*", t1, t2),
                        ("AU", "RK1H", "*", "*", t1, t2),
                        ("AU", "RK2S", "*", "*", t1, t2),
                        ("AU", "GD1S", "*", "*", t1, t2),
                        ("AU", "EIDS", "*", "*", t1, t2),
                        ("AU", "MBWA", "*", "*", t1, t2),
                        ("AU", "MTSU", "*", "*", t1, t2),
                        ("AU", "WB2", "*", "*", t1, t2),
                        ("AU", "WRKA", "*", "*", t1, t2),
                        ("AU", "AS31", "*", "*", t1, t2),
                        ("AU", "RABL", "*", "*", t1, t2),
                        ("AU", "NTNOZ", "*", "*", t1, t2),
                        ("AU", "BOXOZ", "*", "*", t1, t2),
                        ("AU", "LAROZ", "*", "*", t1, t2),
                        ("AU", "AXCOZ", "*", "*", t1, t2),
                        ("AU", "GIRL", "*", "*", t1, t2),
                        ("AU", "MEEK", "*", "*", t1, t2),
                        ("AU", "CVQOZ", "*", "*", t1, t2),
                        ("AU", "MUKOZ", "*", "*", t1, t2),
                        ("AU", "CTA", "*", "*", t1, t2),
                        ("AU", "EIDS", "*", "*", t1, t2),
                        ("AU", "RMQ", "*", "*", t1, t2),
                        ("AU", "QLP", "*", "*", t1, t2),
                        ("AU", "INKA", "*", "*", t1, t2),
                        ("AU", "OOD", "*", "*", t1, t2)]
                        
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], au_parser, ev)
                
                ###########################################################################
                
                # get S1 network
                if ev['starttime'].year >= 2013:
                    bulk = [("S1", "AUDHS", "*", "*", t1, t2),
                            ("S1", "AUNHS", "*", "*", t1, t2),
                            ("S1", "AUAYR", "*", "*", t1, t2),
                            ("S1", "AUMOU", "*", "*", t1, t2),
                            ("S1", "AUCSH", "*", "*", t1, t2),
                            ("S1", "AUNRC", "*", "*", t1, t2),
                            ("S1", "AUKAR", "*", "*", t1, t2),
                            ("S1", "AUCAR", "*", "*", t1, t2),
                            ("S1", "AUKAT", "*", "*", t1, t2)]  
                    
                    st = get_iris_event_data(bulk, folder, ev['timestr'][:16], s1_parser, ev)
                """    
                ###########################################################################
                """
                # get IU network
                t1 = ev['starttime'] - 60
                t2 = t1 + 2100 
                bulk = [("IU", "MBWA", "*", "[BH][HN]*", t1, t2),
                        ("IU", "CTAO", "*", "[BH][HN]*", t1, t2),
                        ("IU", "PMG", "*", "[BH][HN]*", t1, t2)]  
                
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], iu_parser, ev)
                
                ###########################################################################
                
                # get II network
                t1 = ev['starttime']
                t2 = t1 + 2100 
                bulk = [("II", "WRAB", "10", "[BH][HN]*", t1, t2)]  
                '''
                for b in bulk:
                    fname = '.'.join((ev['timestr'][:16],'II',b[1],'mseed')).replace(':','.')
                    fpath = path.join('mseed_dump', fname)
                    try:
                        st = get_arclink_event_data(b, fpath, ii_parser, ev)
                    except:
                        b = 1
                '''
                timetupple = (t1.year,
                              t1.month,
                              t1.day,
                              t1.hour,
                              t1.minute)
                
                st = get_iris_event_data(bulk, folder, ev['timestr'][:16], ii_parser, ev)
                """
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
                
                #st = get_iris_event_data(bulk, folder, ev['timestr'][:16])
                
                for b in bulk:
                    fname = '.'.join((ev['timestr'][:16],'IA',b[1],'mseed')).replace(':','.')
                    fpath = path.join('mseed_dump', fname)
                    try:
                        st = get_arclink_event_data(b, fpath, ge_parser, ev)
                    except:
                        b = 1
                        
                ###########################################################################
                '''
                # get GE network
                t1 = ev['starttime'] - 60
                t2 = t1 + 1500
                bulk = [("GE", "SAUI", "*", "*", t1, t2),
                        ("GE", "SOEI", "*", "*", t1, t2),
                        ("GE", "GENI", "*", "*", t1, t2),
                        ("GE", "PMG",  "*", "*", t1, t2),
                        ("GE", "MMRI", "*", "*", t1, t2),
                        ("GE", "BNDI", "*", "*", t1, t2),
                        ("GE", "PLAI", "*", "*", t1, t2),
                        ("GE", "JAGI", "*", "*", t1, t2),
                        ("GE", "FAKI", "*", "*", t1, t2)]  
                
                #print fname
                #st = get_iris_event_data(bulk, folder, ev['timestr'][:16], ge_parser, ev)
                
                #try:
                
                for b in bulk:
                    fname = '.'.join((ev['timestr'][:16],'GE',b[1],'mseed')).replace(':','.')
                    fpath = path.join('mseed_dump', fname)
                    #st = get_arclink_event_data(b, fpath, ge_parser, ev)
                    try:
                        st = get_arclink_event_data(b, fpath, ge1_parser, ev)
                        #st = get_iris_event_data(b, folder, ev['timestr'][:16], ge1_parser, ev)
                    except:
                        try:
                            st = get_arclink_event_data(b, fpath, ge2_parser, ev)
                            #st = get_iris_event_data(b, folder, ev['timestr'][:16], ge2_parser, ev)
                        except:
                            b = 1
                
                #except:
                #    print('No GEOFON data...')
        else:
           print('Not in polygon:', ev['time'], ev['lon'], ev['lat'])