# -*- coding: utf-8 -*-
"""
Created on Tue Jul 04 16:21:19 2017

@author: u56903
"""

def get_iris_event_data(bulk, folder, timestr):     
    #from obspy.core.utcdatetime import UTCDateTime
    from obspy.clients.fdsn.client import Client
    #from obspy.fdsn import Client
    from os import path
    
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
            st = fdsn_client.get_waveforms(network=b[0], station=b[1], location=b[2],
                                           channel=b[3], starttime=b[4], endtime=b[5],
                                           attach_response=True)
            
            st = st.merge(method=0, fill_value='interpolate')
            sta += st 
            
            fname = '.'.join((timestr, b[0], b[1], 'mseed'))
            fpath = path.join(folder, fname.replace(':','.'))
            
            print 'Writing file:', fpath                
            st.write(fpath, format="MSEED")
            #return st
        except:
            print 'No data for', b[0], b[1]
            
    return sta
    
    
    
def get_arclink_event_data(bulk, fname):     
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.arclink.client import Client
    import obspy.clients.arclink
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
            
    client = Client(user='trevor.allen@ga.gov.au')
    st = client.get_waveforms(bulk[0], bulk[1], bulk[2], bulk[3], bulk[4], bulk[5])
    st = st.merge(method=0, fill_value='interpolate')
        
    print 'Writing file:', fname                   
    st.write(fname, format="MSEED")
    
    return st

    
################################################################################
# start main code here
################################################################################
from obspy.core.utcdatetime import UTCDateTime
import datetime as dt

################################################################################
# load usgs event list
################################################################################

usgscsv = 'neic_banda_png_gt_6.csv'
usgscsv = '2018-02-26_MW7.5_png.csv'

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
from os import path

# read dataless seed volumes
#parser = Parser(path.join('dataless', 'AU.dataless'))
#parser = Parser(path.join('..','..','..','Networks','AU', 'AU.seed'))
folder = 'mseed_dump'

evdict = parse_usgs_events(usgscsv)
for ev in evdict:
    # get AU network    
    t1 = ev['starttime'] + 60
    t2 = t1 + 2100
    '''
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
            ("AU", "CTA", "*", "*", t1, t2)]
            
    st = get_iris_event_data(bulk, folder, ev['timestr'][:16])

    ###########################################################################
    
    # get S network
    bulk = [("S", "AUDHS", "*", "*", t1, t2),
            ("S", "AUNHS", "*", "*", t1, t2),
            ("S", "AUKAT", "*", "*", t1, t2)]  
    
    st = get_iris_event_data(bulk, folder, ev['timestr'][:16])
    '''
    ###########################################################################
    
    # get GE network
    t1 = ev['starttime'] - 60
    t2 = t1 + 1500
    bulk = [("GE", "SAUI", "*", "BH*", t1, t2),
            ("GE", "SOEI", "*", "BH*", t1, t2),
            ("GE", "GENI", "*", "BH*", t1, t2),
            ("GE", "FAKI", "*", "BH*", t1, t2)]  
    
    #print fname
    
    try:
        for b in bulk:
            fname = '.'.join((ev['timestr'][:16],'GE',b[1],'mseed')).replace(':','.')
            fpath = path.join('mseed_dump', fname)
            st = get_arclink_event_data(b, fpath)
    except:
        print 'No GEOFON data...'
    
    
       




















