# -*- coding: utf-8 -*-
"""
Created on Fri Jan 08 15:29:15 2021

@author: u56903
"""

from obspy import read_events
from numpy import array, unique
from os import path

# exampleisc_1900-1960
'''
utc = events.events[-1]['picks'][0]['time']

lon = events.events[-20]['origins'][0]['longitude']

mag = events.events[-2]['magnitudes'][0]['mag']
mag_type = events.events[-2]['magnitudes'][0]['magnitude_type']

sta_amp = events.events[-2]['amplitudes'][0]['generic_amplitude'] # units of 'm'
sta_per = events.events[-2]['amplitudes'][0]['period']

repi = events.events[-2]['origins'][0]['arrivals'][0]['distance']
pick_id = events.events[-2]['origins'][0]['arrivals'][0]['pick_id']
'''

iscfiles = ['isc_1900-1960.xml', 'isc_1960-1970.xml', 'isc_1970-1975.xml', 'isc_1975-1980.xml', \
            'isc_1980-1985.xml', 'isc_1985-1990.xml']

ev_data = []
idxs = []
j = 0
c = 0 # initiate total event counter
# loop through files
for iscfile in iscfiles:
    print(iscfile)
    
    # read quakeML
    events = read_events(path.join('isc_xmls', iscfile), format="QUAKEML")
    
    # loop thru events in file
    for i, event in enumerate(events.events):
        c += 1
        
        # get basic event info
        ot = event['origins'][0]['time']
        lon = event['origins'][0]['longitude']
        lat = event['origins'][0]['latitude']
        
        # get magnitudes
        if 'magnitudes' in event.keys():
            idxs.append(i)
            mags = event['magnitudes']
            
            # loop through mags
            mags_tmp = []
            mtypes = []
            mauth = []
            for mag in mags:
                mags_tmp.append(mag['mag'])
                mtypes.append(mag['magnitude_type'])
                mauth.append(mag['creation_info']['author'])
                            
        # get amplitudes
        if 'amplitudes' in event.keys():
            amps = event['amplitudes']
            
            # loop through mags
            amps_tmp = []
            pers = []
            amptypes = []
            stas = []
            for amp in amps:
                if amp['type'].startswith('S') or amp['type'].startswith('s') or amp['type'].startswith('L'):
                    amps_tmp.append(amp['generic_amplitude'])
                    pers.append(amp['period'])
                    amptypes.append(amp['type'])
                    stas.append(amp['waveform_id']['station_code'])
                    
        if not amps_tmp == []:
            j += 1
            print(j, event['origins'][0]['time'])
                
            tmp_dict = {'ot': ot, 'lon':lon, 'lat':lat, 'mags':mags_tmp, 'mtypes':mtypes, 'mauth':mauth, \
                        'amps':amps_tmp, 'pers':pers, 'amptypes':amptypes, 'stas':stas}
            ev_data.append(tmp_dict)
    