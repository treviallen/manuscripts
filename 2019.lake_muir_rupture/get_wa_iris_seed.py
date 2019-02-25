# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 09:14:29 2018

@author: u56903
"""
import datetime as dt 
from mapping_tools import distance
from data_fmt_tools import get_iris_stn_dataless_seed, get_station_distance, \
                           get_iris_data

austationlist_la = ['RKGY', 'AS31', 'BBOO', 'CARL', 'GIRL', 'INKA', 'KLBR', 'KMBL', \
               'LCRK', 'MEEK', 'MORW', 'MTKN', 'MULG', 'FORT', 'NAPP', 'NWAO', \
               'SDAN', 'STKA', 'TWOA', 'HTT', 'WHYH', 'WB10', 'WRKA', 'BLDU', \
               'YAPP', 'LCRK', 'PSA00']
               
austationlist_sm = ['RKGY', 'CARL', 'GIRL', 'KLBR', 'KMBL', \
                    'MEEK', 'MORW', 'MTKN', 'BLDU']

sstationlist = ['AUALB', 'AUBUS', 'AUKUL', 'AUKAL', 'AUMAZ', 'AUHAR', 'AUROX', \
                'AUCAS', 'AUJCS']
               
iustationlist = ['NWAO', 'MBWA']             

#stationlist = ['BLDU']

'''
# for testing, get Moe data
eqdt = datetime.datetime(2012,06,19,10,53)
eqla = -38.304
eqlo = 146.200
sta = 'TOO'
'''
# get singel station data
td_start = -120
td_end = 1800.

ggcatcsv = 'lake_muir_events.csv'

lines =open(ggcatcsv).readlines()[1:]
    
#initiate event directory
evdict = []
    
for line in lines:
        dat = line.strip().split(',')
        try:
            dateTime = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%f')
        except:
            dateTime = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S')
            
        tdict = {'datetime': dateTime, \
             'year': dateTime.year, 'month': dateTime.month, \
             'day': dateTime.day, 'hour': dateTime.hour, \
             'minute': dateTime.minute, 'second': dateTime.second, \
             'lat': float(dat[2]), 'lon': float(dat[3]), \
             'dep': float(dat[4]), 'mag': float(dat[1])} 
    
        evdict.append(tdict)
        evdict

outtxt = ''

# get dataless seed volumes
print 'Parsing dataless seed volumes...'
au_dataless = get_iris_stn_dataless_seed('AU')
s_dataless = get_iris_stn_dataless_seed('S')
iu_dataless = get_iris_stn_dataless_seed('IU')

# now do loop
for evnum, ev in enumerate(evdict[0:1]): 
    
    eqdt = ev['datetime']
    eqlo = ev['lon']
    eqla = ev['lat']
    eqmag = ev['mag']
    eqdp = ev['dep']
    
    if eqdt.year == 2018 and eqdt.month >= 9:

        if eqmag >= 4.0:
            for sta in austationlist_la:
                # now get data stream & file name and dump to file
                st, msfile = get_iris_data((eqdt.year,eqdt.month,eqdt.day,eqdt.hour,eqdt.minute), sta, 'AU')
                
                try:
                    # get station data
                    stalocs, stations, channel, repi, azim = get_station_distance(st, au_dataless, eqlo, eqla)
                    outtxt += ','.join((eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                                        str(eqlo), str(eqla), str(eqdp), str(eqmag), str(repi[0]), str(azim[0]), \
                                        str(st[0].stats.sampling_rate), msfile)) + '\n'
                except:
                    print 'cannot locate station information', sta    
                    
        elif eqmag < 4.0:
            for sta in austationlist_sm:
                # now get data stream & file name and dump to file
                st, msfile = get_iris_data((eqdt.year,eqdt.month,eqdt.day,eqdt.hour,eqdt.minute), sta, 'AU')
                
                try:
                    # get station data
                    stalocs, stations, channel, repi, azim = get_station_distance(st, au_dataless, eqlo, eqla)
                    outtxt += ','.join((eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                                        str(eqlo), str(eqla), str(eqdp), str(eqmag), str(repi[0]), str(azim[0]), \
                                        str(st[0].stats.sampling_rate), msfile)) + '\n'
                except:
                    print 'cannot locate station information', sta
                     
        # now get seismometers in schools data
        for sta in sstationlist:
            # now get data stream & file name and dump to file
            st, msfile = get_iris_data((eqdt.year,eqdt.month,eqdt.day,eqdt.hour,eqdt.minute), sta, 'S')
            
            try:
                # get station data
                stalocs, stations, channel, repi, azim = get_station_distance(st, s_dataless, eqlo, eqla)
                outtxt += ','.join((eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                                    str(eqlo), str(eqla), str(eqdp), str(eqmag), str(repi[0]), str(azim[0]), \
                                    str(st[0].stats.sampling_rate), msfile)) + '\n'
            except:
                print 'cannot locate station information', sta
                
        # now get IU data
        for sta in iustationlist:
            # now get data stream & file name and dump to file
            st, msfile = get_iris_data((eqdt.year,eqdt.month,eqdt.day,eqdt.hour,eqdt.minute), sta, 'IU')
            
            try:
                # get station data
                stalocs, stations, channel, repi, azim = get_station_distance(st, iu_dataless, eqlo, eqla)
                outtxt += ','.join((eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                                    str(eqlo), str(eqla), str(eqdp), str(eqmag), str(repi[0]), str(azim[0]), \
                                    str(st[0].stats.sampling_rate), msfile)) + '\n'
            except:
                print 'cannot locate station information', sta
                
f = open('wa_records_iris.csv', 'wb')
f.write(outtxt)
f.close()