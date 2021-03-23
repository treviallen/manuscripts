from misc_tools import listdir_extension
from io_catalogues import parse_ga_event_query
from data_fmt_tools import return_sta_data, parse_iris_stationlist
from mapping_tools import distance
from obspy import read
from os import path
from datetime import datetime, timedelta 

##############################################################################
# parse eq list
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')
gadat = parse_ga_event_query('au_ge_4.4_earthquakes_export_edit.csv')

##############################################################################
# get mseed files
##############################################################################

folder = 'cwb_legacy'
mseedfiles = listdir_extension(folder, 'mseed')

iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-noarray.txt')
network = 'AU'


##############################################################################
# loop thru events
##############################################################################
mindist = 200
maxdist = 2200

# loop through mseed files
for mseedfile in mseedfiles:
    # get event date-time from filename
    evdt = datetime.strptime(mseedfile[3:-6], '%Y%m%d%H%M')
    
    # loop through events
    for gad in gadat:
        if evdt > gad['datetime'] - timedelta(minutes=2) and \
           evdt < gad['datetime'] + timedelta(minutes=2):
           
           # load mseed
           st = read(path.join('cwb_legacy', mseedfile))
           
           # loop thru traces & calculate distance
           for tr in st:
               for isl in iris_sta_list:
                   # check if in distance range
                   
                   if isl['sta'] == str(tr.stats.station):
                       repi = distance(gad['lat'], gad['lon'], isl['lat'], isl['lon'])[0]
                       
                       if repi >= mindist and repi <= maxdist:
                           print gad['datetime'], str(tr.stats.station)
                           print repi
                           
                           # make filename
                           filedt = datetime.strftime((gad['datetime']-timedelta(minutes=4)), '%Y-%m-%dT%H.%M')
                           
                           outmseed = path.join('cwb_legacy', '.'.join((filedt, str(tr.stats.network), \
                                                str(tr.stats.station), 'mseed')))
                                                
                           tr.write(outmseed, format='MSEED')
                           
                           
    