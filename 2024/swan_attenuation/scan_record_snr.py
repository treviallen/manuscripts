#import pickle

#wadat = pickle.load(open('wa_recs.pkl','rb'))
#bkndat = pickle.load(open('wa_eqSource.pkl','rb'))

from obspy import read, Trace, Stream, UTCDateTime
from obspy.taup import TauPyModel
from readwaves import return_data, readeqwave, readbkn
from mapping_tools import distance, km2deg
from io_catalogues import parse_ga_event_query
#from response import get_response_info
from misc_tools import listdir_extension
from data_fmt_tools import return_sta_data, fix_stream_channels
from datetime import datetime, timedelta
from os import path, getcwd, remove
from numpy import asarray, where, array, unique
import matplotlib.pyplot as plt
plt.ion()
import matplotlib as mpl
mpl.style.use('classic')


from data_fmt_tools import get_stn_dataless_seed, get_station_distance, \
                           get_sta_cwb_data, remove_low_sample_data

###############################################################################
# plotting function
###############################################################################

def do_picks(st, eqdt):
    fig = plt.figure(1, figsize=(16, 11))
    
    # get reference time for record
    reftime = st[0].stats.starttime-UTCDateTime(eqdt)
    
    i = 0
    channels = []
    # loop through traces
    #print(st)
    for j, tr in enumerate(st):
        if i < 3:
           if tr.stats.channel[1] != 'N' and len(st) > 3:
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
           
           elif tr.stats.channel[1] != 'N' and len(st) == 2:
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
               
           elif len(st) == 3:
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
           
           elif len(st) < 3 and tr.stats.channel[1] != 'N':
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
    
    # labels last
    plt.xlabel('Time Since Earthquake (s)')
    
    plt.suptitle(' '.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), tr.stats.station)))
    
    # now traces are plotted, pick phases
    picks = asarray(plt.ginput(3,timeout=-1))

    # get x indicies for fft
    x1 = int((picks[0,0]-reftime) * tr.stats.sampling_rate)
    x2 = int((picks[1,0]-reftime) * tr.stats.sampling_rate)
    x3 = int((picks[2,0]-reftime) * tr.stats.sampling_rate)
    
    plt.close(fig)
    
    # check channels length
    if len(channels) == 2:
        channels.append('')
    if len(channels) == 1:
        channels.append('')
        channels.append('')
    
    return picks, x1, x2, x3, channels
    
###############################################################################
# acc plotting function
###############################################################################

def do_acc_picks(st, eqdt):
    fig = plt.figure(1, figsize=(16, 11))
    
    # get reference time for record
    reftime = st[0].stats.starttime-UTCDateTime(eqdt)
    
    i = 0
    channels = []
    # loop through traces
    #print(st)
    for j, tr in enumerate(st):
        if i < 3:
           if tr.stats.channel[1] == 'N' and len(st) > 3:
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
           
           elif tr.stats.channel[1] == 'N' and len(st) == 2:
               i += 1
               ax = plt.subplot(3, 1, i)
               plt_trace(tr, plt, ax, reftime)
               channels.append(tr.stats.channel)
               
    
    # labels last
    plt.xlabel('Time Since Earthquake (s)')
    
    plt.suptitle(' '.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), tr.stats.station)))
    
    # now traces are plotted, pick phases
    picks = asarray(plt.ginput(3,timeout=-1))

    # get x indicies for fft
    x1 = int((picks[0,0]-reftime) * tr.stats.sampling_rate)
    x2 = int((picks[1,0]-reftime) * tr.stats.sampling_rate)
    x3 = int((picks[2,0]-reftime) * tr.stats.sampling_rate)
    
    plt.close(fig)
    
    # check channels length
    if len(channels) == 2:
        channels.append('')
    if len(channels) == 1:
        channels.append('')
        channels.append('')
    
    return picks, x1, x2, x3, channels


# plot individual traces
def plt_trace(tr, plt, ax, reftime):    
    
    # set time relative to earthquake origin
    times = reftime + tr.times()
    
    # fiter record
    tr_filt = tr.copy()
    tr_filt.filter('bandpass', freqmin=0.7, freqmax=tr.stats.sampling_rate*0.45, corners=2, zerophase=True)
    
    # get absolute arrival times
    #pArrival = UTCDateTime(eqdt) + pTravelTime
    #sArrival = UTCDateTime(eqdt) + sTravelTime
    
    plt.plot(times, tr_filt.data, 'b-', lw=0.5, label='Data')
    #ax = plt.gca()
    ylims = ax.get_ylim()
    
    # set x lims based on distance
    if rngkm < 10.:
        plt.xlim([pTravelTime-5, pTravelTime+12])
        #plt.ylim([min(tr_filt.data)/10, max(tr_filt.data)/10]) 
    elif rngkm < 20.:
        plt.xlim([pTravelTime-10, pTravelTime+50])
    elif rngkm >= 20. and rngkm < 100.:
        plt.xlim([pTravelTime-20, pTravelTime+50])
    elif rngkm >= 100. and rngkm < 200:
        plt.xlim([pTravelTime-20, pTravelTime+100])
    elif rngkm >= 200. and rngkm < 400:
        plt.xlim([pTravelTime-30, pTravelTime+175])
    elif rngkm >= 400. and rngkm < 700:
        plt.xlim([pTravelTime-60, pTravelTime+250])
    elif rngkm >= 700. and rngkm < 1000:
        plt.xlim([pTravelTime-60, pTravelTime+350])
    elif rngkm >= 1000. and rngkm < 1300:
        plt.xlim([pTravelTime-60, pTravelTime+300])
    else:
        plt.xlim([pTravelTime-60, pTravelTime+750])
        
    if tr.stats.station == 'PIG4' or tr.stats.station == 'CVQOZ':
        plt.xlim([pTravelTime-120, pTravelTime+180])
    
    # set y lims
    xlims = ax.get_xlim()
    idx = where((times > xlims[0]) & (times < xlims[1]))[0]
    
    if len(idx) > 0:
        maxy = max(abs(tr_filt.data[idx]))
        
        plt.ylim([-1*maxy, maxy])
        #print(maxy)
    
    # plt theoretical arrivals
    plt.plot([pTravelTime, pTravelTime], ylims, 'r--', label='P Phase')
    plt.plot([sTravelTime, sTravelTime], ylims, 'g--', label='S Phase')
    plt.plot([sTravelTime2, sTravelTime2], ylims, 'k--', label='S Phase * 2')
    #if rngkm > 12:
    #    plt.plot([rgTravelTime, rgTravelTime], ylims, 'k--', label='Rg Phase')
    plt.plot([lgTravelTime, lgTravelTime], ylims, 'm--', label='Lg Phase')
    
    
    plt.ylabel(' '.join((tr.stats.channel,'-',str('%0.0f' % rngkm),'km')), fontsize=15)
    plt.legend(fontsize=10)

###############################################################################
# set velocity mdel
###############################################################################
model = TauPyModel(model="iasp91")

################################################################################
# parse events
################################################################################
# id,isotime,latitude,longitude,depth,magnitude

evdict = []
lines = open('swan.dd.dec2023.ge2.5.csv').readlines()[1:]

for line in lines:
    dat = line.strip().split(',')
    tmpdict = {'id': dat[0], 'datetime': UTCDateTime(dat[1]),
               'lat': float(dat[2]), 'lon': float(dat[3]), 
               'dep': float(dat[4]), 'mag': float(dat[5])}
    evdict.append(tmpdict)

################################################################################
# parse agmd recs
################################################################################

agmdevs = []
lines = open('20231019_aecom_record_list.csv').readlines()[1:]

# get unique event
evdt = []
for line in lines:
    dat = line.strip().split(',')
    evdt.append(UTCDateTime(dat[0]))
    
evdt = unique(array(evdt))

#re-loop thru events
for ut in evdt:
    for line in lines:
        dat = line.strip().split(',')
        if UTCDateTime(dat[0]) == ut:
            tmpdict = {'id': str(dat[0]), 'datetime': UTCDateTime(dat[0]),
                       'lat': float(dat[2]), 'lon': float(dat[1]), 
                       'dep': float(dat[3]), 'mag': float(dat[5])}
    agmdevs.append(tmpdict)

##############################################################################
# parse eq list
##############################################################################
gaevdict = parse_ga_event_query('2015_present_earthquakes_export.csv')

##############################################################################
# read ga sta list
##############################################################################
'''
iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-noarray.txt')
network = 'AU'
iris_sta_list = parse_iris_stationlist('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/iu-gmap-stations-autrim.txt')
network = 'IU'
'''
###############################################################################
# parse records
###############################################################################
outtxt = ''
"""records = 'preferred_records_edit.csv'
mseedfiles = open(records).readlines(dat[9])[0:]"""

mseedfiles = listdir_extension('auspass_dump', 'mseed')
#mseedfiles = listdir_extension('mseed_test', 'mseed')

m = 1
for mseedfile in mseedfiles:
    
###############################################################################
# first check to see if the pick file exists
###############################################################################
    recfile = path.split(mseedfile)[-1][:-5]+'picks'
    pick_path = path.join('record_picks',recfile)
    
    # check if pick file exists
    if not path.isfile(pick_path):
        
        # fix stream channels
        #fix_stream_channels(path.join('iris_dump', mseedfile))
    
        # read mseed
        #st = read(path.join('mseed_dump', mseedfile))
        st = read(path.join('auspass_dump', mseedfile))
        
        # remove junk channels
        cannotMerge = False
        try:
            st = remove_low_sample_data(st)
        except:
            cannotMerge = True
            
        if len(st) == 0:
            cannotMerge = True
        st = st.merge()
        ###############################################################################
        # associate event and get fe
        ###############################################################################
    
        evFound = False
        if cannotMerge == False:
            # look in RP's file first
            for evnum, ev in enumerate(evdict): 
                #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
                if st[0].stats.starttime+timedelta(seconds=60) < UTCDateTime(ev['datetime']) \
                   and st[0].stats.starttime+timedelta(seconds=180) > UTCDateTime(ev['datetime']):
                    evFound = True
                    eqlo = ev['lon']
                    eqla = ev['lat']
                    eqmag = ev['mag']
                    eqdp = ev['dep']
                    eqdt = ev['datetime']
            
        if evFound == False and cannotMerge == False:
            # look in AGMD list
            deltaT = st[0].stats.endtime-st[0].stats.starttime
            
            for evnum, ev in enumerate(agmdevs): 
                #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
                if deltaT > 300.:
                    if st[0].stats.starttime+timedelta(seconds=60) < UTCDateTime(ev['datetime']) \
                       and st[0].stats.starttime+timedelta(seconds=180) > UTCDateTime(ev['datetime']):
                        evFound = True
                        eqlo = ev['lon']
                        eqla = ev['lat']
                        eqmag = ev['mag']
                        eqdp = ev['dep']
                        eqdt = ev['datetime']
                        
                else:
                    if UTCDateTime(ev['datetime']) > st[0].stats.starttime-timedelta(seconds=30)  \
                       and UTCDateTime(ev['datetime']) < st[0].stats.endtime:
                        evFound = True
                        eqlo = ev['lon']
                        eqla = ev['lat']
                        eqmag = ev['mag']
                        eqdp = ev['dep']
                        eqdt = ev['datetime']
        
        if cannotMerge == False:
            for evnum, ev in enumerate(gaevdict): 
                #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
                if st[0].stats.starttime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
                   and st[0].stats.starttime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
                    evFound = True
                    eqlo = ev['lon']
                    eqla = ev['lat']
                    eqmag = ev['mag']
                    eqdp = ev['dep']
                    eqdt = ev['datetime']
                    
        # if event still not found, look at GA list
        if evFound == False and cannotMerge == False:
            # look in AGMD list
            deltaT = st[0].stats.endtime-st[0].stats.starttime
            
            for evnum, ev in enumerate(agmdevs): 
                #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
                if deltaT > 300.:
                    if st[0].stats.starttime+timedelta(seconds=60) < UTCDateTime(ev['datetime']) \
                       and st[0].stats.starttime+timedelta(seconds=180) > UTCDateTime(ev['datetime']):
                        evFound = True
                        eqlo = ev['lon']
                        eqla = ev['lat']
                        eqmag = ev['mag']
                        eqdp = ev['dep']
                        eqdt = ev['datetime']
                        
                else:
                    if UTCDateTime(ev['datetime']) > st[0].stats.starttime-timedelta(seconds=30)  \
                       and UTCDateTime(ev['datetime']) < st[0].stats.endtime:
                        evFound = True
                        eqlo = ev['lon']
                        eqla = ev['lat']
                        eqmag = ev['mag']
                        eqdp = ev['dep']
                        eqdt = ev['datetime']
        
        if evFound == True and cannotMerge == False:
            # get station details
            print('Getting picks for', mseedfile)
            
            sta_data = return_sta_data(st[0].stats.station)
            rngkm, azim, baz = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])
            rngdeg = km2deg(rngkm)
            
            if rngkm < 500. and eqla < -28. and eqlo < 125. and sta_data['sta'].startswith('NWAO'):
            
                # get arrivals
                if eqdp < 0:
                    arrival_dep = 0.
                else:
                    arrival_dep = eqdp
                    
                arrivals = model.get_travel_times(source_depth_in_km=arrival_dep, distance_in_degree=rngdeg)
                
                # find P and S
                p = []
                s = []
                for a in arrivals:
                    if a.name.upper() == 'P':
                        p.append(a.time)
                    if a.name.upper() == 'S':
                        s.append(a.time)
                        
                pTravelTime = p[0]
                sTravelTime = s[0]
                
                # estimate Lg Arrival (from Goulet et al 2014 P25-26)
                lgTravelTime = sTravelTime + 8.71 * 0.026*rngkm
                
                sTravelTime2 = 4 * sTravelTime
                
                # estimate Rg
                rgTravelTime = rngkm/3.05
                
                ###############################################################################
                # plot
                ###############################################################################
                # plot seismos only
                picks, x1, x2, x3, channels = do_picks(st, eqdt)
                
                # check acc
                if rngkm < 800 and x1 > x3:
                    doAcc = False
                    for tr in st:
                        if tr.stats.channel[1] == 'N':
                           doAcc = True
                    if doAcc == True:
                        picks, x1, x2, x3, channels = do_acc_picks(st, eqdt)
                        
                # capture data proc file
                tr = st[0]
                outtxt = ','.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                                   str(eqlo), str(eqla), str(eqdp), str(eqmag), str('%0.2f' % rngkm), str('%0.1f' % azim), \
                                   str(st[0].stats.sampling_rate), channels[0], channels[1], channels[2], \
                                   str('%0.3f' % picks[0,0]), str('%0.3f' % picks[1,0]), str('%0.3f' % picks[2,0]), \
                                   str(x1), str(x2), str(x3), \
                                   path.join(getcwd(), 'mseed_dump', mseedfile)))
                
                # only write if x3 > x1
                if x3 > x1:
                    # save processing file
                    
                    outfile = path.join('record_picks',recfile)
                    f = open(outfile, 'w')
                    f.write(outtxt)
                    f.close()
                 
                # write junk file so don't reveiw again
                else:
                    outfile = path.join('record_picks',recfile)
                    f = open(outfile, 'w')
                    f.write('junk')
                    f.close()
                    
            else:
                print('    R > 900 km')    
            
        elif evFound == False:
           print('Cannot associate event for:', mseedfile)
           print(st[0])
           #remove(path.join('mseed_dump', mseedfile))
           
        m += 1
        
