#import pickle

#wadat = pickle.load(open('wa_recs.pkl','rb'))
#bkndat = pickle.load(open('wa_eqSource.pkl','rb'))

from obspy import read, Trace, Stream, UTCDateTime
from obspy.taup import TauPyModel
from readwaves import return_data, readeqwave, readbkn
from mapping_tools import distance, km2deg
from response import get_response_info
from misc_tools import listdir_extension
from data_fmt_tools import return_sta_data
from datetime import datetime, timedelta
from os import path, getcwd, remove
from numpy import asarray
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
    tr_filt.filter('bandpass', freqmin=0.2, freqmax=tr.stats.sampling_rate*0.45, corners=2, zerophase=True)
    
    # get absolute arrival times
    #pArrival = UTCDateTime(eqdt) + pTravelTime
    #sArrival = UTCDateTime(eqdt) + sTravelTime
    
    plt.plot(times, tr_filt.data, 'b-', lw=0.5, label='Data')
    #ax = plt.gca()
    ylims = ax.get_ylim()
    
    # set x lims based on distance
    if rngkm < 20.:
        plt.xlim([pTravelTime-10, pTravelTime+60])
    elif rngkm >= 20. and rngkm < 100.:
        plt.xlim([pTravelTime-20, pTravelTime+200])
    elif rngkm >= 100. and rngkm < 200:
        plt.xlim([pTravelTime-20, pTravelTime+300])
    elif rngkm >= 200. and rngkm < 400:
        plt.xlim([pTravelTime-30, pTravelTime+400])
    elif rngkm >= 400. and rngkm < 700:
        plt.xlim([pTravelTime-60, pTravelTime+600])
    elif rngkm >= 700. and rngkm < 1000:
        plt.xlim([pTravelTime-60, pTravelTime+900])
    elif rngkm >= 1000. and rngkm < 1300:
        plt.xlim([pTravelTime-60, pTravelTime+1200])
    else:
        plt.xlim([pTravelTime-60, pTravelTime+1500])
    
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

##########################################################################################
# parse eq epicentres
##########################################################################################

def parse_usgs_events(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
    lines = open(usgscsv).readlines()[1:]
    
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'datetime': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict

usgscsv = '20200821_merged_events.csv'
#usgscsv = '2019-02-26_event.csv'
# parse catalogue
evdict = parse_usgs_events(usgscsv)
###############################################################################
# parse records
###############################################################################
outtxt = ''
"""records = 'preferred_records_edit.csv'
mseedfiles = open(records).readlines(dat[9])[0:]"""

#mseedfiles = listdir_extension('mseed_dump', 'mseed')
mseedfiles = listdir_extension('mseed_2020', 'mseed')

m = 1
for mseedfile in mseedfiles:
    
    
###############################################################################
# first check to see if the pick file exists
###############################################################################
    recfile = path.split(mseedfile)[-1][:-5]+'picks'
    pick_path = path.join('record_picks',recfile)
    
    # check if pick file exists
    if not path.isfile(pick_path):
        
        # read mseed
        #st = read(path.join('mseed_dump', mseedfile))
        st = read(path.join('mseed_2020', mseedfile))
        
        # remove junk channels
        st = remove_low_sample_data(st)
        
        ###############################################################################
        # associate event and get distance
        ###############################################################################
    
        evFound = False
        for evnum, ev in enumerate(evdict): 

            if st[0].stats.starttime > UTCDateTime(ev['datetime']-timedelta(seconds=901)) \
               and st[0].stats.starttime < UTCDateTime(ev['datetime']+timedelta(seconds=120)):
                evFound = True
                eqlo = ev['lon']
                eqla = ev['lat']
                eqmag = ev['mag']
                eqdp = ev['dep']
                eqdt = ev['datetime']
        
        if evFound == True:
            # get station details
            print('Getting picks for', mseedfile)
            
            sta_data = return_sta_data(st[0].stats.station)
            rngkm, azim, baz = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])
            rngdeg = km2deg(rngkm)
            
            if rngkm < 1800.:
            
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
                print('    R > 1800 km')    
            
        elif evFound == False:
           print('Cannot associate event for:', mseedfile)
           #remove(path.join('mseed_dump', mseedfile))
           
        m += 1
        
