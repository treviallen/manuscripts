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
from numpy import asarray, ceil, log10, array
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
    
    
    
    i = 0
    channels = []
    # loop through traces
    #print(st)
    for j, tr in enumerate(st):
    	  # get reference time for record
        reftime = tr.stats.starttime-UTCDateTime(eqdt)
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
    
    # set x lims based on distance
    if rngkm < 10.:
        #plt.xlim([pTravelTime-10, pTravelTime+30])
        x1, x2 = [pTravelTime-10, pTravelTime+30]
    elif rngkm < 20.:
        #plt.xlim([pTravelTime-10, pTravelTime+60])
        x1, x2 = [pTravelTime-10, pTravelTime+60]
    elif rngkm >= 20. and rngkm < 100.:
        #plt.xlim([pTravelTime-40, pTravelTime+150])
        x1, x2 = [pTravelTime-40, pTravelTime+150]
    elif rngkm >= 100. and rngkm < 200:
        #plt.xlim([pTravelTime-20, pTravelTime+250])
        x1, x2 = [pTravelTime-20, pTravelTime+250]
    elif rngkm >= 200. and rngkm < 400:
        #plt.xlim([pTravelTime-30, pTravelTime+400])
        x1, x2 = [pTravelTime-30, pTravelTime+400]
    elif rngkm >= 400. and rngkm < 700:
        #plt.xlim([pTravelTime-60, pTravelTime+600])
        x1, x2 = [pTravelTime-60, pTravelTime+600]
    elif rngkm >= 700. and rngkm < 1000:
        #plt.xlim([pTravelTime-60, pTravelTime+900])
        x1, x2 = [pTravelTime-60, pTravelTime+900]
    elif rngkm >= 1000. and rngkm < 1300:
        #plt.xlim([pTravelTime-60, pTravelTime+1200])
        x1, x2 = [pTravelTime-60, pTravelTime+1200]
    else:
        #plt.xlim([pTravelTime-60, pTravelTime+1500])
        x1, x2 = [pTravelTime-60, pTravelTime+1500]
        
    if tr.stats.station == 'PIG4' or tr.stats.station == 'CVQOZ':
        #plt.xlim([pTravelTime-120, pTravelTime+180])
        x1, x2 = [pTravelTime-120, pTravelTime+180]
    
    plt.xlim([x1, x2])
    ylims = array(ax.get_ylim())
    
    # now get indexes
    xi1 = int(round((x1 * tr.stats.sampling_rate)))
    xi2 = int(round((x2 * tr.stats.sampling_rate)))
    '''
    try:
        maxy = max(abs(tr_filt.data[xi1:xi2]))
        	
        loy = abs(maxy/ylims[0])
        if loy < 0.2:
            ylims[0] = -0.8*10**(ceil(log10(maxy)))
            
        hiy = abs(maxy/ylims[1])
        if hiy < 0.2:
            ylims[1] = 0.8*10**(ceil(log10(maxy)))
        
        #print(ylims, loy, hiy)  
    except:
        print('Bad Channel: '+tr.stats.channel)  
    '''
    plt.ylim(ylims)
        
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

##############################################################################
# parse eq list
##############################################################################
#gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')
evdict = parse_ga_event_query('au_ge_4.4_earthquakes_export_edit.csv')

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

mseedfiles = listdir_extension('iris_dump', 'mseed')

#mseedfiles = listdir_extension('iris_dump', 'ms')

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
        st = read(path.join('iris_dump', mseedfile))
        
        if st[0].stats.station == 'WYMK4':
            for i in range(0,len(st)):
                st[i].stats.station = 'WYKM4'
                st[i].stats.network = '3B'
                print('WYKM4')
            st.write(path.join('iris_dump', mseedfile), format="MSEED")
            st = read(path.join('iris_dump', mseedfile)) 
            
        if st[0].stats.station == 'WYMK6':
            for i in range(0,len(st)):
                st[i].stats.station = 'WYKM6'
                st[i].stats.network = '3B'
                print('WYKM6')
            st.write(path.join('iris_dump', mseedfile), format="MSEED")
            st = read(path.join('iris_dump', mseedfile))    
        
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
            for evnum, ev in enumerate(evdict): 
                """
                if ev['datetime'] > UTCDateTime(2019,7,14,5,38) and ev['datetime'] < UTCDateTime(2019,7,14,5,40):
                    print('Broome')
                    toff = 2500
                    '''
                    elif ev['datetime'] > UTCDateTime(2019,7,14,5,54) and ev['datetime'] < UTCDateTime(2019,7,14,5,56):
                        print('Broome AS')
                        toff = 3360
                    '''
                else:
                """
                toff = 601
            
            for evnum, ev in enumerate(evdict):         
                if st[0].stats.starttime > UTCDateTime(ev['datetime']-timedelta(seconds=toff)) \
                   and st[0].stats.starttime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
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
            print(sta_data)
            print(eqla, eqlo)
            rngkm, azim, baz = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])
            rngdeg = km2deg(rngkm)
            
            if rngkm < 2250.:
            
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
                if st[0].stats.network == 'M8':
                    print(rngkm)
                print('    R > 2500 km')    
            
        elif evFound == False:
           print('Cannot fassociate event for:', mseedfile)
           #print(st[0])
           #remove(path.join('mseed_dump', mseedfile))
           
        m += 1
        
