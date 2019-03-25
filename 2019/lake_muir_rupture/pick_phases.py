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
from os import path, getcwd
from numpy import asarray, nan, isnan
import matplotlib.pyplot as plt
plt.ion()

from data_fmt_tools import get_stn_dataless_seed, get_station_distance, \
                           get_sta_cwb_data, remove_low_sample_data

###############################################################################
# plotting function
###############################################################################

def do_picks(tr, pTravelTime, sTravelTime, lgTravelTime, rgTravelTime):
    fig = plt.figure(1, figsize=(15, 4))
    
    # get reference time for record
    print UTCDateTime(eqdt)
    reftime = tr.stats.starttime-UTCDateTime(eqdt)
    #reftime=0
    
    # set time relative to earthquake origin
    times = reftime + tr.times()
    
    # fiter record
    tr_filt = tr.copy()
    tr_filt = tr_filt.taper(0.001, type='hann', max_length=None, side='both')
    tr_filt.filter('bandpass', freqmin=0.5, freqmax=10, corners=2, zerophase=True)
    
    # get absolute arrival times
    #pArrival = UTCDateTime(eqdt) + pTravelTime
    #sArrival = UTCDateTime(eqdt) + sTravelTime
    
    plt.plot(times, tr_filt.data, 'b-', lw=0.5, label='Data')
    ax = plt.gca()
    ylims = ax.get_ylim()
    
    # set x lims based on distance
    if rngkm < 20.:
        plt.xlim([pTravelTime-3, pTravelTime+12])
    elif rngkm >= 20. and rngkm < 100.:
        plt.xlim([pTravelTime-7, pTravelTime+40])
    elif rngkm >= 100. and rngkm < 200:
        plt.xlim([pTravelTime-10, pTravelTime+120])
    elif rngkm >= 200. and rngkm < 400:
        plt.xlim([pTravelTime-10, pTravelTime+180])
    elif rngkm >= 400. and rngkm < 700:
        plt.xlim([pTravelTime-20, pTravelTime+300])
    elif rngkm >= 700. and rngkm < 1200:
        plt.xlim([pTravelTime-20, pTravelTime+500])
    elif rngkm >= 1200. and rngkm < 1500:
        plt.xlim([pTravelTime-30, pTravelTime+900])
    else:
        plt.xlim([pTravelTime-40, pTravelTime+1200])
    
    # plt theoretical arrivals
    plt.plot([pTravelTime, pTravelTime], ylims, 'r--', label='P Phase')
    plt.plot([sTravelTime, sTravelTime], ylims, 'g--', label='S Phase')
    if rngkm > 12:
        plt.plot([rgTravelTime, rgTravelTime], ylims, 'k--', label='Rg Phase')
    plt.plot([lgTravelTime, lgTravelTime], ylims, 'm--', label='Lg Phase')
    
    plt.title(' '.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
                        tr.stats.station, tr.stats.channel)))
    plt.xlabel('Time Since Earthquake (s)')
    plt.legend(fontsize=10)
    
    # select p, s and Rg
    picks = asarray(plt.ginput(3,timeout=-1))

    # get x indicies for fft
    x1 = int((picks[0,0]-reftime) * tr.stats.sampling_rate)
    x2 = int((picks[1,0]-reftime) * tr.stats.sampling_rate)
    x3 = int((picks[2,0]-reftime) * tr.stats.sampling_rate)
    
    plt.close(fig)
    
    return picks, x1, x2, x3

###############################################################################
# set velocity mdel
###############################################################################
model = TauPyModel(model="iasp91")

###############################################################################
# parse event file
###############################################################################

wacatcsv = 'lake_muir_events.csv'

lines =open(wacatcsv).readlines()[1:]
    
#initiate event directory
evdict = []
    
for line in lines:
    dat = line.strip().split(',')
    try:
        dateTime = datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%f')
    except:
        dateTime = datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S')
        
    tdict = {'datetime': dateTime, \
             'year': dateTime.year, 'month': dateTime.month, \
             'day': dateTime.day, 'hour': dateTime.hour, \
             'minute': dateTime.minute, 'second': dateTime.second, \
             'lat': float(dat[2]), 'lon': float(dat[3]), \
             'dep': float(dat[4]), 'mag': float(dat[1])} 

    evdict.append(tdict)

###############################################################################
# parse records
###############################################################################
outtxt = ''
"""records = 'preferred_records_edit.csv'
mseedfiles = open(records).readlines(dat[9])[0:]"""
'''
ranges = []
mseedfiles = []
lines =  open('preferred_records_edit.csv').readlines()[:]
for row in lines:
		mseedfiles.append(row.strip().split(',')[9])
		ranges.append(float(row.strip().split(',')[5]))
'''

mseedfiles = listdir_extension('data', 'mseed')
m = 0
for mseedfile in mseedfiles[m:]:
    
    # read mseed
    st = read(path.join('data', mseedfile))
    
    # remove junk channels
    st = remove_low_sample_data(st)
    
###############################################################################
# associate event and get distance
###############################################################################

    evFound = False
    for evnum, ev in enumerate(evdict): 
	
        '''
        if st[0].stats.starttime > UTCDateTime(ev['datetime']-timedelta(seconds=360)) \
           and st[0].stats.starttime < UTCDateTime(ev['datetime']+timedelta(seconds=60)):
        '''
        print UTCDateTime(ev['datetime'])
        print st[0].stats.starttime
        if UTCDateTime(ev['datetime']) > st[0].stats.starttime  \
           and UTCDateTime(ev['datetime']) < st[0].stats.endtime:
            evFound = True
            eqlo = ev['lon']
            eqla = ev['lat']
            eqmag = ev['mag']
            eqdp = ev['dep']
            eqdt = ev['datetime']            
    
    if evFound == True:
        # get station details
        print 'Getting picks for', mseedfile
        
        sta_data = return_sta_data(st[0].stats.station)
        rngkm, azim, baz = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])
        rngdeg = km2deg(rngkm)
        
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
        
        # estimate Rg
        rgTravelTime = rngkm/3.05
        
###############################################################################
# plot
###############################################################################
        
        # plot verticals only
        x1 = nan
        
        for tr in st:
            if len(st.traces) > 3:
                print tr.stats.channel
                if tr.stats.channel.endswith('Z') and tr.stats.channel[1].encode('ascii','ignore') != 'P' \
                   and tr.stats.channel[1].encode('ascii','ignore') != 'N':
                    picks, x1, x2, x3 = do_picks(tr, pTravelTime, sTravelTime, lgTravelTime, rgTravelTime)
                    
            elif len(st.traces) == 3:
                if tr.stats.channel.endswith('Z') and tr.stats.channel[1].encode('ascii','ignore') != 'P':
                    picks, x1, x2, x3 = do_picks(tr, pTravelTime, sTravelTime, lgTravelTime, rgTravelTime)
                    
            elif len(st.traces) < 3:
                if tr.stats.channel.endswith('Z') and tr.stats.channel[1].encode('ascii','ignore') != 'P' \
                   and tr.stats.channel[1].encode('ascii','ignore') != 'N':
                    picks, x1, x2, x3 = do_picks(tr, pTravelTime, sTravelTime, lgTravelTime, rgTravelTime)
            
        # capture data proc file
        if not isnan(x1):
            outtxt = ','.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
    					          str(eqlo), str(eqla), str(eqdp), str(eqmag), str('%0.2f' % rngkm), str('%0.2f' % azim), \
    						       str(st[0].stats.sampling_rate), \
                                str('%0.3f' % picks[0,0]), str('%0.3f' % picks[1,0]), str('%0.3f' % picks[2,0]), \
                                str(x1), str(x2), str(x3), \
                                path.join(getcwd(), mseedfile), eqdt.strftime('%Y-%m-%dT%H:%M:%S.%f')))
            
            # only write if x3 > x1
            if x3 > x1:
                # save processing file
                recfile = path.split(mseedfile)[-1][:-5]+'picks'
                outfile = path.join('picks',recfile)
                f = open(outfile, 'wb')
                f.write(outtxt)
                f.close()
        
    elif evFound == False:
       print 'Cannot associate event for:', mseedfile
       
    m += 1
    
'''        
        outtxt += ','.join((tr.stats.starttime.strftime('%Y-%m-%dT%H:%M:%S.%f'), \
					           str(eqlo), str(eqla), str(eqdp), str(eqmag), str(repi), str(azim), \
						        str(st[0].stats.sampling_rate), mseedfile)) + '\n'
'''
###############################################################################
# set arrival times
###############################################################################

                            

    