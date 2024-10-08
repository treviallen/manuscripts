from obspy import read, UTCDateTime
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from data_fmt_tools import remove_low_sample_data, return_sta_data, remove_acceleration_data
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay
from io_catalogues import parse_ga_event_query
from os import path, chmod, stat, getcwd
from numpy import array, arange, sqrt, pi, exp, log, logspace, interp, nan, where, isnan, nanmean, unique
from datetime import datetime, timedelta
import pickle
import warnings
warnings.filterwarnings("ignore")
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.style.use('classic')

# script to remove non trillium chans
def remove_htt(st):
    cnt = 0
    
    for tr in st:
        if tr.stats.channel == 'BHZ':
            cnt += 1
        
    # remove nontrillium 
    if cnt == 2:
        for tr in st:
           if not tr.stats.location == '10':
               st.remove(tr)
               
    return st

#plt.ion()
def parse_pickfile(pickfile):
    from numpy import nan
    
    # parse pick file
    line = open(path.join(folder, pickfile)).read()
    
    data = line.strip().split(',')
    
    if not data[0].startswith('junk'):
        starttime = UTCDateTime(datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f'))
        origintime = UTCDateTime(datetime.strptime(data[1], '%Y-%m-%dT%H:%M:%S.%f'))
        recdate = datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f')
        evdate = datetime.strptime(data[1], '%Y-%m-%dT%H:%M:%S.%f')
        
        # recalculate distance with improved locations
        #recalc dists!
        
        pickDat = {'starttime': starttime, 'origintime': origintime, \
                   'ev':data[1][0:16].replace(':','.'), 'evdt':UTCDateTime(evdate),
                   'eqlo': float(data[2]), 'eqla': float(data[3]),
                   'eqdp': float(data[4]), 'mag': float(data[5]), 'rhyp': float(data[6]),
                   'azim': float(data[7]), 'sps': float(data[8]), \
                   'ch1': data[9], 'ch2': data[10], 'ch3': data[11], 
                   'ppk': float(data[12]), 'spk': float(data[13]), 'epk': float(data[14]), \
                   'pidx': int(data[15]), 'sidx': int(data[16]), 'eidx': int(data[17]), 'mseed_path': data[18]}
                   	
    else:
        pickDat = {'mag': nan}
        
    return pickDat

def calc_disp_spectra(tr, corfftr, corffti, freq):
    if tr.stats.channel[1] == 'N':
        dispamp = sqrt((tr.stats.delta)**2 * (corfftr**2 + corffti**2)) / ((2 * pi * freq)**2)
    else:
        dispamp = sqrt((tr.stats.delta)**2 * (corfftr**2 + corffti**2)) / (2 * pi * freq)
        
    return dispamp

def get_smoothed_fft_spectra(freqs, disp_amps):
    
    # smooth spectra
    if len(freqs) > 10000:
        sw = 201
    elif len(freqs) > 5000:
        sw = 101
    elif len(freqs) > 1000:
        sw = 51
    elif len(freqs) > 500:
        sw = 21
    else:
        sw = 5
    
    # smooth spectra
    if sw > 11:
        smoothed_disp = exp(savitzky_golay(log(disp_amps), sw, 3))
    else:
        smoothed_disp = exp(savitzky_golay(log(disp_amps), sw, 2))
        
    # log-log interpolate to "interp_freqs"
    smoothed_interp_disp = exp(interp(log(interp_freqs), log(freqs), log(smoothed_disp), \
                               left=nan, right=nan))
    
    return smoothed_disp, smoothed_interp_disp

def response_corrected_fft(tr, pickDat):
    import matplotlib.pyplot as plt
    from numpy import fft, sqrt
    
    seedid = tr.get_id()
    channel = tr.stats.channel
    start_time = tr.stats.starttime
    
    # check if HSR AU data
    use_stationlist = False
    if tr.stats.network == 'AU' and tr.stats.channel.startswith('HH'):
        use_stationlist = True
    elif tr.stats.network == 'AU' and tr.stats.channel.startswith('BH'):
        use_stationlist = True
    elif tr.stats.network == 'AU' and tr.stats.start_time.year >= 2017:
        use_stationlist = True
    elif tr.stats.network == 'OA' and tr.stats.channel.startswith('HH'):
        use_stationlist = True
    elif tr.stats.network == '2O' and tr.stats.channel.startswith('HH'):
        use_stationlist = True
    elif tr.stats.network == 'AU' and tr.stats.channel.startswith('EH'):
        use_stationlist = True
    elif tr.stats.network == 'AU' and tr.stats.channel.startswith('SH'):
        use_stationlist = True
    elif tr.stats.network == '' and tr.stats.channel.startswith('EN'): # for DRS
        use_stationlist = True  
    elif tr.stats.network == '7M':
        use_stationlist = True                       
    elif tr.stats.network == '5J':
        use_stationlist = True       
    elif tr.stats.station == 'AS32' or tr.stats.station == 'ARPS' or tr.stats.station == 'ARPS' or tr.stats.network == 'MEL': 
        use_stationlist = True
    #print('use_stationlist', use_stationlist) 
       
    if use_stationlist == True:
        #recdate = datetime.strptime(ev['timestr'], "%Y-%m-%dT%H:%M:%S.%fZ")
        recdate = pickDat['origintime']
        #print(seedid, channel)
        nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
              = get_response_info(tr.stats.station, recdate.datetime, tr.stats.channel, tr.stats.network)
        
        # get fft of trace
        freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
        mi = (len(freq)/2)
        mi = int(round(mi))
        
        # get response for given frequencies
        real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                            gain, inst_ty)
        
        # deconvolve response
        corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
        
        '''
        # from gmprocess
        dt = trace.stats.delta
        spec = abs(np.fft.rfft(trace.data, n=nfft)) * dt
        freqs = np.fft.rfftfreq(nfft, dt)
        return spec, freqs
        '''
        
        dispamp = calc_disp_spectra(tr, corfftr, corffti, freq)
        
        staloc = {'latitude':stla, 'longitude':stlo}
        
    # just use IRIS dataless volume - much easier, but less transparent!
    else:
        if tr.stats.network == 'AU':
            try:
                paz = au_parser.get_paz(seedid,start_time)
                staloc = au_parser.get_coordinates(seedid,start_time)
            except:
                paz = cwb_parser.get_paz(seedid,start_time)
                staloc = cwb_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == 'GE':
            paz = ge_parser.get_paz(seedid,start_time)
            staloc = ge_parser.get_coordinates(seedid,start_time)
        
        elif tr.stats.network == 'IU':
            paz = iu_parser.get_paz(seedid,start_time)
            staloc = iu_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == 'II':
            paz = ii_parser.get_paz(seedid,start_time)
            staloc = ii_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == 'S1':
            paz = s1_parser.get_paz(seedid,start_time)
            staloc = s1_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '1K':
            paz = d1k_parser.get_paz(seedid,start_time)
            staloc = d1k_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '1H':
            paz = d1h_parser.get_paz(seedid,start_time)
            staloc = d1h_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '1P':
            paz = d1p_parser.get_paz(seedid,start_time)
            staloc = d1p_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '1Q':
            paz = d1q_parser.get_paz(seedid,start_time)
            staloc = d1q_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '6F':
            paz = d6f_parser.get_paz(seedid,start_time)
            staloc = d6f_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7B':
            paz = d7b_parser.get_paz(seedid,start_time)
            staloc = d7b_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7D':
            paz = d7d_parser.get_paz(seedid,start_time)
            staloc = d7d_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7G':
            paz = d7g_parser.get_paz(seedid,start_time)
            staloc = d7g_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7H':
            paz = d7h_parser.get_paz(seedid,start_time)
            staloc = d7h_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7I':
            paz = d7i_parser.get_paz(seedid,start_time)
            staloc = d7i_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7J':
            paz = d7j_parser.get_paz(seedid,start_time)
            staloc = d7j_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7K':
            paz = d7k_parser.get_paz(seedid,start_time)
            staloc = d7k_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7M':
            try:
                paz = d7m_parser.get_paz(seedid,start_time)
                staloc = d7m_parser.get_coordinates(seedid,start_time)
            except:
                paz = d7n_parser.get_paz(seedid,start_time)
                staloc = d7n_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7S':
            paz = d7s_parser.get_paz(seedid,start_time)
            staloc = d7s_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '7T':
            paz = d7t_parser.get_paz(seedid,start_time)
            staloc = d7t_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '8K':
            paz = d8k_parser.get_paz(seedid,start_time)
            staloc = d8k_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '2P':
            paz = d2p_parser.get_response(seedid,start_time)
            staloc = d2p_parser.get_coordinates(seedid,start_time)
        '''
        elif tr.stats.network == 'G':
            paz = g_parser.get_paz(seedid,start_time)
            staloc = g_parser.get_coordinates(seedid,start_time)
        '''
        # simulate response
        if tr.stats.network == '2P':
            tr.remove_response(inventory=d2p_parser)
        else:
            if tr.stats.channel.endswith('SHZ') or tr.stats.channel.endswith('EHZ'):
                tr = tr.simulate(paz_remove=paz, water_level=10) #  testing water level for SP instruments
            else:
                tr = tr.simulate(paz_remove=paz)
        
        # get fft
        freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
        mi = (len(freq)/2)
        mi = int(round(mi))
        
        corfftr = wavfft.real
        corffti = wavfft.imag
                
        dispamp = calc_disp_spectra(tr, corfftr, corffti, freq)
                
    return freq[1:mi], dispamp[1:mi]

def retry_stationlist_fft(tr, pickDat):
    seedid = tr.get_id()
    channel = tr.stats.channel
    start_time = tr.stats.starttime
    
    recdate = pickDat['origintime']
    
    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
          = get_response_info(tr.stats.station, recdate.datetime, tr.stats.channel, tr.stats.network)
    
    # get fft of trace
    freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
    mi = (len(freq)/2)
    mi = int(round(mi))
    
    # get response for given frequencies
    real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                        gain, inst_ty)
    
    # deconvolve response
    corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
    
    dispamp = calc_disp_spectra(tr, corfftr, corffti, freq)
    
    staloc = {'latitude':stla, 'longitude':stlo}
    	
    return freq[1:mi], dispamp[1:mi]
   	
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
               'dep': float(dat[4]), 'mag': float(dat[5]),
               'magType': 'ML'}
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
            # get pref mag and mag type
            if dat[7].startswith('ML'):
                mag = float(dat[5])
                magType = 'ML'
            else:
                 mag = float(dat[6])
                 magType = 'MW'
                 
            tmpdict = {'id': str(dat[0]), 'datetime': UTCDateTime(dat[0]),
                       'lat': float(dat[2]), 'lon': float(dat[1]), 
                       'dep': float(dat[3]), 'mag': mag, 'magType': 'MW'}
    agmdevs.append(tmpdict)

def get_ev_deets(fft_datetime):
    evFound = False
    mag = nan
    eqmag = nan
    eqmagType = ''
    magType = ''
    
    for evnum, ev in enumerate(evdict): 
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=60)) \
            and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=180)):
        
            evFound = True
            eqlo = ev['lon']
            eqla = ev['lat']
            eqmag = ev['mag']
            eqmagType = ev['magType']
            eqdp = ev['dep']
            eqdt = ev['datetime']
            print(ev)
            gaid = ev['gaid']
            
            mag = eqmag
            magType = eqmagType
            
    if evFound == False:
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
                    eqmagType = ev['magType']
                    eqdp = ev['dep']
                    eqdt = ev['datetime']
                    gaid = ''
                    
            else:
                if UTCDateTime(ev['datetime']) > st[0].stats.starttime-timedelta(seconds=30)  \
                   and UTCDateTime(ev['datetime']) < st[0].stats.endtime:
                    evFound = True
                    eqlo = ev['lon']
                    eqla = ev['lat']
                    eqmag = ev['mag']
                    eqmagType = ev['magType']
                    eqdp = ev['dep']
                    eqdt = ev['datetime']
                    gaid = ''
                    
                    
        mag = eqmag
        magType = eqmagType
           
    return mag, magType


################################################################################
# get pick files
################################################################################
folder = 'record_picks'
#folder = 'new_picks' # for testing
pickfiles = listdir_extension(folder, 'picks')

################################################################################
# set some defaults
################################################################################

interp_freqs = logspace(-1.,2,176)[:-26] # from 0.1-30 Hz

################################################################################
# parse AU dataless
################################################################################

# read dataless seed volumes
print('\nReading dataless seed volumes...')
from obspy.io.xseed import Parser
from obspy import read_inventory

if getcwd().startswith('/nas'):
    au_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU.IRIS.dataless')
    cwb_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU.cwb.dataless')
    s1_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/S1/S1.IRIS.dataless')
    #ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE1.IRIS.dataless')
    #ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE2.IRIS.dataless')
    iu_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/IU.IRIS.dataless')
    ii_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/II/II.IRIS.dataless')
    
else:
    au_parser = Parser('/Users/trev/Documents/Networks/AU/AU.IRIS.dataless')
    cwb_parser = Parser('/Users/trev/Documents/Networks/AU/AU.cwb.dataless')
    s1_parser = Parser('/Users/trev/Documents/Networks/S1/S1.IRIS.dataless')
    iu_parser = Parser('/Users/trev/Documents/Networks/IU/IU.IRIS.dataless')
    #g_parser = Parser('/Users/trev/Documents/Networks/G/G.IRIS.dataless')
    ii_parser = Parser('/Users/trev/Documents/Networks/II/II.IRIS.dataless')
    d1h_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/1H_EAL2_2010.dataless')
    d1k_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/1K_ALFREX_2013.dataless')
    d1p_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/1P_BASS_2011.dataless')
    d1q_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/1Q_AQT_2016.dataless')
    d6f_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/6F_BILBY_2008.dataless')
    d7b_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7B_SKIPPY_1993.dataless')
    d7d_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7D_KIMBA97_1997.dataless')
    d7g_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7G_WACRATON_2000.dataless')
    d7h_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7H_TIGGERBB_2001.dataless')
    d7i_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7I_TASMAL_2003.dataless')
    d7j_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7J_CAPRAL_2005.dataless')
    d7k_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7K_SOC_2007.dataless')
    d7m_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7M_COPA_2014.dataless')
    d7n_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7M_MALTLACHLAN_1998.dataless')
    d7s_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7S_SETA_2006.dataless')
    d7t_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/7T_SEAL2_2007.dataless')
    d8k_parser = Parser('/Users/trev/Documents/Networks/AUSPASS/8K_CAPRICORNHPS_2014.dataless')
    d2p_parser = read_inventory('/Users/trev/Documents/Networks/AUSPASS/2p-inventory-edit.xml')

################################################################################
# loop through pick files
################################################################################
records = [] 
f = 0 
start_idx = 0
#pickfiles = ['2023-01-05T05.08.AU.ONGER.picks']
for p, pf in enumerate(pickfiles[start_idx:]):
    skipRec = False
    tr = nan
    recDat = {}
    
    pickDat = parse_pickfile(pf)
    
    if isnan(pickDat['mag']) == False: # and pf == '1997-03-05T06.15.00.AD.WHY.picks':
        
        channels = []
        if not pickDat['ch1'] == '':
            channels.append(pickDat['ch1'])
        if not pickDat['ch2'] == '':
            channels.append(pickDat['ch2'])
        if not pickDat['ch3'] == '':
            channels.append(pickDat['ch3'])
        
        try:
            # look and file size and skip if too big    
            #print(pickDat['mseed_path'])
            #fileStat = stat(pickDat['mseed_path'])
            #print('filesize', fileStat.st_size
            #mseedfile = path.join('iris_dump', path.split(pickDat['mseed_path'])[-1])
            st = read(pickDat['mseed_path'])
            fullmseed = pickDat['mseed_path']
        except:
            try:
                mseedfile = path.split(pickDat['mseed_path'])[-1]
                st = read(path.join('auspass_dump', mseedfile))
                fullmseed = path.join('auspass_dump', mseedfile)
            except:
                try:
                    # fix year in path
                    pathsplit = pickDat['mseed_path'].split(path.sep)
                    newp1 = path.sep.join(pathsplit[:-4])
                    newp2 = path.sep.join(pathsplit[-3:])
                    newpath = path.join(newp1,'2023',newp2).replace('mseed_dump','iris_dump')
                    st = read(newpath)
                    fullmseed = newpath
                except:
                    print('Skipping: '+pickDat['mseed_path'])
                    skipRec = True
        
        if skipRec == False:
            
            # first - fix S1 network
            if st[0].stats.network == 'S':
                for i in range(0, len(st)):
                    st[i].stats.network = 'S1'
                    
                st.write(fullmseed, format='MSEED')
            
            # split trace containing gaps into contiguous unmasked traces
            st = st.split()
            
            # remove low sample rate data
            new_st = remove_low_sample_data(st)
            
            # remove HTT stations
            if new_st[0].stats.starttime > UTCDateTime(2018,12,12):
                new_st = remove_htt(new_st)
            
            # remove acceleration data
            #new_st = remove_acceleration_data(new_st)
            
            # purge unused traces
            for tr in new_st:
                removeTrace = True
                for chan in channels:
                    if chan == tr.stats.channel:
                        removeTrace = False
                if removeTrace == True:
                    new_st = new_st.remove(tr)
                
            new_st.merge()
            
            if len(new_st) > 3:
                new_st = new_st[0:3]
            
            #st_filt = new_st.copy()
            #sidx = int(round(0.05*st_filt[-1].stats.npts))
            #eidx = int(round(0.95*st_filt[-1].stats.npts))
            
            #st_filt.filter('bandpass', freqmin=0.5, freqmax=10, corners=2, zerophase=True)
                    
            print('\n'+str(p)+' Reading mseed file:', path.split(pickDat['mseed_path'])[-1])
                
            #####################################################################
            # loop thru traces
            #####################################################################
            chandict = []
            for tr in new_st:
                # only do Z channel
                if tr.stats.channel.endswith('Z'):
                    traceDat = {}
                    seedid = tr.get_id()
                    channel = tr.stats.channel
                    start_time = tr.stats.starttime 
                    
                    # demean and taper data  
                    tr_proc = tr.copy()
                    tr_proc.detrend(type='demean')
                    
                    # EN? dodgy stn channel from parsing JUMP data
                    if tr.stats.channel.startswith('BH') or tr.stats.channel.startswith('HH') \
                       or tr.stats.channel[1] == 'N':
                        lofreq = 0.075
                    else:
                        lofreq = 0.2
                    lofreq=0.2
                    hifreq = 0.475 * tr.stats.sampling_rate
                    
                    # get picks
                    xdat = range(0, tr.stats.npts)
                    #plt.plot(xdat, tr.data, 'b-', lw=0.5)
                    pidx = pickDat['pidx']
                    sidx = pickDat['sidx']
                    eidx = pickDat['eidx']
                    #####################################################################
                    # now do ffts
                    #####################################################################
                    
                    '''
                    # get noise fft of trace
                    '''
                    f += 1
                    sttr = tr.stats.starttime + 2. # should add 1-2 secs instead of proportion 
                    ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 5. # allow buffer
                    
                    if sttr > ettr or ettr-sttr < 10.:
                        sttr = tr.stats.starttime + 1. # should add 1-2 secs instead of proportion 
                        ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 2. # allow buffer
                        
                    if sttr > ettr or ettr-sttr < 10.:
                        sttr = tr.stats.starttime + 0.5 # should add 1-2 secs instead of proportion 
                        ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 1. # allow buffer
                    
                    noise_window = ettr-sttr
                    ntr_trim = tr_proc.copy()
                    ntr_trim.trim(sttr, ettr)
                    ntr_trim.taper(0.02, type='hann', max_length=None, side='both')
                    
                    # get instrument corrected spectra
                    try:
                        freqs, n_disp_amp = response_corrected_fft(ntr_trim, pickDat)
                    except:
                        freqs, n_disp_amp = retry_stationlist_fft(ntr_trim, pickDat)
                        
                    smoothed_disp, smoothed_interp_disp = get_smoothed_fft_spectra(freqs, n_disp_amp)
                    
                    traceDat = {'hi_freq_filt': hifreq, 'lo_freq_filt': lofreq, 
                                'noise_spec': smoothed_interp_disp, 'freqs': interp_freqs,
                                'sample_rate': tr.stats.sampling_rate}
                        
                    #plt.loglog(freqs, n_disp_amp, 'b-', lw=0.3)
                                  
                    '''
                    # get p/s-wave fft of trace
                    '''
                    
                    sttr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 10. # allow buffer
                    ettr = tr.stats.starttime + pickDat['eidx'] * tr.stats.delta + 10. # allow buffer
                    
                    pstr_trim = tr_proc.copy()
                    pstr_trim.trim(sttr, ettr)
                    pstr_trim.taper(0.02, type='hann', max_length=None, side='both')
                    
                    # get instrument corrected spectra
                    try:
                        freqs, ps_disp_amp = response_corrected_fft(pstr_trim, pickDat)
                    except:
                        freqs, ps_disp_amp = retry_stationlist_fft(pstr_trim, pickDat)
                        
                    smoothed_disp, smoothed_interp_disp = get_smoothed_fft_spectra(freqs, ps_disp_amp)
                    
                    traceDat['p-swave_spec'] = smoothed_interp_disp
                    #plt.loglog(freqs, ps_disp_amp, 'g-', lw=0.3)
                    
                    '''
                    # get s-wave fft of trace
                    '''
                    
                    sttr = tr.stats.starttime + pickDat['sidx'] * tr.stats.delta - 5. # allow buffer
                    ettr = tr.stats.starttime + pickDat['eidx'] * tr.stats.delta + 5. # allow buffer
                    
                    str_trim = tr_proc.copy()
                    str_trim.trim(sttr, ettr)
                    str_trim.taper(0.02, type='hann', max_length=None, side='both')
                    
                    # get instrument corrected spectra
                    try:
                        freqs, s_disp_amp = response_corrected_fft(str_trim, pickDat)
                    except:
                        freqs, s_disp_amp = retry_stationlist_fft(str_trim, pickDat)
                        
                    smoothed_disp, smoothed_interp_disp = get_smoothed_fft_spectra(freqs, s_disp_amp)
                    
                    traceDat['swave_spec'] = smoothed_interp_disp
                    
                    #plt.loglog(freqs, s_disp_amp, 'r-', lw=0.3)
                    #plt.xlim([0.05, 10])
                    #plt.show()
                    
                    '''
                    # get SN-Ratio
                    '''
                    sn_ratio = traceDat['swave_spec'] / traceDat['noise_spec']
                    
                    # now set frequency limits - use 1 Hz as centre
                    sn_thresh = 5.
                    
                    # find nan ratios
                    nanidx = where(isnan(sn_ratio))[0]
                    
                    # if limited noise window - set default value
                    #get mean ratio
                    mean_sn_ratio = nanmean(sn_ratio)
                    if mean_sn_ratio >= 100. or noise_window < 4.:
                        fidx = where((interp_freqs[nanidx] >= 0.25) & (interp_freqs[nanidx] < 5.))[0]
                        sn_ratio[fidx] = 9999.
                        fidx = where(interp_freqs[nanidx] < 0.25)[0]
                        sn_ratio[fidx] = 0.
                        nanidx = where(isnan(sn_ratio))[0]
                        sn_ratio[nanidx] = 0.
                        
                    else:
                        sn_ratio[nanidx] = 0.
                    
                    traceDat['sn_ratio'] = sn_ratio
                    
                    # set hi freq limit
                    fidx = where((interp_freqs >= 1) & (sn_ratio < sn_thresh))[0]
                    if len(fidx) == 0:
                        traceDat['hif_limit'] = interp_freqs[-1]
                    else:
                        traceDat['hif_limit'] = interp_freqs[fidx[0]-1]
                    
                    # set lo freq limit
                    fidx = where((interp_freqs <= 1) & (sn_ratio < sn_thresh))[0]
                    if len(fidx) == 0:
                        traceDat['lof_limit'] = interp_freqs[0]
                    else:
                        traceDat['lof_limit'] = interp_freqs[fidx[-1]+1]
                    
                    # do a manual check for good data with no noisepickDat['evdt']
                    '''
                    if sn_ratio[15] > 100.:
                        traceDat['lof_limit'] = interp_freqs[0]
                        fidx = where((interp_freqs >= 3) & (sn_ratio < sn_thresh))[0]
                        if len(fidx) == 0:
                            traceDat['hif_limit'] = interp_freqs[-1]
                        else:
                            traceDat['hif_limit'] = interp_freqs[fidx[0]-1]
                    '''
                    #####################################################################
                    # add trace data to recDat
                    #####################################################################
                    recDat[channel] = traceDat
                    #recDat = traceDat
                    #chandict.append(channel.encode('ascii', 'ignore').lower())
                    chandict.append(channel)
                        
            #####################################################################
            # populate record dictionary
            #####################################################################
            recDat['channels'] = chandict
            #recDat['sta'] = tr.stats.station.encode('ascii','ignore')
            recDat['sta'] = tr.stats.station
            recDat['net'] = tr.stats.network
            recDat['location'] = tr.stats.location
            recDat['sampling_rate'] = tr.stats.sampling_rate
            recDat['ev'] = pickDat['ev']
            recDat['eqlo'] = pickDat['eqlo']
            recDat['eqla'] = pickDat['eqla']
            recDat['eqdp'] = pickDat['eqdp']
            mag, magType, gaid = get_ev_deets(UTCDateTime(pickDat['evdt']))
            recDat['mag'] = mag
            recDat['magType'] = magType
            recDat['gaid'] = gaid
            print(mag, magType)
            #recDat['place'] = evName
            recDat['rhyp'] = pickDat['rhyp']
            recDat['eqla'] = pickDat['eqla']
            
            
            # get sta data
            staDat = return_sta_data(tr.stats.station)
            recDat['stlo'] = staDat['stlo']
            recDat['stla'] = staDat['stla']
            
            recDat['mseed_path'] = pickDat['mseed_path']
            records.append(recDat)
            
            '''
            # plt spectra
            plt.loglog(interp_freqs, recDat['p-swave_spec'])
            plt.loglog(interp_freqs, recDat['noise_spec'])
            plt.show()
            '''
            '''
            else:
                print('Large file size:', pickDat['mseed_path'])
            '''
#plt.show()        
#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('fft_data.pkl', 'wb')
pickle.dump(records, pklfile, protocol=-1)
pklfile.close()


# records[1][channels[0]]  


# set permissions for all to execute
#chmod('nac_fft_data.pkl', 0o777)
        







