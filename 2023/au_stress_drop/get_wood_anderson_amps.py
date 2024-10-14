from obspy import read, UTCDateTime, Trace
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from data_fmt_tools import remove_low_sample_data, return_sta_data, remove_acceleration_data, fix_stream_channels_bb2sp
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay
from io_catalogues import parse_ga_event_query
from mapping_tools import distance
from os import path, chmod, stat, getcwd
from numpy import arange, sqrt, pi, exp, log, logspace, interp, nan, where, isnan, nanmean
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

def calc_disp_wave(tr, corfftr, corffti, freq):
    from numpy import fft, pi
    
    if tr.stats.channel[1] == 'N':
        dispfftr = corfftr / ((2 * pi * abs(freq))**2)
        dispffti = corffti / ((2 * pi * abs(freq))**2)
    else:
        dispfftr = corfftr / (2 * pi * abs(freq))
        dispffti = corffti / (2 * pi * abs(freq))
        
    n = len(corfftr)
    freq[0] = 1.0
        
    dispfftr[0] = 0
    dispffti[0] = 0
    freq[0] = 0
    complex_array = dispfftr + 1j*dispffti

    idisp = fft.ifft(complex_array,n)

    return idisp

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

def get_wa_amp(displacement_trace, sensitivity, start_idx, end_idx):
    from obspy.signal.invsim import corn_freq_2_paz
    
    corn_freq = 1.25
    
    # get W-A paz as per Uhrhammer & Collins (1990)    
    if sensitivity == 2800:
        damping = 0.8
    elif sensitivity == 2080:
        damping = 0.7
    
    paz = corn_freq_2_paz(corn_freq, damp=damping)
    
    # claculate the displacement amplitude from Wood-Anderson seismometer
    paz_wa = {'sensitivity': sensitivity, 'zeros': paz['zeros'], 'gain': 1,
              'poles': paz['poles']}
    disp_wa = displacement_trace.copy().simulate(paz_remove=None,
                                                 paz_simulate=paz_wa, water_level=10)

    try:
        ampl = max(abs(disp_wa.data[start_idx:end_idx]))
    except:
        ampl = nan

    return ampl * 1000.

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
        
        # change channel to EHZ - bit of a kluge     
        if pazfile == 'NULL' and netid == 'OZ' and tr.stats.channel == 'HHZ':
            #fix_stream_channels_bb2sp('iris_dump/'+mseedfile)
            print('    Changing channel code ...')
            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                 = get_response_info(tr.stats.station, recdate.datetime, 'EHZ', tr.stats.network)
         
        if pazfile == 'NULL' and tr.stats.network == '':
            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
              = get_response_info(tr.stats.station, recdate.datetime, tr.stats.channel, 'OZ')
            
        if tr.stats.network == '' or tr.stats.network == 'AB':
            if pazfile == 'NULL':
                nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                  = get_response_info(tr.stats.station, recdate.datetime, 'EHZ', 'OZ')        
        
        # get fft of trace
        freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
        mi = (len(freq)/2)
        mi = int(round(mi))
        
        # get response for given frequencies
        real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                            gain, inst_ty)
        
        # deconvolve response
        corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
        
        dispwave = calc_disp_wave(tr, corfftr, corffti, freq)
        
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
        elif tr.stats.network == 'M8':
            paz = dm8_parser.get_response(seedid,start_time)
            staloc = dm8_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '5C':
            paz = d5c_parser.get_response(seedid,start_time)
            staloc = d5c_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == '3O':
            paz = d3o_parser.get_response(seedid,start_time)
            staloc = d3o_parser.get_coordinates(seedid,start_time)
        elif tr.stats.network == 'AM':
            paz = dam_parser.get_response(seedid,start_time)
            staloc = dam_parser.get_coordinates(seedid,start_time)
        
        '''
        elif tr.stats.network == 'G':
            paz = g_parser.get_paz(seedid,start_time)
            staloc = g_parser.get_coordinates(seedid,start_time)
        '''
        # simulate response
        if tr.stats.network == '2P':
            tr.remove_response(inventory=d2p_parser)
        elif tr.stats.network == '5C':
            tr.remove_response(inventory=d5c_parser)
        elif tr.stats.network == 'M8':
            tr.remove_response(inventory=dm8_parser)
        elif tr.stats.network == 'AM':
            tr.remove_response(inventory=dam_parser)
        elif tr.stats.network == '3O':
            tr.remove_response(inventory=d3o_parser)
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
                
        dispwave = calc_disp_wave(tr, corfftr, corffti, freq)
                
    return dispwave

def retry_stationlist_fft(tr, pickDat):
    seedid = tr.get_id()
    channel = tr.stats.channel
    start_time = tr.stats.starttime
    
    recdate = pickDat['origintime']
    
    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
          = get_response_info(tr.stats.station, recdate.datetime, tr.stats.channel, tr.stats.network)
    
    # change channel to EHZ - bit of a kluge     
    if pazfile == 'NULL' and netid == 'OZ' and tr.stats.channel == 'HHZ':
        #fix_stream_channels_bb2sp('iris_dump/'+mseedfile)
        print('    Changing channel code ...')
        nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
             = get_response_info(tr.stats.station, recdate.datetime, 'EHZ', tr.stats.network)
     
    if pazfile == 'NULL' and tr.stats.network == '':
        nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
          = get_response_info(tr.stats.station, recdate.datetime, tr.stats.channel, 'OZ')
          
    if tr.stats.network == '' or tr.stats.network == 'AB':
            if pazfile == 'NULL':
                nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                  = get_response_info(tr.stats.station, recdate.datetime, 'EHZ', 'OZ')
                  
    # get fft of trace
    freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
    mi = (len(freq)/2)
    mi = int(round(mi))
    
    # get response for given frequencies
    real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                        gain, inst_ty)
    
    # deconvolve response
    corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
    
    dispwave = calc_disp_wave(tr, corfftr, corffti, freq)
    
    staloc = {'latitude':stla, 'longitude':stlo}
    	
    return dispwave
    	
################################################################################
# find eevents - should have done this in pick files
################################################################################

evdict = parse_ga_event_query('au_ge_4.4_earthquakes_export_edit.csv')

def get_ev_deets(fft_datetime):
    for evnum, ev in enumerate(evdict): 
        #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
               magtype = ev['magType']
               evname = ev['description']
               mag = ev['mag'] # pref_mag
            
               evdat = {'mag': ev['mag'], 'magtype': ev['magType'], 'place': ev['description'],
                        'eqla': ev['lat'], 'eqlo': ev['lon'], 'eqdep': ev['dep']}
               
    return evdat


################################################################################
# get pick files
################################################################################
folder = 'record_picks'
#folder = 'new_picks' # for testing
pickfiles = listdir_extension(folder, 'picks')

################################################################################
# set some defaults
################################################################################

interp_freqs = logspace(-1.5,2,176)[:-26] # from 0.03-30 Hz

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
    #print('test parsers')
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
    dm8_parser = read_inventory('/Users/trev/Documents/Networks/AUSPASS/m8-inventory.xml')
    d2p_parser = read_inventory('/Users/trev/Documents/Networks/AUSPASS/2p-inventory-edit.xml')
    d5c_parser = read_inventory('/Users/trev/Documents/Networks/AUSPASS/5c-inventory.xml')
    d3o_parser = read_inventory('/Users/trev/Documents/Networks/AUSPASS/3o-inventory.xml')
    dam_parser = read_inventory('/Users/trev/Documents/Networks/AM/R7AF5.xml')

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
        except:
            try:
                mseedfile = path.split(pickDat['mseed_path'])[-1]
                st = read(path.join('iris_dump', mseedfile))
            except:
                print('Skipping: '+pickDat['mseed_path'])
                skipRec = True
        
        if skipRec == False:
            
            # split trace containing gaps into contiguous unmasked traces
            st = st.split()
            
            # remove low sample rate data
            new_st = remove_low_sample_data(st)
            
            # remove HTT stations
            if new_st[0].stats.starttime > UTCDateTime(2018,12,12) and new_st[0].stats.starttime < UTCDateTime(2023,4,30):
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
                if tr.stats.channel.endswith('Z') and tr.stats.sampling_rate >= 10.0:
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
                        lofreq = 0.02
                    else:
                        lofreq = 0.2
                    hifreq = 0.475 * tr.stats.sampling_rate
                    
                    # get picks
                    xdat = range(0, tr.stats.npts)
                    #plt.plot(xdat, tr.data, 'b-', lw=0.5)
                    pidx = pickDat['pidx']
                    sidx = pickDat['sidx']
                    eidx = pickDat['eidx']
                    
                    #####################################################################
                    # now correct to disp wave
                    #####################################################################
                    
                    # get instrument corrected spectra
                    try:
                        dispwave = response_corrected_fft(tr, pickDat)
                    except:
                        dispwave = retry_stationlist_fft(tr, pickDat)
                    
                    disp_tr = Trace(data=dispwave.real)
                    disp_tr.times = tr.times
                    disp_tr.stats.sampling_rate = tr.stats.sampling_rate
                    disp_tr.stats.starttime = tr.stats.starttime
                    
                    #####################################################################
                    # loop thru WA coeffs and filters
                    #####################################################################
                    
                    # loop through magnification and damping
                    wa_sensitivities = [2080, 2800]
                    hipass_filters = [lofreq, 0.2, 0.5, 0.75]
                    
                    for wa_sensitivity in wa_sensitivities:
                        for hpf in hipass_filters:
                            
                            # pre-filter disp wave
                            tr_filt = disp_tr.copy()
                            tr_filt.filter('highpass', freq=hpf, corners=4, zerophase=True)
    
                            # get W-A amp
                            start_idx = pickDat['sidx']
                            end_idx = pickDat['eidx']
                            wa_amp = get_wa_amp(tr_filt, wa_sensitivity, start_idx, end_idx)
                    
                            tr_key = '_'.join(('wa_amp', str(wa_sensitivity), str(hpf)))
            
                            traceDat[tr_key] = wa_amp
                    
                    #####################################################################
                    # add trace data to recDat
                    #####################################################################
                    recDat['wa_data'] = traceDat
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
            evdat = get_ev_deets(UTCDateTime(pickDat['evdt']))
            recDat['mag'] = evdat['mag']
            recDat['magType'] = evdat['magtype']
            recDat['eqlo'] = evdat['eqlo']
            recDat['eqla'] = evdat['eqla']
            recDat['eqdp'] = evdat['eqdep']
            recDat['place'] = evdat['place']
            recDat['evdt'] = pickDat['evdt']
            
            # get sta data
            staDat = return_sta_data(tr.stats.station)
            recDat['stlo'] = staDat['stlo']
            recDat['stla'] = staDat['stla']
            
            # calc new distance from event list
            recDat['repi'] = distance(evdat['eqla'], evdat['eqlo'], staDat['stla'], staDat['stlo'])[0]
            recDat['rhyp'] = sqrt(recDat['repi']**2 + recDat['eqdp']**2)
            
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

pklfile = open('wa_data.pkl', 'wb')
pickle.dump(records, pklfile, protocol=-1)
pklfile.close()


# records[1][channels[0]]  


# set permissions for all to execute
#chmod('nac_fft_data.pkl', 0o777)
        







