from obspy import UTCDateTime
from obspy import read #, read_inventory
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay
from os import path, chmod, stat
from numpy import sqrt, pi, exp, log, logspace, interp, nan
from datetime import datetime
from data_fmt_tools import get_iris_stn_dataless_seed, get_station_distance, \
                           get_iris_data, remove_low_sample_data
import pickle

#plt.ion()

def get_fft_spectra(tr_trim, sttr, ettr, channel):
    tr_trim.taper(0.02, type='hann', max_length=None, side='both')
    freq, wavfft = calc_fft(tr_trim.data, tr_trim.stats.sampling_rate)
    mxidx = (len(freq)/2)
    
    # get fft freqs and amps
    freqs = abs(freq[1:mxidx])
    spec_amps = abs(wavfft.real[1:mxidx])
    
    # get displacement spectra
    if channel[1] == 'N': # if accelerometer
        disp_amps = spec_amps / (2.*(pi)*freqs)**2
    else: # if seismometer
        disp_amps = spec_amps / (2.*(pi)*freqs)
        
    # smooth spectra
    smoothed_disp = exp(savitzky_golay(log(disp_amps), 21, 3))
    
    # log-log interpolate to "interp_freqs"
    smoothed_interp_disp = exp(interp(log(interp_freqs), log(freqs), log(smoothed_disp), \
                               left=nan, right=nan))
    
    return tr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp

################################################################################
# get pick files
################################################################################
#folder = 'test_picks'
folder = 'picks'
pickfiles = listdir_extension(folder, 'picks')

################################################################################
# set some defaults
################################################################################

interp_freqs = logspace(-1,2,31)[:-3] # from 0.1-50 Hz

################################################################################
# parse AU dataless
################################################################################

# read dataless seed volumes
print '\nReading dataless seed volumes...'
from obspy.io.xseed import Parser
'''
au_parser = Parser(path.join('..','..','..','Networks','AU', 'AU.IRIS.dataless'))
iu_parser = Parser(path.join('..','..','..','Networks','IU', 'IU.IRIS.dataless')) 
s_parser = Parser(path.join('..','..','..','Networks','S', 'S.IRIS.dataless')) # may need to newer file here
'''
print 'Parsing dataless seed volumes...'
au_dataless = get_iris_stn_dataless_seed('AU')
s_dataless = get_iris_stn_dataless_seed('S')
iu_dataless = get_iris_stn_dataless_seed('IU')


################################################################################
# loop through pick files
################################################################################
records = [] 
f = 0                 
for p, pf in enumerate(pickfiles[0:]):
    # parse pick file
    #pf = '2002-06-20T03.54.AU.CA3.A.picks'
    line = open(path.join(folder, pf)).read()
    
    data = line.strip().split(',')
    
    utcDT = UTCDateTime(datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f'))
    recdate = datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f')
    
    pickDat = {'datetime': utcDT, 'eqlo': float(data[1]), 'eqla': float(data[2]),
               'eqdp': float(data[3]), 'mag': float(data[4]), 'rhyp': float(data[5]),
               'azim': float(data[6]), 'sps': float(data[7]), 'ppk': float(data[8]),
               'spk': float(data[9]), 'rgpk': float(data[10]), 'pidx': int(data[11]),
               'sidx': int(data[12]), 'rgidx': int(data[13]), 'mseed_path': data[14]}
    
    # look and file size and skip if too big
    fileStat = stat(pickDat['mseed_path'])
    #print 'filesize', fileStat.st_size
    
    # if file less than 10MB
    if fileStat.st_size < 1300000000L:
        #try:
        st = read(pickDat['mseed_path'])
        
        # split trace containing gaps into contiguous unmasked traces
        st = st.split()
        
        # remove low sample rate data
        new_st = remove_low_sample_data(st)
        new_st.merge()
        
        st_filt = new_st.copy()
        sidx = int(round(0.05*st_filt[-1].stats.npts))
        eidx = int(round(0.90*st_filt[-1].stats.npts))
        
        st_filt.filter('bandpass', freqmin=0.5, freqmax=10, corners=2, zerophase=True)
                
        print '\nReading mseed file:', path.split(pickDat['mseed_path'])[-1]
                
        #####################################################################
        # remove instrument response
        #####################################################################
        recDat = {'pick_file': pickDat}
        chandict = []
        addRecDat = False
        for tr in new_st:
            # if not network, set to AU
            if len(tr.stats.network) == 0:
                tr.stats.network = 'AU'
                
            seedid = tr.get_id()
            channel = tr.stats.channel
            start_time = tr.stats.starttime 
            
            # demean and taper data
            tr = tr.detrend(type='demean')
            tr = tr.taper(0.02, type='hann', max_length=None, side='both')
            
            # EN? dodgy stn channel from parsing JUMP data
            if tr.stats.channel.startswith('BH') or tr.stats.channel.startswith('HH') \
               or tr.stats.channel[1] == 'N':
                lofreq = 0.075
            else:
                lofreq = 0.2
            
            hifreq = 0.45 * tr.stats.sampling_rate
            
            '''
            tr = tr.filter('bandpass', freqmin=lofreq, freqmax=hifreq, \
                            corners=2, zerophase=True)
            '''
            #####################################################################
            # get instrument response info
            #####################################################################
            # skip out if no station data available
            try:
                
                # either use dataless or stationlist
                try:
                    paz = []
                    # get poles & zeros from dataless seed
                    try:
                        if tr.stats.network == 'IU':
                            paz = iu_parser.get_paz(seedid, start_time)
                        elif tr.stats.network == 'S':
                            paz = s_parser.get_paz(seedid, start_time)
                        else:
                            paz = au_parser.get_paz(seedid, start_time)
                            
                    # dodgy station codes assigned in conversion from wfdisc files
                    except:
                        
                        if channel.startswith('EH'):
                            # retry getting info
                            seedid = seedid.replace('EH', 'HH')
                            paz = au_parser.get_paz(seedid, start_time)
                            print 'Resetting channel name...'
                    
                    # remove response - should trigger exception if paz not available
                    tr = tr.simulate(paz_remove=paz)
                    
                    # HNZ responses in dataless seed files are dodgy - trigger exception
                    if channel.startswith('HN') or channel.startswith('BN'):
                        
                        # force these stations to use stationlist.dat
                        tr.stats.station = tr.stats.station.encode('ascii', 'ignore')
                        if tr.stats.station == 'CA3' or tr.stats.station == 'CMC' or tr.stats.station == 'DAM' \
                           or tr.stats.station == 'FAW' or tr.stats.station == 'CMC' or tr.stats.station == 'DAM' \
                           or tr.stats.station == 'GOK' or tr.stats.station.startswith('BK') \
                           or tr.stats.station == 'PIG2' or tr.stats.station == 'PIG' or tr.stats.station == 'EPS' \
                           or tr.stats.station == 'KPK' or tr.stats.station == 'DOW' or tr.stats.station.startswith == 'LM':
                               print 'forceExcept' 
                               forceExcept = forceExcept # ugly, but should force into except loop
                        
                        print '    Ignoring:', seedid
                        addRecDat = False
                        
                    
                    
                    else:
                        addRecDat = True
                                            
                        print 'Using dataless seed:', seedid
                
                # use stationlist file instead
                except:
                                   
                    # fix some channels on the fly
                    if tr.stats.station == 'FORT' or tr.stats.station == 'KMBL' \
                       or tr.stats.station == 'BLDU' or tr.stats.station == 'FORT':
                        tr.stats.channel = tr.stats.channel.encode('ascii','ignore').replace('EH', 'HH')
                        channel = tr.stats.channel
                        
                    seedid = tr.get_id()
                    print 'Using stationlist data:', seedid
                    
                    # get data from stationlist.dat
                    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                          = get_response_info(tr.stats.station, recdate, tr.stats.channel)
                          
                    if not stlo == -12345:
                    
                        # get fft of trace
                        freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
                        
                        # get response for given frequencies
                        real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                                            gain, inst_ty)
                        
                        # deconvolve response                                              
                        corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
                        
                        # make new instrument corrected velocity trace
                        pgv, ivel = get_cor_velocity(corfftr, corffti, freq, inst_ty)
                        
                        # overwrite tr.data
                        tr.data = ivel.real
                        
                        addRecDat = True
                        
                    
                #####################################################################
                # now do ffts
                #####################################################################
                # checks for errors in FFT
                #try:
                if addRecDat == True:
                    # remove post instrument correction low-freq noise
                    tr_filt = tr.copy()
                    tr_filt.filter('highpass', freq=0.05, corners=2, zerophase=True)
                    
                    '''
                    # get noise fft of trace
                    '''
                    f += 1
                    sttr = tr.stats.starttime + tr.stats.npts * tr.stats.delta * 0.02
                    st_offset = tr.stats.npts * tr.stats.delta * 0.02
                    ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta * 0.98 # allow buffer
                    
                    ntr_trim = tr_filt.copy()
                    ntr_trim.trim(sttr, ettr)
                    
                    ntr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp \
                        = get_fft_spectra(ntr_trim, sttr, ettr, channel)
                    
                    traceDat = {'hi_freq_filt': hifreq, 'lo_freq_filt': lofreq, 
                                'noise_spec': smoothed_interp_disp, 'freqs': interp_freqs,
                                'sample_rate': tr.stats.sampling_rate}
                    
                    
                    '''
                    # get s-wave fft of trace
                    '''
                    
                    sttr = tr.stats.starttime + pickDat['sidx'] * tr.stats.delta * 0.98 # allow buffer
                    ettr = tr.stats.starttime + pickDat['rgidx'] * tr.stats.delta * 1.02 # allow buffer
                    
                    str_trim = tr_filt.copy()
                    str_trim.trim(sttr, ettr)
                
                    
                    str_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp \
                        = get_fft_spectra(str_trim, sttr, ettr, channel)
                    
                    traceDat['swave_spec'] = smoothed_interp_disp
                    
                                  
                    '''
                    # get p/s-wave fft of trace
                    '''
                    
                    sttr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta * 0.98 # allow buffer
                    ettr = tr.stats.starttime + pickDat['rgidx'] * tr.stats.delta * 1.02 # allow buffer
                    
                    pstr_trim = tr_filt.copy()
                    pstr_trim.trim(sttr, ettr)
                    
                    pstr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp \
                        = get_fft_spectra(pstr_trim, sttr, ettr, channel)
                    
                    traceDat['p-swave_spec'] = smoothed_interp_disp
                    
                    #####################################################################
                    # add trace data to recDat
                    #####################################################################
                    recDat[channel.encode('ascii', 'ignore').lower()] = traceDat
                    chandict.append(channel.encode('ascii', 'ignore').lower())
                    
                    '''
                    except:
                        print 'Not calculating FFT for:', seedid
                        addRecDat = False
                    '''           
            
            except:
                #print 'Cannot find response info for:', seedid
                addRecDat = False
                
        #####################################################################
        # populate record dictionary
        #####################################################################
        if addRecDat == True:
            recDat['channels'] = chandict
            recDat['sta'] = tr.stats.station.encode('ascii','ignore')
            records.append(recDat)

    else:
        print 'Large file size:', pickDat['mseed_path']
#plt.show()        
#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('wa_fft_data.pkl', 'wb')
pickle.dump(records, pklfile, protocol=-1)
pklfile.close()

# set permissions for all to execute
chmod('wa_fft_data.pkl', 0o777)
        







