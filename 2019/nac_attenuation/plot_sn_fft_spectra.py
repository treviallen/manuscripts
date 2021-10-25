from obspy import UTCDateTime
from obspy import read #, read_inventory
#from obspy.core.stream import Stream
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from data_fmt_tools import remove_low_sample_data
#from plotting import trim_wave
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay
#from write_data import write_response_spectra
#from mapping_tools import distance
from os import path, chmod, stat, getcwd
from numpy import sqrt, pi, exp, log, logspace, interp, nan, where, isnan, arange
from datetime import datetime
import matplotlib.pyplot as plt
import pickle
import matplotlib as mpl
mpl.style.use('classic')


#plt.ion()

def get_fft_spectra(tr_trim, sttr, ettr, channel):
    
    tr_trim.taper(0.02, type='hann', max_length=None, side='both')
    freq, wavfft = calc_fft(tr_trim.data, tr_trim.stats.sampling_rate)
    mxidx = (len(freq)/2)
    
    mxidx = int(round(mxidx))
    # get fft freqs and amps
    freqs = abs(freq[1:mxidx])
    spec_amps = abs(wavfft.real[1:mxidx])
    
    # get displacement spectra
    if channel[1] == 'N': # if accelerometer
        disp_amps = spec_amps / (2.*(pi)*freqs)**2
    else: # if seismometer
        disp_amps = spec_amps / (2.*(pi)*freqs)
        
    if len(freq) > 10000:
        sw = 201
    elif len(freq) > 5000:
        sw = 101
    elif len(freq) > 1000:
        sw = 51
    elif len(freq) > 500:
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
                               
    return tr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp

################################################################################
# get pick files
################################################################################
#folder = 'test_picks'
folder = 'record_picks'
pickfiles = listdir_extension(folder, 'picks')

################################################################################
# set some defaults
################################################################################

interp_freqs = logspace(-1,2,31)[:-3] # from 0.1-50 Hz

################################################################################
# parse AU dataless
################################################################################
'''
# read dataless seed volumes
print('\nReading dataless seed volumes...')
from obspy.io.xseed import Parser
if getcwd().startswith('/nas'):
    au_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU.IRIS.dataless')
    s1_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/S1/S1.IRIS.dataless')
    ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE1.IRIS.dataless')
    #ge_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/GE/GE2.IRIS.dataless')
    iu_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/IU/IU.IRIS.dataless')
    ii_parser = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/II/II.IRIS.dataless')
    
else:
    au_parser = Parser('/Users/trev/Documents/Networks/AU/AU.IRIS.dataless')
    s1_parser = Parser('/Users/trev/Documents/Networks/S1/S1.IRIS.dataless')
    iu_parser = Parser('/Users/trev/Documents/Networks/IU/IU.IRIS.dataless')
    ge_parser = Parser('/Users/trev/Documents/Networks/GE/GE.IRIS.dataless')
    ii_parser = Parser('/Users/trev/Documents/Networks/II/II.IRIS.dataless')
'''
fig = plt.figure(1, figsize=(8,8))
cmap = plt.cm.get_cmap('Paired',12)
cols = (cmap(arange(12)))
################################################################################
# loop through pick files
################################################################################
f = 0                 
pf  = pickfiles[400]
pf  = pickfiles[800]
pf = '1996-12-30T19.41.AU.KAKA.picks'
pf = '2001-10-07T02.21.AU.KAKA.picks' # used for paper
#pf  = pickfiles[35]
#pf = '1995-12-25T04.43.AU.DRS.picks'
#pf = '1995-12-25T04.43.AU.DPH.picks'

skipRec = False
recDat = {}

# parse pick file
line = open(path.join(folder, pf)).read()

data = line.strip().split(',')

if not data[0].startswith('junk'):
    starttime = UTCDateTime(datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f'))
    origintime = UTCDateTime(datetime.strptime(data[1], '%Y-%m-%dT%H:%M:%S.%f'))
    recdate = datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f')
    
    pickDat = {'starttime': starttime, 'origintime': origintime, \
               'eqlo': float(data[2]), 'eqla': float(data[3]),
               'eqdp': float(data[4]), 'mag': float(data[5]), 'rhyp': float(data[6]),
               'azim': float(data[7]), 'sps': float(data[8]), \
               'ch1': data[9], 'ch2': data[10], 'ch3': data[11], 
               'ppk': float(data[12]), 'spk': float(data[13]), 'epk': float(data[14]), \
               'pidx': int(data[15]), 'sidx': int(data[16]), 'eidx': int(data[17]), 'mseed_path': data[18]}
    
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
        st = read(pickDat['mseed_path'])
    except:
        print('Skipping: '+pickDat['mseed_path'])
        skipRec = True
    
    if skipRec == False:
        
        
        # split trace containing gaps into contiguous unmasked traces
        st = st.split()
        
        # remove low sample rate data
        new_st = remove_low_sample_data(st)
        
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
        
        st_filt = new_st.copy()
        #sidx = int(round(0.05*st_filt[-1].stats.npts))
        #eidx = int(round(0.95*st_filt[-1].stats.npts))
        
        #st_filt.filter('bandpass', freqmin=0.5, freqmax=10, corners=2, zerophase=True)
                
        print('\nReading mseed file:', path.split(pickDat['mseed_path'])[-1])
            
        #####################################################################
        # loop thru traces
        #####################################################################
        chandict = []
        for tr in new_st:
            traceDat = {}
            seedid = tr.get_id()
            channel = tr.stats.channel
            start_time = tr.stats.starttime 
            
            # demean and taper data
            tr = tr.detrend(type='demean')
            #tr = tr.taper(0.02, type='hann', max_length=None, side='both')
            
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
            xdat = range(0, tr.stats.npts)
            #plt.plot(xdat, tr.data, 'b-', lw=0.5)
            pidx = pickDat['pidx']
            eidx = pickDat['eidx']
            #####################################################################
            # now do ffts
            #####################################################################
            # checks for errors in FFT
            #try:
            # remove post instrument correction low-freq noise
            tr_filt = tr.copy()
            tr_filt.filter('highpass', freq=lofreq, corners=2, zerophase=True)
            
            '''
            # get noise fft of trace
            '''
            f += 1
            sttr = tr.stats.starttime + 10. # should add 1-2 secs instead of proportion 
            ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 10. # allow buffer
            
            if sttr > ettr:
                sttr = tr.stats.starttime + 2. # should add 1-2 secs instead of proportion 
                ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 2. # allow buffer
                
            if sttr > ettr:
                sttr = tr.stats.starttime + 1. # should add 1-2 secs instead of proportion 
                ettr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 1. # allow buffer
            
            ntr_trim = tr_filt.copy()
            ntr_trim.trim(sttr, ettr)
            
            ntr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp \
                = get_fft_spectra(ntr_trim, sttr, ettr, channel)
            
            if channel.endswith('Z'):
                plt.loglog(freqs, disp_amps, '-', lw=0.5, c=cols[0], label='Noise (Raw)')
                plt.loglog(freqs, smoothed_disp, '-', lw=1.5, c=cols[1], zorder=3, label='Noise (Smoothed)')

            
            traceDat = {'hi_freq_filt': hifreq, 'lo_freq_filt': lofreq, 
                        'noise_spec': smoothed_interp_disp, 'freqs': interp_freqs,
                        'sample_rate': tr.stats.sampling_rate}
                
            #plt.plot(xdat[0:eidx], tr.data[0:eidx], 'g-', lw=0.5)
                          
            '''
            # get p/s-wave fft of trace
            '''
            
            sttr = tr.stats.starttime + pickDat['pidx'] * tr.stats.delta - 10. # allow buffer
            ettr = tr.stats.starttime + pickDat['eidx'] * tr.stats.delta + 10. # allow buffer
            
            pstr_trim = tr_filt.copy()
            pstr_trim.trim(sttr, ettr)
            
            pstr_trim, freqs, disp_amps, smoothed_disp, smoothed_interp_disp \
                = get_fft_spectra(pstr_trim, sttr, ettr, channel)
            
            traceDat['p-swave_spec'] = smoothed_interp_disp
            
            if channel.endswith('Z'):
                plt.loglog(freqs, disp_amps, '-', lw=0.5, c=cols[4], label='Signal (Raw)')
                plt.loglog(freqs, smoothed_disp, '-', lw=1.5, c=cols[5], label='Signal (Smoothed)')
            
            '''
            # get SN-Ratio
            '''
            traceDat['sn_ratio'] = traceDat['p-swave_spec'] / traceDat['noise_spec']
            
            # now set frequency limits - use 1 Hz as centre
            sn_thresh = 4.
            
            # find nan ratios
            nanidx = where(isnan(traceDat['sn_ratio']))[0]
            sn_ratio = traceDat['sn_ratio']
            sn_ratio[nanidx] = 0.
            
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
            
            # do a manual check for good data with no noise
            if sn_ratio[15] > 100.:
                traceDat['lof_limit'] = interp_freqs[0]
                fidx = where((interp_freqs >= 3) & (sn_ratio < sn_thresh))[0]
                if len(fidx) == 0:
                    traceDat['hif_limit'] = interp_freqs[-1]
                else:
                    traceDat['hif_limit'] = interp_freqs[fidx[0]-1]
                    
            # now check if short period and set lof_limit to 0.4
            if channel.startswith('E') or channel.startswith('S'):
                if traceDat['lof_limit'] < 0.4:
                    traceDat['lof_limit'] = 0.4
            
            if channel.endswith('Z'):
                ll = traceDat['lof_limit']
                hl = traceDat['hif_limit']
                plt.fill([hl, hl, ll, ll, hl], [1E-2, 1E8, 1E8, 1E-2, 1E-2], 'w', ec='0.2', label='Used Frequencies', alpha=0.25)
                plt.fill([0.01, 0.01, ll, ll, 0.01], [1E-2, 1E8, 1E8, 1E-2, 1E-2], '0.8', ec='0.8', label='Rejected Frequencies')
                plt.fill([100, 100, hl, hl, 100], [1E-2, 1E8, 1E8, 1E-2, 1E-2], '0.8', ec='0.8')
            
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
        recDat['ev'] = data[1][0:16].replace(':','.')
        recDat['mseed_path'] = pickDat['mseed_path']
        
        
        '''
        else:
            print('Large file size:', pickDat['mseed_path'])
        '''
plt.legend(loc=1, fontsize=13)
plt.xlim([0.07, 30])
plt.ylim([0.1, 2E6])
plt.grid(which='both')
plt.xlabel('Frequency (Hz)', fontsize=16)
plt.ylabel('Fourier Displacement (Count-s)', fontsize=16)
plt.savefig('sn_limits_example.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig_5.eps', fmt='eps', dpi=300, bbox_inches='tight')
plt.show()
      







