from obspy import read, UTCDateTime
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from data_fmt_tools import remove_low_sample_data, return_sta_data, remove_acceleration_data
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay
from io_catalogues import parse_ga_event_query
from os import path, chmod, stat, getcwd
from numpy import arange, sqrt, pi, exp, log, logspace, interp, nan, where, isnan, nanmean
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pickle

def parse_pickfile(pickfile):
    from numpy import nan
    
    # parse pick file
    folder = 'record_picks'
    line = open(path.join(folder, pickfile)).read()
    
    data = line.strip().split(',')
    
    if not data[0].startswith('junk'):
        starttime = UTCDateTime(datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f'))
        origintime = UTCDateTime(datetime.strptime(data[1], '%Y-%m-%dT%H:%M:%S.%f'))
        recdate = datetime.strptime(data[0], '%Y-%m-%dT%H:%M:%S.%f')
        evdate = datetime.strptime(data[1], '%Y-%m-%dT%H:%M:%S.%f')
        
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

def stationlist_fft(tr, pickDat):
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
    
    if inst_ty == 'N':
        dispamp = sqrt((tr.stats.delta)**2 * (corfftr**2 + corffti**2)) / ((2 * pi * freq)**2)
    else:
        dispamp = sqrt((tr.stats.delta)**2 * (corfftr**2 + corffti**2)) / (2 * pi * freq)
    
    staloc = {'latitude':stla, 'longitude':stlo}
    	
    return freq[1:mi], dispamp[1:mi]
    	
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
        
################################################################################
# start here
################################################################################
    
pickfile = '2023-01-05T05.08.AU.ONGER.picks'
pickDat = parse_pickfile(pickfile)

mseedfile = path.split(pickDat['mseed_path'])[-1]
st = read(path.join('iris_dump', mseedfile))

for tr in st:
    # only do Z channel
    if tr.stats.channel.endswith('Z'):
        #traceDat = {}
        seedid = tr.get_id()
        channel = tr.stats.channel
        start_time = tr.stats.starttime 
        
        # demean and taper data  
        tr_proc = tr.copy()
        tr_proc.detrend(type='demean')

pidx = pickDat['pidx']
sidx = pickDat['sidx']
eidx = pickDat['eidx']

################################################################################
# check wave picks
################################################################################

# plot all sata
fig = plt.figure(1, figsize=(10, 4))

plt.plot(tr_proc.times(), tr_proc.data, 'b-', lw=0.5)
plt.plot(tr_proc.times()[sidx:eidx], tr_proc.data[sidx:eidx], 'r-', lw=0.5)
	
plt.show()

################################################################################
# correct response and plt spectra
################################################################################
fig = plt.figure(2, figsize=(8, 8))
	
sttr = tr.stats.starttime + (pickDat['sidx']-5) * tr.stats.delta #- 5. # allow buffer
ettr = tr.stats.starttime + (pickDat['eidx']+5) * tr.stats.delta #+ 5. # allow buffer

str_trim = tr_proc.copy()
str_trim.trim(sttr, ettr)
str_trim.taper(0.02, type='hann', max_length=None, side='both')

# get instrument corrected spectra
freqs, s_disp_amp = stationlist_fft(str_trim, pickDat)

interp_freqs =  logspace(-1,2,241)[32:-24]
smoothed_disp, smoothed_interp_disp = get_smoothed_fft_spectra(freqs, s_disp_amp)

# plot
plt.loglog(freqs, s_disp_amp, 'b-', lw=0.5)
plt.loglog(interp_freqs, smoothed_interp_disp, 'r-', lw=2.5)


# load pickle and compare
recs = pickle.load(open('ongere_fft_data.pkl', 'rb' ))

rec = recs[0]
channel = pickDat['ch3']
plt.loglog(rec[channel]['freqs'], rec[channel]['swave_spec'], 'g-', lw=1.0)

plt.show()