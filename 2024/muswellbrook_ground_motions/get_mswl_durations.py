# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 17:02:41 2019

@author: u56903
"""
def use_stationlist_fft(tr, measure='acc'):
    from numpy import fft, pi
    
    seedid = tr.get_id()
    channel = tr.stats.channel
    start_time = tr.stats.starttime
    
    recdate = str(start_time.year)+str(start_time.month).zfill(2)+str(start_time.day).zfill(2)
    
    #print(tr.stats.station, recdate.datetime, tr.stats.channel, tr.stats.network)        
    
    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
          = get_response_info(tr.stats.station, start_time.datetime, tr.stats.channel, tr.stats.network)
    
    # get fft of trace
    freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
    mi = (len(freq)/2)
    mi = int(round(mi))
    
    # get response for given frequencies
    real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                        gain, inst_ty)
    
    # deconvolve response
    corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
    n = len(corfftr)
    
    complex_array = corfftr + 1j*corffti
    
    
    '''
    if measure == 'vel': # why does this not work?!?!
        complex_array *= (2 * pi * freq)
        #complex_array[0] = 0. + 0j
        
    elif measure == 'disp':
        complex_array *= (2 * pi * freq)**2
        complex_array[0] = 0. + 0j
    '''
    itr = fft.ifft(complex_array,n)
    
    staloc = {'latitude':stla, 'longitude':stlo}
    	
    return itr

def set_ylims(plt, ax):
    from numpy import array
    ylims = ax.get_ylim()
    
    absmax = max(abs(array(ylims)))
    
    plt.ylim([-1*absmax, absmax])


# a function to annotate max and min ground motions
def annotate_maxmin(plt, ax ,tr):
    import numpy as np
    
    times = tr.times()
    starttime = tr.stats.starttime.datetime

    maxval = str("%0.3f" % max(tr.data))
    minval = str("%0.3f" % min(tr.data))
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    ymin = -1*max(np.abs(ylim))
    ymax = max(np.abs(ylim))
    ax.set_ylim((ymin, ymax))
    plt.annotate('max = '+maxval, (0,0), fontsize=10, \
                 xytext=((xlim[1]-xlim[0])*0.02, ymax*0.80), xycoords='data')
    plt.annotate('min = '+minval, (0,0), fontsize=10, \
                 xytext=((xlim[1]-xlim[0])*0.02, ymin*0.85), xycoords='data')


# copied from https://github.com/eng-tools/eqsig/blob/master/eqsig/im.py

def calc_sig_dur_vals(motion, dt, start=0.05, end=0.95, se=False):
    import numpy as np
    """
    Computes the significant duration using cumulative acceleration according to Trifunac and Brady (1975).

    Parameters
    ----------
    motion: array-like
        acceleration time series
    dt: float
        time step
    start: float, default=0.05
        threshold to start the duration
    end: float, default=0.95
        threshold to end the duration
    se: bool, default=False
        If true then return the start and end times

    Returns
    -------
    tuple (start_time, end_time)
    """

    cum_acc2 = np.cumsum(motion ** 2)
    ind2 = np.where((cum_acc2 > start * cum_acc2[-1]) & (cum_acc2 < end * cum_acc2[-1]))
    start_time = ind2[0][0] * dt
    end_time = ind2[0][-1] * dt

    if se:
        return start_time, end_time
    return end_time - start_time
                 
#from obspy.io.xseed import Parser
from data_fmt_tools import remove_low_sample_data
from obspy import read, UTCDateTime
from scipy.constants import g
from os import getcwd
from sys import argv
from datetime import timedelta
from response import get_response_info, paz_response, deconvolve_instrument
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
import warnings
warnings.filterwarnings("ignore")

# read file
mseed = '/Users/trev/Documents/Earthquake_Data/20240906.Muswellbrook/2024-09-06T19-58-00.YW.MSWL1.ga2024rqpnyt.mseed'
mseed = '/Users/trev/Documents/Earthquake_Data/20240906.Muswellbrook/2024-11-12T01-12-56.YW.MSWL6.ga2024widwze.mseed'

mseed = argv[1]    

st = read(mseed)
# remove low sample rate data
#new_st = remove_low_sample_data(st)
#st.merge()

print('\n'+str(st[0].stats.starttime))

for tr in st:
    if tr.stats.channel.startswith('HN'):

        seedid = tr.get_id()
        channel = tr.stats.channel
        start_time = tr.stats.starttime
        
        start_off = 60
        end_off = 80
        #end_off = 200
        
        #paz = iu_parser.get_paz(seedid, start_time)
        #paz = au_parser.get_response(seedid,start_time)
        #staloc = au_parser.get_coordinates(seedid,start_time)
        
        # demean and taper data
        tr = tr.detrend(type='demean')
        tr = tr.taper(0.02, type='hann', max_length=None, side='both')
        #print(len(tr.data))
        tr.data = use_stationlist_fft(tr, measure='acc')
        
        tr.filter("bandpass", freqmin=0.15, freqmax=100.) 
        
        durn = calc_sig_dur_vals(tr.data, tr.stats.delta, start=0.05, end=0.95, se=False)
        
        printtxt = ','.join((tr.stats.station, tr.stats.channel, str('%0.2f' % durn)))
        print(printtxt)
        
    elif len(st) == 1:
        tr = st[0].copy()
        
        seedid = tr.get_id()
        channel = tr.stats.channel
        start_time = tr.stats.starttime
        
        #crash
        # demean and taper data
        tr = tr.detrend(type='demean')
        tr = tr.taper(0.02, type='hann', max_length=None, side='both')
        
        tr.filter("bandpass", freqmin=0.2, freqmax=100.) 
        
        #print(len(tr.data))
        tr.data = use_stationlist_fft(tr, measure='vel')
        tr.differentiate()
        
        durn = calc_sig_dur_vals(tr.data, tr.stats.delta, start=0.05, end=0.95, se=False)
        
        printtxt = ','.join((tr.stats.station, tr.stats.channel, str('%0.2f' % durn)))
        print(printtxt)