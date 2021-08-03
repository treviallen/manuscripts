#from process_waves import common_resp, common_resp
from response import get_response_info, paz_response, deconvolve_instrument
from spectral_analysis import calc_fft, get_cor_velocity
from obspy import read
from numpy import arange, fft
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.constants import g
mpl.style.use('classic')

##############################################################################
# start main    
##############################################################################

# set file
mseedfile = '/Users/trev/Documents/Earthquake_Data/20180916.Lake_Muir/2018-09-16T04.50.00.AU.RKGY.mseed'
mseedfile = '2018-11-08T21.05.AU.LM04.mseed'
plt_channels = ['HNE', 'HNN', 'HNZ']
#plt_channels = ['HHE', 'HHN', 'HHZ']

nwid = len(plt_channels)

fig = plt.figure(1, figsize=(18,10))

# read mseed
st = read(mseedfile)

# loop thru traces
for tr in st:
    for i, channel in enumerate(plt_channels):
        if tr.stats.channel == channel:
            c = i+1
            seedid = tr.get_id()
            #channel = tr.stats.channel
            start_time = tr.stats.starttime
            
            print(seedid, channel)
            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                  = get_response_info(tr.stats.station, tr.stats.starttime.datetime, tr.stats.channel)

            # get fft of trace
            freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
            
            # get response for given frequencies
            real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                                gain, inst_ty)
            
            # deconvolve response
            corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
            
            # make new instrument corrected velocity trace
            complex_array = corfftr + 1j*corffti
            n = len(corfftr)
            iwave = fft.ifft(complex_array,n)
            tr.data = iwave.real
                
            if channel[1] == 'H':
                tr_vel = tr.copy()
                
                tr_vel.filter('bandpass', freqmin=0.2, freqmax=50., corners=2, zerophase=True)
                
                # trim trace
                tr_vel.trim(starttime=tr.stats.starttime+360,endtime=tr.stats.starttime+460)
                
                tr_disp = tr_vel.copy()
                tr_acc = tr_vel.copy()
                
                tr_disp.integrate()
                tr_acc.differentiate()
                
            elif channel[1] == 'N':
                tr_acc = tr.copy()
                
                tr_acc.filter('bandpass', freqmin=0.4, freqmax=50., corners=2, zerophase=True)
                
                tr_acc.trim(starttime=tr.stats.starttime+118,endtime=tr.stats.starttime+152)
                
                tr_disp = tr_acc.copy()
                tr_vel = tr_acc.copy()
                
                tr_disp.integrate().integrate()
                tr_vel.integrate()
                
            # now make plots
            times = tr_vel.times()
            starttime = tr_vel.stats.starttime.datetime
            
            dt_times = []
            for time in times:
                dt_times.append(starttime + timedelta(seconds=time))
                
            # make acc plot
            ax = plt.subplot(nwid, 3, c)
            plt.plot(dt_times, (100*tr_acc.data/g), 'r-', lw=0.5, label=seedid)
            
            ticks = ax.get_xticks()[::2]
            ax.set_xticks(ticks)
            labels = [mpl.dates.num2date(t).strftime('%H:%M:%S') for t in ticks]
            ax.set_xticklabels(labels)
            plt.ylim([-13, 13])
            
            '''
            idx=list(arange(1,10,2))
            for label in ax.xaxis.get_ticklabels()[idx]:
                label.set_visible(False)
            '''
            if c == 1:
                plt.ylabel('Acceleration (%g)', fontsize=16)
                plt.title('East-West', fontsize=18)
            elif c == 2:
                plt.title('North-South', fontsize=18)
            elif c == 3:
                plt.title('Vertical', fontsize=18)
            
            # make vel plot
            ax = plt.subplot(nwid, 3, nwid+c)
            plt.plot(dt_times, tr_vel.data*1000., 'r-', lw=0.5, label=seedid)
            
            ax.set_xticks(ticks)
            labels = [mpl.dates.num2date(t).strftime('%H:%M:%S') for t in ticks]
            ax.set_xticklabels(labels)
            plt.ylim([-22, 22])
            
            if c == 1:
                plt.ylabel('Velocity (mm/s)', fontsize=16)
            
            # make vel plot
            ax = plt.subplot(nwid, 3, 2*nwid+c)
            plt.plot(dt_times, tr_disp.data*1000., 'r-', lw=0.5, label=seedid)
            
            ax.set_xticks(ticks)
            labels = [mpl.dates.num2date(t).strftime('%H:%M:%S') for t in ticks]
            ax.set_xticklabels(labels)
            plt.ylim([-2, 3])
            
            if c == 1:
                plt.ylabel('Displacement (mm)', fontsize=16)
            
            plt.xlabel('UTC Time', fontsize=16)
            
plt.savefig('20181108.lake_muir_LM04.png', fmt='png', dpi=150, bbox_inches='tight')
plt.show()

