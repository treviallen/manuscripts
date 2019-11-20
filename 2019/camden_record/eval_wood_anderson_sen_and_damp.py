# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:33:23 2019

@author: u56903
"""

import response, readwaves, spectral_analysis, plotting, write_data, spatial_tools
import calculate_magnitudes
import numpy as np
from sys import argv
from os import getcwd, path
from datetime import datetime

wavfile = argv[1]

def common_read(allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, sacseed):
    
    #print(len(alldata[:,17])
    # select component
    chan, chan_no = readwaves.select_channel(allsta, comps)

    # unify kelunji channel names
#    tmpchan = chan.split(' Trans')
#    chan = ''
#    for i in range(0,len(tmpchan)):
#        chan = chan + tmpchan[i]

    # get channel data of interest
    sta = allsta[chan_no].strip()
    sps = int(allsps[chan_no])
    
    if sacseed == True:
        chan_dat = alldata[:,chan_no]
    else:
        chan_dat = alldata[chan_no] # edited as of 2017-09-26 to read new eqwave
    
    tmprecdate = allrecdate[chan_no]
    recdate = datetime.strptime(tmprecdate, "%Y%m%d%H%M%S")

    # remove trailing nans
    nantrue = True
    while nantrue == True:
        if chan_dat[-1] == np.nan:
            chan_dat = chan_dat[0:-1]
        else:
            nantrue = False

    chan_dat = chan_dat[0:allnsamp[chan_no]]
    chan_dat = chan_dat.reshape(1,len(chan_dat))

    # remove DC offset
    chan_dat = readwaves.remove_dc_offset(chan_dat)

    # check to see if response exists
    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
        = response.get_response_info(sta,recdate,chan)

    # if info does not exist, ask for user input
    if nat_freq == -12345 and pazfile == 'NULL':
        inst_ty, mindate, maxdate, stlo, stla, netid, nat_freq, damping, \
            sen, recsen, gain, pazfile = response.enter_response_info(sta, chan, sps)

        # write new data to file
        response.write_response_info(sta, inst_ty, mindate, maxdate, stlo, \
            stla, netid, nat_freq, damping, sen, recsen, gain, chan, pazfile)
    
    return sta, inst_ty, sps, recdate, nat_freq, damping, sen, recsen, gain, \
           chan, chan_no, chan_dat, stlo, stla, pazfile, alldata, netid

"""****************************************************************************
Common FFT functions
****************************************************************************"""
def common_fft(chan_dat, inst_ty, sps, seltask):
    # view wave and trim if necessary
    print(chan_dat)
    start_index, stop_index = plotting.trim_wave(chan_dat, sps, inst_ty, False)
    #taper_dat = chan_dat[0, start_index:stop_index]
    taper_dat = chan_dat[start_index:stop_index]
    taper_dat = taper_dat.reshape(1, len(taper_dat))

    # taper time-domain before fft
    taper_dat, costaper = spectral_analysis.cosine_taper(taper_dat)
    
    # do acausal filter in time domain using 4th order Butterworth filter
#    lofreq, hifreq, filtwave = spectral_analysis.acausal_butter_filter(inst_ty, sps, taper_dat, seltask)
#    filtwave = taper_dat.reshape(1, len(filtwave))

    # get FFT of data stream for full record
#    freq, wavfft = spectral_analysis.calc_fft(filtwave, sps)
    freq, wavfft = spectral_analysis.calc_fft(taper_dat[0], sps)

    # filter FFT using 4th order Butterworth filter
    lofreq, hifreq, wavfft = spectral_analysis.butter_filter_user(inst_ty, sps, freq, wavfft, seltask)

    # get time offset for filenames
    dt = start_index / float(sps)

    return freq, lofreq, hifreq, wavfft, dt

 #   return freq, wavfft

"""****************************************************************************
Common functions to remove FAP response
****************************************************************************"""

def common_resp(freq, nat_freq, damping, sen, recsen, gain, wavfft, inst_ty):

    # calculate FAP instrument response
    if pazfile == 'NULL':
        real_resp, imag_resp = response.fap_response(freq, nat_freq, damping, \
                                                     sen, recsen, gain, inst_ty)

    # or use PAZ file
    else:
        real_resp, imag_resp = response.paz_response(freq, pazfile, sen, recsen, \
                                                     gain, inst_ty)

    # deconvolve instrument response from record and return corrected FFT
    corfftr, corffti = response.deconvolve_instrument(real_resp, imag_resp, wavfft)

    return corfftr, corffti, real_resp, imag_resp

"""****************************************************************************
Start main code
****************************************************************************"""
seltask = '4'
try:
    # try importing obspy modules
    from obspy.core import read

    print('\nObsPy Installed')

    # try testing if SAC or miniSEED
    try:
        st = read(wavfile)
        sacseed = True

    except:
        sacseed = False

except:
    # print(statement if cannot access obspy
    print('\nCannot import ObsPy modules!')

# if SAC or SEED read waves
if sacseed == True:
    allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readseed(st)
    fmt = 'obspy'

# esle read text file format
else:
    # check text file format
    fmt = readwaves.check_file_fmt(wavfile)
    if fmt == 'eqw':
        # read the eqWave text file
        try:
            allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, cntpvolt, allsen, allgain = readwaves.readeqwave(wavfile)
        except:
            allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readeqwave(wavfile)
    elif fmt == 'nmx':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readnmx(wavfile)
    elif fmt == 'tspair':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readtspair(wavfile)
    elif fmt == 'sm':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readseismac(wavfile)
    else:
        '\nFile format not recognised!'

# do common read functions
sta, inst_ty, sps, recdate, nat_freq, damping, sen, recsen, gain, chan, \
        chan_no, chan_dat, stlo, stla, pazfile, alldata, netid = \
        common_read(allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, sacseed)

# do common fft functions
chan_dat = chan_dat[0]
freq, lofreq, hifreq, wavfft, dt = common_fft(chan_dat, inst_ty, sps, seltask)

# do common instrument deconvolution functions
corfftr, corffti, real_resp, imag_resp = common_resp(freq, nat_freq, damping, sen, \
                                         recsen, gain, wavfft, inst_ty)

# get filename
filename, evdate = write_data.get_filename(sta, recdate, chan, inst_ty, dt, sps)

# get source-to-site distance
rhyp, azim, eqla, eqlo, eqdep, eqmag, evdate = spatial_tools.get_eq_distance(stlo, stla, evdate)
    
# loop through magnification and damping
waMagnification = [2080, 2800]
damping = [0.7, 0.8]

wafile = 'eval_wood-anderson_params.csv'
for wam in waMagnification:
    for damp in damping:
        # convolve with Wood-Anderson instrument and get displacement wave
        wadisp = response.convolve_WoodAnderson(freq, corfftr, corffti, inst_ty, ampfact=wam, damping=damp)
        
        # plot W-A displacement time history and get windo for ML
        watrim = plotting.plot_WoodAnderson(wadisp, sps, filename, chan_no)
        
        # calculate magnitudes
        logA = np.log10(max(abs(watrim)))
        mlDict = calculate_magnitudes.main(filename, logA, rhyp, eqdep)
        
        # write to file
        lines = open(wafile).read()
    
        # append new line                      
        newline = ','.join((sta, str(wam), str(damp), str('%0.1f' % mlDict['rhyp']), \
                            str('%0.1f' % mlDict['repi']), str('%0.4f' % mlDict['MLM92'])))
        
        f = open(wafile, 'w')
        f.write(lines + '\n' + newline)
        f.close()
        





















