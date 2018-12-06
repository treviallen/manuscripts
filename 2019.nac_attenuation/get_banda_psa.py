from obspy import UTCDateTime
from obspy import read #, read_inventory
from obspy.core.stream import Stream
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from plotting import trim_wave
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension
from write_data import write_response_spectra
from mapping_tools import distance
from os import path
from numpy import sqrt
from datetime import datetime

def parse_usgs_events(usgscsv):
    lines = open(usgscsv).readlines()[1:]
    
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'time': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict

################################################################################
# loop through earthquakes and get data
################################################################################
usgscsv = 'neic_banda_png_extra.csv'
usgscsv = 'neic_banda_png_5-6.csv'
evdict = parse_usgs_events(usgscsv)

# read dataless seed volumes
from obspy.io.xseed import Parser
#au_parser = Parser('AU.dataless')
au_parser = Parser(path.join('..','..','..','Networks','AU', 'AU.dataless'))
s_parser = Parser(path.join('..','..','..','Networks','S', 'S.dataless'))
ge_parser = Parser(path.join('..','..','..','Networks','GE', 'GE.dataless'))
ms_parser = Parser(path.join('..','..','..','Networks','MS', 'MS.dataless'))
my_parser = Parser(path.join('..','..','..','Networks','MY', 'MY.dataless'))
tm_parser = Parser(path.join('..','..','..','Networks','TM', 'TM.dataless'))

                    
from obspy.clients.arclink.client import Client
arclink_client = Client(user='trevor.allen@ga.gov.au')

from obspy.clients.fdsn.client import Client
iris_client = Client("IRIS")

# get mseed datafiles
extension = 'mseed'
folder = 'mseed_dump_new'
#folder = 'mseed_jump'
mseedfiles = listdir_extension(folder, extension)

sstxt = ''

for ev in evdict:
    for msf in mseedfiles:
        newst = Stream()
        filename = path.split(msf)[-1]
        if filename.startswith(ev['timestr'][:16].replace(':', '.')):
            print filename
            st = read(path.join(folder, msf))
            st = st.merge(method=0, fill_value='interpolate')
            
            # check to see if we delete traces
            if len(st) > 3:
                hh = False
                eh = False
                for tr in st:
                    if tr.stats.channel.startswith('HH'):
                        hh = True
                    elif tr.stats.channel.startswith('EH'):
                        eh = True
                        
                # now delete traces
                if hh == True:
                    for tr in st:
                        if tr.stats.channel.startswith('HH') == False:
                            print 'Deleting channel', tr.stats.station, tr.stats.channel
                            st.remove(tr)
                
                elif eh == True:
                    for tr in st:
                        if tr.stats.channel.startswith('EH') == False:
                            print 'Deleting channel', tr.stats.station, tr.stats.channel
                            st.remove(tr)
            
            #####################################################################
            # now start grunt work
            #####################################################################
            try:
                for tr in st:
                    
                    if tr.stats.channel.startswith('BN') == False: # ignore accelerometers for now
                        '''
                        # if using GEOFON data, need to request response
                        if tr.stats.network == 'GE':
                            paz = arclink_client.get_paz(tr.stats.network, tr.stats.station, '', \
                                                         tr.stats.channel, ev['starttime'])
                                                         
                            meta = arclink_client.get_metadata(tr.stats.network, tr.stats.station, '', \
                                                               tr.stats.channel, ev['starttime'])
                            staloc = meta['coordinates']
                        
                        else:
                        '''
                        #####################################################################
                        # get instrument response
                        #####################################################################
                        
                        seedid=tr.get_id()
                        channel = tr.stats.channel
                        start_time = tr.stats.starttime
                        
                        print filename, tr.stats.channel
                        
                        if tr.stats.network == 'AU':
                            paz = au_parser.get_paz(seedid,start_time)
                            staloc = au_parser.get_coordinates(seedid,start_time)
                        elif tr.stats.network == 'GE':
                            paz = ge_parser.get_paz(seedid,start_time)
                            staloc = ge_parser.get_coordinates(seedid,start_time)
                        elif tr.stats.network == 'MY':
                            paz = my_parser.get_paz(seedid,start_time)
                            staloc = my_parser.get_coordinates(seedid,start_time)
                        elif tr.stats.network == 'MS':
                            paz = ms_parser.get_paz(seedid,start_time)
                            staloc = ms_parser.get_coordinates(seedid,start_time)
                        elif tr.stats.network == 'TM':
                            paz = tm_parser.get_paz(seedid,start_time)
                            staloc = tm_parser.get_coordinates(seedid,start_time)
                        else:
                            paz = s_parser.get_paz(seedid,start_time)
                            staloc = s_parser.get_coordinates(seedid,start_time)
                        
                        #####################################################################
                        # remove instrument response
                        #####################################################################
                        
                        tr = tr.detrend(type='demean')
                        tr = tr.taper(0.05, type='hann', max_length=None, side='both')
                        
                        # EN? dodgy stn channel from parsing JUMP data
                        if tr.stats.channel.startswith('BH') or tr.stats.channel.startswith('HH') \
                           or tr.stats.channel.startswith('EN'):
                            lofreq = 0.075
                        else:
                            lofreq = 0.2
                        
                        hifreq = min([12., tr.stats.sampling_rate/2.])
                        
                        tr = tr.filter('bandpass', freqmin=lofreq, freqmax=hifreq, \
                                        corners=2, zerophase=True)
                        
                        # check if jump site and if so, use stationlist data - IRIS dataless does not work!
                        if tr.stats.station == 'DPH' or tr.stats.station == 'CN1H' \
                           or tr.stats.station == 'DRS' or tr.stats.station == 'CN2S' \
                           or tr.stats.station == 'AS31':
                            
                            recdate = datetime.strptime(ev['timestr'], "%Y-%m-%dT%H:%M:%S.%fZ")
                            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                                  = get_response_info(tr.stats.station, recdate, tr.stats.channel)
                            
                            # get fft of trace
                            freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
                            
                            # get response for given frequencies
                            real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                                                gain, inst_ty)
                            
                            # deconvolve response                                              
                            corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
                            
                            # make new instrument corrected velocity trace
                            pgv, ivel = get_cor_velocity(corfftr, corffti, freq, inst_ty)
                            
                            tr.data = ivel.real
                            
                        # just use IRIS dataless volume - much easier, but less transparent!
                        else:
                            tr = tr.simulate(paz_remove=paz)
                            
                            
                        '''
                        # else if using FDSN data
                        else:
                            pre_filt = [0.001, 0.005, 45, 50]
                            tr.remove_response(pre_filt=pre_filt, taper=True, taper_fraction=0.05, output="VEL", zero_mean=True)
                        '''
                        #####################################################################
                        # get response spectra
                        #####################################################################
                        
                        # prep for response spectra
                        freq, wavfft = calc_fft(tr.data, tr.stats.sampling_rate)
                        
                        # prep for psa
                        iacc = prep_psa_simple(wavfft.real, wavfft.imag, freq, 'B')
                        
                        # calc response spectra
                        h = 5.0 # %
                        minT = 0.1
                        maxT = 10
                        T, psa, pga = calc_response_spectra(iacc, tr.stats.sampling_rate, h, minT, maxT)
                        pgv = max(abs(tr.data))
                        
                        repi = distance(ev['lat'], ev['lon'], staloc['latitude'], staloc['longitude'])[0]
                        rhyp = sqrt(repi**2 + ev['dep'])
                        azim = distance(ev['lat'], ev['lon'], staloc['latitude'], staloc['longitude'])[1]
                        
                        psafilename = msf[:-6] + '.' + tr.stats.channel
                        
                        write_response_spectra(tr.stats.station, ev['time'], tr.stats.sampling_rate, \
                                               T, psa, pga, pgv, psafilename, staloc['latitude'], staloc['longitude'], \
                                               ev['lat'], ev['lon'], ev['dep'], ev['mag'], rhyp, azim, lofreq, hifreq)
                        
                        newst += tr
                        
                # trim data
                #start_index, stop_index = trim_wave(newst[0].data, tr.stats.sampling_rate, 'S', False)
            
                # save corrected waves
                newfname = '.'.join((msf[:-6], 'cor', 'mseed'))
                newfpath = path.join('mseed_cor', newfname)
                newst.write(newfpath, format='MSEED')
                
                # add to start/stop list
                #sstxt += ','.join((newfpath, str(start_index), str(stop_index))) + '\n'
            except:
                print 'Failed:', filename

