from obspy import UTCDateTime
from obspy import read #, read_inventory
from obspy.core.stream import Stream
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
#from plotting import trim_wave
from response import get_response_info, paz_response, deconvolve_instrument
from data_fmt_tools import remove_low_sample_data
from misc_tools import listdir_extension
from write_data import write_response_spectra
from mapping_tools import distance
from os import path, getcwd
from numpy import sqrt
from datetime import datetime
import matplotlib.pyplot as plt
plt.ion()

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
# get pick files
################################################################################
#folder = 'test_picks'
folder = 'record_picks'
folder = 'record_picks_new'
pickfiles = listdir_extension(folder, 'picks')

################################################################################
# loop through earthquakes and get data
################################################################################
usgscsv = '20190625_merged_events.csv'
evdict = parse_usgs_events(usgscsv)

# read dataless seed volumes
from obspy.io.xseed import Parser
#au_parser = Parser('AU.dataless')
print('Reading dataless seed volumes...')
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
    ge_parser = Parser('/Users/trev/Documents/Networks/GE/GE1.IRIS.dataless')
    #ge2_parser = Parser('/Users/trev/Documents/Networks/GE/GE2.IRIS.dataless')

                    
from obspy.clients.arclink.client import Client
arclink_client = Client(user='trevor.allen@ga.gov.au')

from obspy.clients.fdsn.client import Client
iris_client = Client("IRIS")

'''
# get mseed datafiles
extension = 'mseed'
folder = 'record_picks'
#folder = 'mseed_dump'
mseedfiles = listdir_extension(folder, extension)
'''
sstxt = ''

# loop thru pickfiles
#for p, pf in enumerate(pickfiles[0:2900][::-1]):
    #for p, pf in enumerate(pickfiles[3000:]):
for p, pf in enumerate(pickfiles):
    # parse pick file
    line = open(path.join(folder, pf)).read()
    
    if not line.startswith('junk'):
    
        data = line.strip().split(',')
        
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
    
        channels = [pickDat['ch1'], pickDat['ch2'], pickDat['ch3']]
        
        #newst = Stream()
        msf = path.join('mseed_dump', path.split(pickDat['mseed_path'])[-1])
                
        st = read(msf)
        st = st.merge(method=0, fill_value='interpolate')
        
        # remove junk channels
        st = remove_low_sample_data(st)
        
        # trim to picks to +/1 10 sec
        #taperTime = len(st[0].data)*0.05 / pickDat['sps']
        #if taperTime > 10:
        staperTime = 60.
        etaperTime = 300.
    
        #print taperTime
        if pickDat['origintime'] + pickDat['ppk'] - staperTime > st[0].stats.starttime:
            starttime = pickDat['origintime'] + pickDat['ppk'] - staperTime
        else:
            starttime = st[0].stats.starttime
        
        if pickDat['origintime'] + pickDat['epk'] + etaperTime > st[0].stats.endtime:
            endtime = st[0].stats.endtime
        else:
            endtime = pickDat['origintime'] + pickDat['epk'] + etaperTime
        
        '''
        # for plotting
        fig = plt.figure(p+1, figsize=(12, 4))
        tr = st[0]
        reftime = tr.stats.starttime-origintime
        plt.plot(reftime+tr.times(), tr.data, 'b-', lw=0.5)
        ptime = pickDat['ppk']
        etime = pickDat['epk']
        '''
        
        
        '''
        # plot to test
        trt = st_trim[0]
        treftime = trt.stats.starttime-origintime    
        plt.plot(treftime+trt.times(), trt.data, 'r-', lw=0.5)
        plt.plot([ptime, ptime], [max(tr.data), min(tr.data)], 'k--')
        plt.plot([etime, etime], [max(tr.data), min(tr.data)], 'k--')
        plt.title(pf)
        plt.show()
        '''
        
        #####################################################################
        # now start grunt work
        #####################################################################
        try:
            # now trim
            st_trim = st.copy()    
            st_trim = st.trim(starttime, endtime)
        
            for tr in st_trim:
                doChan = False
                
                
                psa_path = path.join('psa', path.split(msf[:-6])[-1] + '.' + tr.stats.channel + '.psa')
                    
                # first check that we want to do this channel
                for channel in channels:
                    if tr.stats.channel == channel:
                        doChan = True
                        print(psa_path)
                        if channel.endswith('Z'):
                            doChan = False
                    
                    # chack if psa file exists and set to false
                    if path.isfile(psa_path):
                        doChan = False
                        print('PSA file exists '+str(p))
                    
                    if doChan == True:
                        #if tr.stats.channel.startswith('BN') == False: # ignore accelerometers for now
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
                        #print(seedid)
                        
                        #####################################################################
                        # remove instrument response
                        #####################################################################
                        
                        tr = tr.detrend(type='demean')
                        tr = tr.taper(0.05, type='hann', max_length=30., side='both')
                        
                        # EN? dodgy stn channel from parsing JUMP data
                        if tr.stats.channel.startswith('BH') or tr.stats.channel.startswith('HH') \
                           or tr.stats.channel.startswith('EN'):
                            lofreq = 0.075
                        else:
                            lofreq = 0.2
                        
                        hifreq = min([12., 0.45*tr.stats.sampling_rate])
                        
                        tr = tr.filter('bandpass', freqmin=lofreq, freqmax=hifreq, \
                                        corners=2, zerophase=True)
                        
                        # check if jump site and if so, use stationlist data - IRIS dataless does not work!
                        #tr = tr.simulate(paz_remove=paz)
                        print(tr.stats.station)
                        if tr.stats.station == 'AS31' or tr.stats.station.startswith('PH0'):
                            
                            #recdate = datetime.strptime(ev['timestr'], "%Y-%m-%dT%H:%M:%S.%fZ")
                            recdate = pickDat['origintime']
                            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                                  = get_response_info(tr.stats.station, recdate, tr.stats.channel)
                            
                            print(pazfile, tr.stats.station)
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
                            
                            staloc = {'latitude':stla, 'longitude':stlo}
                            
                        # just use IRIS dataless volume - much easier, but less transparent!
                        else:
                            if tr.stats.network == 'AU':
                                paz = au_parser.get_paz(seedid,start_time)
                                staloc = au_parser.get_coordinates(seedid,start_time)
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
                            print(tr.stats.station+'2')
                            
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
                        
                        repi = distance(pickDat['eqla'], pickDat['eqlo'], staloc['latitude'], staloc['longitude'])[0]
                        rhyp = sqrt(repi**2 + pickDat['eqdp'])
                        azim = distance(pickDat['eqla'], pickDat['eqlo'], staloc['latitude'], staloc['longitude'])[1]
                        
                        psafilename = path.split(msf[:-6])[-1] + '.' + tr.stats.channel
                        
                        write_response_spectra(tr.stats.station, pickDat['origintime'], tr.stats.sampling_rate, \
                                               T, psa, pga, pgv, psafilename, staloc['latitude'], staloc['longitude'], \
                                               pickDat['eqla'], pickDat['eqlo'], pickDat['eqdp'], pickDat['mag'], rhyp, azim, lofreq, hifreq)
                        
        except:
            print('Failed: ' + pickDat['mseed_path'])
    
