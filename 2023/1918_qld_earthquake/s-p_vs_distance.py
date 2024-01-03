from obspy import read, UTCDateTime
from obspy.taup import TauPyModel
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from data_fmt_tools import remove_low_sample_data, return_sta_data, remove_acceleration_data
from response import get_response_info, paz_response, deconvolve_instrument
from misc_tools import listdir_extension, savitzky_golay, get_binned_stats
from mapping_tools import distance, km2deg
from mag_tools import get_au_ml_zone
from io_catalogues import parse_ga_event_query
from os import path, chmod, stat, getcwd
from numpy import arange, sqrt, pi, exp, log, logspace, interp, nan, where, isnan, nanmean, array, polyfit
from datetime import datetime, timedelta
import pickle
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

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

###############################################################################
# set velocity mdel
###############################################################################
print('Getting velocity model ...')
plt_repi = arange(0, 2500, 10)

'''
model = TauPyModel(model="ak135")
sp_ak135 = []

for pr in plt_repi:
    arrival_dep = 10
    rngdeg = km2deg(pr)
    
    arrivals = model.get_travel_times(source_depth_in_km=arrival_dep, distance_in_degree=rngdeg)

    # find P and S
    p = []
    s = []
    for a in arrivals:
        if a.name.upper() == 'P':
            pTravelTime = a.time
        if a.name.upper() == 'S':
            sTravelTime = a.time
    
    sp_ak135.append(sTravelTime - pTravelTime)   

pklfile = open('sp_ak135.pkl', 'wb')
pickle.dump(sp_ak135, pklfile, protocol=-1)
pklfile.close()
'''
sp_ak135 = pickle.load(open('sp_ak135.pkl', 'rb' ))     
   
'''
# do IASPEI
model = TauPyModel(model="iasp91")
sp_iasp91 = []

for pr in plt_repi:
    arrival_dep = 10
    rngdeg = km2deg(pr)
    
    arrivals = model.get_travel_times(source_depth_in_km=arrival_dep, distance_in_degree=rngdeg)

    # find P and S
    p = []
    s = []
    gotP = False
    gotS = False
    for a in arrivals:
        if a.name.upper() == 'P' and gotP == False:
            pTravelTime = a.time
            gotP = True
        if a.name.upper() == 'S' and gotS == False:
            sTravelTime = a.time
            gotS = True
    
    sp_iasp91.append(sTravelTime - pTravelTime)   

pklfile = open('sp_iasp91.pkl', 'wb')
pickle.dump(sp_iasp91, pklfile, protocol=-1)
pklfile.close()
'''
sp_iasp91 = pickle.load(open('sp_iasp91.pkl', 'rb' ))        	

################################################################################
# find eevents - should have done this in pick files
################################################################################

#evdict = parse_ga_event_query('au_ge_4.4_earthquakes_export_edit.csv')

def get_ev_deets(fft_datetime):
    for evnum, ev in enumerate(evdict): 
        #ev['datetime'] = UTCDateTime(2009,3,18,5,28,17)
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
               magtype = ev['magType']
               evname = ev['description']
               mag = ev['mag'] # pref_mag
               
    return mag, magtype, evname


################################################################################
# get pick files
################################################################################
folder = '../au_stress_drop/record_picks'
#folder = 'new_picks' # for testing
pickfiles = listdir_extension(folder, 'picks')

	
################################################################################
# loop through pick files
################################################################################

print('Looping thru pick files ...')
sp_diff = []
rhyp = []
dep = []     

start_idx = 0
#pickfiles = ['2023-01-05T05.08.AU.ONGER.picks']
for p, pf in enumerate(pickfiles[start_idx:]):
    pickDat = parse_pickfile(pf)
    
    if not isnan(pickDat['mag']):
        # check ml zones
        reg = get_au_ml_zone([pickDat['eqlo']], [pickDat['eqla']])
        #print(reg)
        
        sp = pickDat['spk'] - pickDat['ppk']
        if sp < 20. and pickDat['rhyp'] > 300.:
            dummy = 0
        elif sp < 30. and pickDat['rhyp'] > 500.:
            dummy = 0
        else:
            if reg == 'WCA':
                
                sp_diff.append(sp)
                rhyp.append(pickDat['rhyp'])
                dep.append(pickDat['eqdp'])
    
sp_diff = array(sp_diff)
rhyp = array(rhyp)
dep = array(dep)
repi = sqrt(rhyp**2 - dep**2)

################################################################################
# fit quadratic
################################################################################

medamp, stdbin, medx, binstrp, nperbin = get_binned_stats(plt_repi, repi, sp_diff)

# fit quadratic

sp_coeff = polyfit(medx, medamp, 2)

from scipy import odr

data = odr.Data(medx, medamp)
odr_obj = odr.ODR(data, odr.quadratic)
out = odr_obj.run()
sp_coeff = out.beta


print(sp_coeff)

sp_fit = sp_coeff[0] * plt_repi**2 + sp_coeff[1] * plt_repi + sp_coeff[2]

################################################################################
# plot
################################################################################

fig = plt.figure(1, figsize=(8,12))
plt.subplot(2,1,1)
plt.plot(repi, sp_diff, '+', c='0.6', lw=0.5, label='WCA Data')

plt.plot(plt_repi, sp_iasp91, '-', c='orange', lw=2, label='IASPEI91')
plt.plot(plt_repi, sp_ak135, '-', c='green', lw=2, label='AK135')
plt.plot(plt_repi, sp_fit, 'k-', lw=2, label='WCA Data Fit')

plt.grid(which='both')
plt.legend(loc=2, numpoints=1)
plt.xlim([0, 2300])

plt.ylabel('S-P Time (s)')

################################################################################
# PLOT RESIDUAL
################################################################################

isapei91_res = sp_fit - sp_iasp91
ak135_res = sp_fit - sp_ak135

plt.subplot(4,1,3)
plt.plot([0, 2300], [0, 0], 'k--', lw=1)
plt.plot(plt_repi, isapei91_res, '-', c='orange', lw=2)
plt.plot(plt_repi, ak135_res, '-', c='green', lw=2)

plt.xlabel('Epicentral Distance (km)')
plt.ylabel('AU - Other Residual (s)')
plt.xlim([0, 2300])
plt.grid(which='both')

plt.savefig('au_s-p_vs_global_velocity.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
