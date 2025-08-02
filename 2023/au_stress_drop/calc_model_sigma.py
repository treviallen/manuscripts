import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack, nanstd
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, get_log_xy_locs
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term, get_regional_term
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from datetime import timedelta
from response import stationlist2dict
import warnings
warnings.filterwarnings("ignore")
import shapefile
from sys import argv
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colorbar
import matplotlib as mpl
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42

cor_dist = 100
def fit_regional_correction(medx, logmedamp):
    idx = where(10**medx >= cor_dist)[0]

    data = odrpack.RealData(medx[idx], logmedamp[idx])

    # fit all as free params
    afit = odrpack.Model(correct_far_field)
    odr = odrpack.ODR(data, afit, beta0=[0.0])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    ffc = out.beta

    return ffc[0]
    
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'M8', 'DU', 'WG', '4N', \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX', 'AM', 'YW', '3B', '1K', \
                 '1Q', '3O', '7F', '6K', '5G', '5C'])

ignore_stas = open('sta_ignore.txt').readlines()
#ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))
chan = recs[0]['channels'][0]
freqs = recs[0][chan]['freqs']

# load final Mw estimates
lines = open('brune_stats.csv').readlines()[1:]
brunedat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'datetime':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'qual':int(float(dat[-1])), 'sd':float(dat[10])}
    
    brunedat.append(tmp)
    
def get_brune_deets(fft_datetime):
    bruneStats = {'mw':nan, 'qual':0}
    for evnum, ev in enumerate(brunedat): 
        #ev['datetime'] = UTCDateTime(2024,2,27,16,4,9)
        
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
            
               bruneStats = ev
               
    return bruneStats

# load atten coeffs
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))
    
###############################################################################
# set datasets
###############################################################################
'''
events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
stations = unique(dictlist2array(recs, 'sta'))
eqlo = dictlist2array(recs, 'eqlo')
eqla = dictlist2array(recs, 'eqla')
stlo = dictlist2array(recs, 'stlo')
stla = dictlist2array(recs, 'stla')
'''

### !!! COMMENT OUT AFTER FIRST RUN !!!!
"""
# set brune mags
print('Setting brune mags ...')
for i, rec in enumerate(recs):
    bruneStats = get_brune_deets(UTCDateTime(rec['evdt']))
    
    recs[i]['mwb'] = bruneStats['mw']
    recs[i]['qual'] = bruneStats['qual']
    recs[i]['sd'] = bruneStats['sd']


###############################################################################
# loop thru freqs
###############################################################################
resDict = []
for f, c in enumerate(coeffs):
    
    print("Coeffs Freq = " +str('%0.3f' % c['freq']))
        
    chan = recs[0]['channels'][0]
    freq = recs[0][chan]['freqs'][f]
    print("Reg Freq = " +str('%0.3f' % freq))
    
    if not freq == c['freq']:
       print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
       crash

    ###############################################################################
    # parse coefs and get model prediction
    ###############################################################################
    rhyps = []
    yres = []
    rmags = []
    revent = []
    rsd = []
    for i, rec in enumerate(recs):
        
        if rec['net'] in keep_nets and rec['qual'] == 1:
            if not rec['sta'] in ignore_stas:
                try:
                    channel = rec['channels'][0]
            
                    # filter by instrument type
                    addData = True
                    if channel.startswith('SH') or channel.startswith('EH'):
                        if rec[channel]['freqs'][f] < 0.9 and rec['pazfile'].endswith('s6000-2hz.paz'):
                            addData = False
                        elif rec[channel]['freqs'][f] < 0.4:
                            addData = False
                    
                    # filer by sample-rate
                    if rec[channel]['freqs'][f] > (0.4 * rec[channel]['sample_rate']):
                        addData = False
                        
                    if rec[channel]['sn_ratio'][f] >= 4. and addData == True:
                        rhyps.append(rec['rhyp'])
            
                        # get mag term
                        #magterm = get_magnitude_term(bruneStats['mw'], c)
                        magterm = get_magnitude_term(rec['mwb'], c)
            
                        # get distance term
                        distterm = get_distance_term(rec['rhyp'], c)
                        
                        #	get distance independent kappa
                        kapterm = get_kappa_term(rec['sta'], c['freq'])
                        
                        #	get regional term
                        regterm = get_regional_term(rec['rhyp'], c, rec['eqdom'])
                        
                        # get total correction
                        ypred = magterm + distterm + kapterm + regterm
            
                        yobs = log10(rec[channel]['swave_spec'][f])
                        yres.append(yobs - ypred)
                        rmags.append(rec['mwb'])
                        revent.append(rec['ev'])
                        rsd.append(rec['sd'])
                except:
                    # do nothing
                    dummy = 0
                
    resDict.append({'yres':array(yres), 'mags':array(rmags), 'ev':array(revent), \
    	              'rhyps':array(rhyps), 'rsd':array(rsd)})
    	
pklfile = open('residual_data.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()
"""
resDict = pickle.load(open('residual_data.pkl', 'rb' ))

###############################################################################
# get within event sigmas
###############################################################################

sigmaDict = []
sigma_be = []
sigma_we = []
for i, resData in enumerate(resDict):
    
    # get unique events
    events = array(resData['ev'])
    uevents = unique(array(resData['ev']))
    yres = array(resData['yres'])
    mags = array(resData['mags'])
    rhyps = array(resData['rhyps'])
    sds = array(resData['rsd'])
    print(i, len(uevents))
    
    yres_evterm = zeros_like(yres)
    
    # loop thru events
    event_terms = []
    event_mags = []
    event_sds = []
    for ue in uevents:
        
        # get records for given event
        ridx = where(events == ue)[0]
        event_mag = mags[ridx][0]
        event_mags.append(event_mag)
        event_sds.append(sds[ridx][0])
        
        # get event terms
        event_term = nanmean(log(10**yres[ridx]))
        event_terms.append(event_term)
        
        # remove event terms
        yres_evterm[ridx] = log(10**yres[ridx]) - event_term
    
    sigma_be.append(nanstd(array(event_terms)))
    sigma_we.append(nanstd(array(yres_evterm)))    
    
    # keep some data for plotting
    if i == 69:
        plt_event_terms1 = event_terms
        plt_event_mags1 = event_mags
        plt_event_sds1 = event_sds
        plt_yres_evterm1 = yres_evterm
        plt_mags1 = mags
        plt_rhyps1 = rhyps
    elif i == 100:
        plt_event_terms2 = event_terms
        plt_event_mags2 = event_mags
        plt_event_sds2 = event_sds
        plt_yres_evterm2 = yres_evterm
        plt_mags2 = mags
        plt_rhyps2 = rhyps

###############################################################################
# plot combined tau/phi
###############################################################################
fig = plt.figure(figsize=(14,8))
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
#plt.tight_layout() 
props = dict(boxstyle='round', facecolor='w', alpha=1)

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']

ax = plt.subplot(2,2,1)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=16)
plt.ylabel('Between-Event\n(ln Residual)', fontsize=16)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-3, 3])
pertxt = 'f = 0.75 Hz'
xpos = get_log_xy_locs([3, 7], 0.04)
ypos = (6*0.94) - 3
plt.text(xpos, ypos, pertxt, ha='left', va='top', fontsize=16, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[0], fontsize=20, va='bottom', ha='right')

###########

ax = plt.subplot(2,2,2)
plt.scatter(plt_event_mags2, plt_event_terms2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=16)
plt.ylabel('Between-Event\n(ln Residual)', fontsize=16)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-3, 3])
pertxt = 'f = 2.0 Hz'
xpos = get_log_xy_locs([3, 7], 0.04)
ypos = (6*0.94) - 3
plt.text(xpos, ypos, pertxt, ha='left', va='top', fontsize=16, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[1], fontsize=20, va='bottom', ha='right')

#############

ax = plt.subplot(2,2,3)
norm2 = mpl.colors.Normalize(vmin=3, vmax=7)
plt.scatter(plt_rhyps1, plt_yres_evterm1, c=plt_mags1, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
plt.plot([0, 2000],[0,0], 'k--', lw=1.)

plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Within-Event\n(ln Residual)', fontsize=16)

plt.xlim([0, 2000])
plt.ylim([-3, 3])
pertxt = 'f = 0.75 Hz'
xpos = get_log_xy_locs([0, 2000], 0.04)
ypos = (6*0.92) - 3
plt.text(xpos, ypos, pertxt, ha='left', va='top', fontsize=16, bbox=props)

xpos = -2000*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[2], fontsize=20, va='bottom', ha='right')

xticks = [0, 500, 700, 1000, 1500, 2000]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

############

ax = plt.subplot(2,2,4)
plt.scatter(plt_rhyps2, plt_yres_evterm2, c=plt_mags2, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
plt.plot([0, 2000],[0,0], 'k--', lw=1.)

plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Within-Event\n(ln Residual)', fontsize=16)

plt.xlim([0, 2000])
plt.ylim([-3, 3])
pertxt = 'f = 2.0 Hz'
xpos = get_log_xy_locs([0, 2000], 0.04)
ypos = (6*0.92) - 3
plt.text(xpos, ypos, pertxt, ha='left', va='top', fontsize=16, bbox=props)

xpos = -2000*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[2], fontsize=20, va='bottom', ha='right')

xticks = [0, 500, 700, 1000, 1500, 2000]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

###########

# set cbars
plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)
cax = fig.add_axes([0.91,0.58,0.02,0.3]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.plasma_r, orientation='vertical', alpha=1, norm=norm1)
ticks = [-1. , -0.5,  0. ,  0.5,  1. ,  1.5]
cb.set_ticks(ticks)
cb.set_ticklabels([str('%0.1f' % 10**x) for x in ticks])
cb.set_label('Stress Drop (MPa)', rotation=270, labelpad=20, fontsize=15)

# set cbars
cax = fig.add_axes([0.91,0.12,0.02,0.3]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.viridis_r, orientation='vertical', alpha=1, norm=norm2)
ticks = [3.0, 4.0, 5.0, 6.0, 7.0]
cb.set_ticks(ticks)
cb.set_ticklabels([str(x) for x in ticks])
cb.set_label('Moment Magnitude', rotation=270, labelpad=20, fontsize=15)

plt.savefig('sigma_model.png', fmt='png', dpi=600, bbox_inches='tight')       
plt.show()
    

# now smooth coeffs
sg_window = 21
sg_poly = 3
smooth_be = savitzky_golay(sigma_be, sg_window, sg_poly)
smooth_we = savitzky_golay(sigma_we, sg_window, sg_poly)

for f, c in enumerate(coeffs):
    coeffs[f]['sigbe'] = smooth_be[f]
    coeffs[f]['sigwe'] = smooth_we[f]
    coeffs[f]['sigto'] = sqrt(smooth_be[f]**2 + smooth_we[f]**2)
     
pklfile = open(coeff_pkl, 'wb')
pickle.dump(coeffs, pklfile, protocol=-1)
pklfile.close()

###############################################################################
# write sigmas
###############################################################################
sigtxt = 'frequency,sigma_be,sigma_we,sigma_t\n'
# write PGM first
for i, f in enumerate(freqs):
    sigtxt += ','.join((str(f), str('%0.4e' % smooth_be[i]), str('%0.4e' % smooth_we[i]), \
                        str('%0.4e' % coeffs[i]['sigto']))) + '\n'
    
f = open('sigma_model.csv', 'w')
f.write(sigtxt)
f.close()
