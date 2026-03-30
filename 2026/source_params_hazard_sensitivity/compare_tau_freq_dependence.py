import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack, nanstd, loadtxt, nanmedian
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
from gmt_tools import cpt2colormap
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

ignore_stas = open('../../2023/au_stress_drop/sta_ignore.txt').readlines()
#ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('../../2023/au_stress_drop/fft_data.pkl', 'rb' ))
#recs = pickle.load(open('fft_data_mag_match.pkl', 'rb' ))
chan = recs[0]['channels'][0]
freqs = recs[0][chan]['freqs']

# load final Mw estimates
lines = open('brune_stats_cluster.csv').readlines()[1:]
brunedat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'datetime':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'qual':int(float(dat[-2])), \
           'lon':float(dat[2]), 'lat':float(dat[3]), 'sd':float(dat[10]), 'cluster':int(float(dat[-1]))}
    
    brunedat.append(tmp)
    
def get_brune_deets(fft_datetime):
    bruneStats = {'mw':nan, 'qual':0}
    for evnum, ev in enumerate(brunedat): 
        #ev['datetime'] = UTCDateTime(2024,2,27,16,4,9)
        
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
            
               bruneStats = ev
    if isnan(bruneStats['mw']):
        print('Event not found: ',fft_datetime)
                       
    return bruneStats

# load atten coeffs
print('run compare_tau_freq_dependence.py ../../2023/au_stress_drop/atten_coeffs_1.3_5km.pkl')
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))
 
###############################################################################
# do ml-mw conversions
###############################################################################
 
s0, s1, s2 = loadtxt('../../2023/au_stress_drop/mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)   

lines = open('ml_mw_bias_cluster.csv').readlines()
cmean = []
for line in lines[1:]:
    dat = line.strip().split(' +- ')[0].split(',')
    cmean.append(float(dat[1]))
    
print(cmean)

'''
cdata = loadtxt('ml_mw_bias_cluster.csv', delimiter=',', skiprows=1) 
cmean = cdata[:,1]
'''

def get_ml2mw_cluster(ml, cluster=0):
    #print(ml)
    return s0 * ml**2 + s1 * ml + s2 - cmean[cluster+1], cmean[cluster+1] # add 1 because first row is all data


# load ml-mw data from "get_cluster_ml_mw_bias.py"
csvfile = '../../2023/au_stress_drop/ml_mw_stats.csv'

print(csvfile)
lines = open(csvfile).readlines()
ml = []
mldate = []

for line in lines[1:]:
    dat = line.strip().split(',')
    mldate.append(dat[0])
    ml.append(float(dat[2]))
    
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
'''
# set brune mags
print('Setting brune mags and MLs ...')
for i, rec in enumerate(recs):
    bruneStats = get_brune_deets(UTCDateTime(rec['evdt']))
    
    recs[i]['mwb'] = bruneStats['mw']
    recs[i]['qual'] = bruneStats['qual']
    recs[i]['sd'] = bruneStats['sd']
    recs[i]['mag_cluster'] = bruneStats['cluster']
    
    # lmatch events
    ml_match = nan
    for k, mld in enumerate(mldate):
        if UTCDateTime(rec['evdt']) == mld:
            ml_match = ml[k]
            #print(ml_match)
    recs[i]['ml'] = ml_match
'''
###############################################################################
# loop thru freqs
###############################################################################
fidxs = arange(49, 135,3)

def get_inter_event_terms(magType, recs):
    '''
    magType: 1 = Brune mags
             2 = ML-MW cluster correction (MW - dMW) where dMW from mean MW - MW(conv)
             3 = Use converted magnitude
    '''
    resDict = []
    for f, c in enumerate(coeffs):
    
        if f in set(fidxs):
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
            revdt = []
            for i, rec in enumerate(recs):
                
                # if qual == 1, "mag" is Mwb
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
                    
                                #magterm = get_magnitude_term(bruneStats['mw'], c)
                                if magType == 1:
                                    prefmw = rec['mag']
                                elif magType == 2:
                                    mwconv, cor = get_ml2mw_cluster(rec['ml2800'], cluster=rec['cluster']) # only use to get correction
                                    prefmw = rec['mag'] - cor 
                                    #print(rec['mag'], cor, prefmw, rec['cluster'],rec['evdt'])
                                    
                                elif magType == 3:
                                    prefmw, cor = get_ml2mw_cluster(rec['ml2800'], -1) # subtract 1 as we want index 0
                                    
                                    
                                # get mag term
                                magterm = get_magnitude_term(prefmw, c)
                                #print(magterm)
                                
                                # get distance term
                                distterm = get_distance_term(rec['rhyp'], c)
                                
                                #	get distance independent kappa
                                kapterm = get_kappa_term(rec['sta'], c['freq'])
                                
                                #	get regional term
                                regterm = get_regional_term(rec['rhyp'], c, rec['eqdom'])
                                
                                # get total correction
                                ypred = magterm + distterm + kapterm + regterm
                                #print(ypred)
                    
                                yobs = log10(rec[channel]['swave_spec'][f])
                                yres.append(yobs - ypred)
                                rmags.append(rec['mag'])
                                revent.append(rec['ev'])
                                revdt.append(rec['evdt'])
                                rsd.append(rec['sd'])
                        except:
                            # do nothing
                            dummy = 0
                            
                        
            resData = {'yres':array(yres), 'mags':array(rmags), 'ev':array(revent), \
                       'rhyps':array(rhyps), 'rsd':array(rsd), 'evdt':array(revdt)}
            
            resDict.append(resData)
                    
    return resDict
    
###############################################################################
# calculate different scenarios
###############################################################################
'''
# ergodic
print('Getting raw residual data ...')
resDict = get_inter_event_terms(1, recs)
pklfile = open('residual_data_fdep.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

#  mw using regional correction
print('Getting cluster mag residual data ...')
resDict = get_inter_event_terms(2, recs)
pklfile = open('residual_ml_cluster_data_fdep.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

# ML2MW no corrections
print('Getting ml2mw residual data ...')
resDict = get_inter_event_terms(3, recs)
pklfile = open('residual_ml_nocluster_data_fdep.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()
'''

resDictRaw = pickle.load(open('residual_data_fdep.pkl', 'rb' ))
resDictCluster = pickle.load(open('residual_ml_cluster_data_fdep.pkl', 'rb' ))
resDictNoCluster = pickle.load(open('residual_ml_nocluster_data_fdep.pkl', 'rb' ))

###############################################################################
# get within event sigmas
###############################################################################
def get_inter_event_residuals(resDict):
    
    sigmaDict = []
    sigma_be = []
    sigma_we = []
    mags_list = []
    sds_list = []
    ev_list = []
    ev_term_list = []

    for i, resData in enumerate(resDict):
        
        # get unique events
        events = array(resData['ev'])
        uevents = unique(array(resData['ev']))
        evdt = array(resData['evdt'])
        uevdt = unique(array(resData['evdt']))
        yres = array(resData['yres'])
        mags = array(resData['mags'])
        rhyps = array(resData['rhyps'])
        sds = array(resData['rsd'])
        #print(sds)
        #print(i, len(uevents))
        
        yres_weterm = zeros_like(yres)
        
        # loop thru events
        event_terms = []
        event_mags = []
        event_sds = []
        uevdt = []
        
        for ue in uevents:
            
            #print(ue)
            # get records for given event
            ridx = where(events == ue)[0]
            event_mag = mags[ridx][0]
            
            #print(event_mag)
            event_mags.append(event_mag)
            event_sds.append(sds[ridx][0])
            uevdt.append(evdt[ridx][0])
            
            # get event terms
            event_term = nanmean(log(10**yres[ridx]))
            event_terms.append(event_term)
            
            # remove event terms
            yres_weterm[ridx] = log(10**yres[ridx]) - event_term
        
        sigma_be.append(nanstd(array(event_terms)))
        sigma_we.append(nanstd(array(yres_weterm)))   
        mags_list.append(array(event_mags))
        sds_list.append(array(event_sds)) 
        ev_list.append(array(uevdt)) 
        ev_term_list.append(array(event_terms))
                
    return sigma_be, sigma_we, mags_list, sds_list, ev_list, ev_term_list
    
###############################################################################
# get stress drop correction
###############################################################################

plt_freqs = freqs[fidxs]

sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, ev_term_list \
    = get_inter_event_residuals(resDictRaw)

sdreg_coeffs = []  
sd_corrected_event_terms = []  
sd_corrected_tau = []
clust_sd_corrected_tau = []
ergodic_tau = []
for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(ev_term_list[i]) == False)[0]
    
    ergodic_tau.append(nanstd(ev_term_list[i][idx]))
    
    sdreg = linregress(log10(sds_list[i][idx]), ev_term_list[i][idx])
    
    sdreg_coeffs.append(sdreg)
    
    # correct event terms for SD
    sd_corrected_event_terms.append(ev_term_list[i][idx] - (sdreg[0] * log10(sds_list[i][idx]) + sdreg[1]))
    sd_corrected_tau.append(nanstd(ev_term_list[i][idx] - (sdreg[0] * log10(sds_list[i][idx]) + sdreg[1])))
    
    ###############################################################################
    # get cluster SD correction
    ###############################################################################
    
    # get clusters
    clusters = []
    for ev in ev_list[i]:
        bruneStats = get_brune_deets(UTCDateTime(ev))
        clusters.append(bruneStats['cluster'])

    uclusters = unique(clusters)
    #clust_logsd = array(clust_logsd)

    mean_cluster_log_sd = []
    for cluster in uclusters:
        idx = where(clusters == cluster)[0]
        
        mean_cluster_log_sd.append(nanmean(log10(sds_list[i][idx])))

    # now get regional correction    
    ev_term = ev_term_list[i].copy()
    sd_clust_cor = zeros_like(ev_term)
    for j, cluster in enumerate(uclusters):
        idx = where(clusters == cluster)[0]
        sd_clust_cor[idx] = ev_term[idx] - (sdreg[0] * mean_cluster_log_sd[j] + sdreg[1])
        
    clust_sd_corrected_tau.append(nanstd(sd_clust_cor))

###############################################################################
# now get alt mag clusterd taus
###############################################################################


sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, ev_term_list \
    = get_inter_event_residuals(resDictNoCluster)

mconv_tau = []
for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(ev_term_list[i]) == False)[0]
    
    mconv_tau.append(nanstd(ev_term_list[i][idx])) 


sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, ev_term_list \
    = get_inter_event_residuals(resDictCluster)

delta_m_tau = []
for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(ev_term_list[i]) == False)[0]
    
    delta_m_tau.append(nanstd(ev_term_list[i][idx])) 
    
###############################################################################
# now plot
###############################################################################
ncols = 6
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)
cs = (cmap(arange(ncols)))

syms = ['o', '^', 's', 'd', 'v', '<', 'h', '>', 'p']

plt.clf()
fig = plt.figure(figsize=(6,4))
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

plt.semilogx(plt_freqs, ergodic_tau, syms[0], ls='-', c=cs[0], lw=2, \
             ms=6, mec=cs[0], mfc='w', mew=2, markevery=2, label='Ergodic')
plt.semilogx(plt_freqs, mconv_tau, syms[1], ls='-', c=cs[1], lw=2, \
             ms=6, mec=cs[1], mfc='w', mew=2, markevery=2, label='$\mathregular{M_{conv}}$')
plt.semilogx(plt_freqs, sd_corrected_tau, syms[2], ls='-', c=cs[2], lw=2, \
             ms=6, mec=cs[2], mfc='w', mew=2, markevery=2, label=r"$\Delta\sigma_i$ adjustment")
plt.semilogx(plt_freqs, delta_m_tau, syms[3], ls='-', c=cs[3], lw=2, \
             ms=6, mec=cs[3], mfc='w', mew=2, markevery=2, label=r"M($\Delta$" + '$\mathregular{M_k}$' + ') adjustment')
plt.semilogx(plt_freqs, clust_sd_corrected_tau, syms[4], ls='-', c=cs[4], lw=2, \
             ms=6, mec=cs[4], mfc='w', mew=2, markevery=2, label=r"$\Delta\sigma_k$ adjustment")

plt.xlabel('Frequency (Hz)', fontsize=15)
plt.ylabel(r"$\tau$", weight='bold', fontsize=18)
plt.legend(loc=2, fontsize=10, ncol=2, numpoints=1)
plt.grid(which='both')
plt.xlim([0.3, 10])
plt.ylim([0.2, 1.0])
plt.savefig('figures/freq_vs_tau.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()
