import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, \
                  argsort, loadtxt
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from mag_tools import m02mw
from calculate_magnitudes import calc_R35, calc_HB87, calc_MLM92, calc_BJ84
from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from ltsfit import lts_linefit
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

shpfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/shapefiles/NSHA13_Background/NSHA13_BACKGROUND_NSHA18_May2016.shp'

sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

zone_code = get_field_data(sf, 'CODE', 'str')
zone_trt = get_field_data(sf, 'TRT', 'str')
used_zones = set(['EBGZ', 'CBGZ', 'NCCZ'])

####################################################################################
# correct attenuation
####################################################################################
print('change me to brune output')    
def parse_filtering_data():
    
    mwdat = []
    # read parameter file
    lines = open('brune_stats.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[9]), 'maxf': float(dat[10]), 'qual':float(dat[11]),
        	      'mw': float(dat[6]), 'sd': float(dat[7])}
    
        mwdat.append(filt)
    
    return mwdat
    

   
####################################################################################
# set def params
####################################################################################


###############################################################################
# load pickle 
###############################################################################

#recs = pickle.load(open('onger_fft_data.pkl', 'rb' ))
recs = pickle.load(open('wa_data.pkl', 'rb' ))

###############################################################################
# set data to use
###############################################################################

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX'])
# get stas to ignore
ignore_stas = open('sta_ignore.txt').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])


####################################################################################
# start main
####################################################################################

# loop through & plot station  data
cs = get_mpl2_colourlist()
events = unique(dictlist2array(recs, 'evdt'))
stations = unique(dictlist2array(recs, 'sta'))
rhyps = dictlist2array(recs, 'rhyp')

# sort by rhyps
rhyp_sort_idx = argsort(rhyps)

# get event filters
mwdat = parse_filtering_data()

events_dict = []

sp = 0
ii = 1	
magcsv = 'EVENT,ML_2800,ML_2080,MW,SD\n'
stacsv = 'EVENT,STA,RHYP,ML_2800,ML_2080,MW\n'

print('need to get mag region')

ml_array = []
mw_array = []
sd_array = []

for e, event in enumerate(events): # [::-1]): #[-2:-1]:
    print(event)
    
    # get upper & lower f for filtering
    qual = 0
    mw = nan
    sd = nan
    for fdat in mwdat:
        if fdat['ev'] == event:
            minf = fdat['minf']
            maxf = fdat['maxf']
            qual = fdat['qual']
            if qual == 1.0:
                mw = fdat['mw']
                sd = fdat['sd']
            else:
                mw = nan
                sd = nan

    #print(','.join((event,str(minf),str(maxf))))

    i = 0
    pltboxes = True
    
    bj84_2800 = []
    bj84_2080 = []
    mlm92_2800 = []
    mlm92_2080 = []
    evstas = []
    
    for rsi in rhyp_sort_idx:
        rec = recs[rsi]
        #for rec in recs:
        if rec['net'] in keep_nets:
            if not rec['sta'] in ignore_stas:
                if len(rec['channels']) > 0:
                    if rec['evdt'] == event and rec['rhyp'] <= 800:
                        print('   '+rec['sta'])
                        
                        # get NSHA23 preferred amp
                        if rec['mag'] < 4.0:
                            wa2800 = rec['wa_data']['wa_amp_2800_0.75']
                            wa2080 = rec['wa_data']['wa_amp_2080_0.75']
                        
                        elif rec['mag'] < 6.0:
                            wa2800 = rec['wa_data']['wa_amp_2800_0.5']
                            wa2080 = rec['wa_data']['wa_amp_2080_0.5']
                            
                        else:
                            try:
                                wa2800 = rec['wa_data']['wa_amp_2800_0.02']
                                wa2080 = rec['wa_data']['wa_amp_2080_0.02']
                            except:
                                wa2800 = rec['wa_data']['wa_amp_2800_0.2']
                                wa2080 = rec['wa_data']['wa_amp_2080_0.2']
                            
                        #bj84_2800.append(calc_BJ84(1, log10(wa2800), rec['rhyp']) + 0.18) # added 0.18 as mean H-V correction - see hv_ratio.png in dropbox
                        m92_2800 = calc_MLM92(0, log10(wa2800), rec['rhyp'])
                        mlm92_2800.append(m92_2800)
                        
                        #bj84_2080.append(calc_BJ84(1, log10(wa2080), rec['rhyp']) + 0.18) # added 0.18 as mean H-V correction - see hv_ratio.png in dropbox
                        m92_2080 = calc_MLM92(0, log10(wa2080), rec['rhyp'])
                        mlm92_2080.append(m92_2080)
                        
                        evstas.append(rec['sta'])
                        
                        stacsv += ','.join((str(event), rec['sta'], str('%0.1f' % rec['rhyp']), \
                                            str('%0.2f' % m92_2800), str('%0.2f' % m92_2080), \
                                            str('%0.2f' % mw))) + '\n'

    if len(mlm92_2080) >= 3:
        ml_2800 = trim_mean(array(mlm92_2800), 0.1)
        ml_2080 = trim_mean(array(mlm92_2080), 0.1)
    else:
        ml_2800 = nan
        ml_2080 = nan

    magcsv += ','.join((str(event), str('%0.2f' % ml_2800), str('%0.2f' % ml_2080), \
                        str('%0.2f' % mw), str('%0.5f' % sd))) + '\n'
    
    #print(evstas)
    #print(mlm92_2080)
    
    ml_array.append(ml_2800)
    mw_array.append(mw)
    sd_array.append(sd)

ml_array = array(ml_array)
mw_array = array(mw_array)
sd_array = array(sd_array)
##########################################################################################
# plot
fig = plt.figure(1, figsize=(10.5,9))

plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
#plt.plot(ml_array, mw_array, 'o', c='darkorange', label='Data')

# plt 2023 simulated (2800) relationship
xplt = arange(3, 6.81, 0.01)
s0, s1, s2 = loadtxt('mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)
yplt = s0 * xplt**2 + s1 * xplt + s2
plt.plot(xplt, yplt, '-', lw=2, c='k', label='2023 Simulated (W-A 2800)')

cm = plt.cm.get_cmap('RdBu_r')
sc = plt.scatter(ml_array, mw_array, c=log10(sd_array), vmin=-1.5, vmax=1.5, s=36, cmap=cm, label='Data')

plt.xlim([2.5,7])
plt.ylim([2.5,7])
plt.xlabel('ML2800', fontsize=18)
plt.ylabel('Brune MW', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc=2, fontsize=14)
plt.grid(which='both')

# make colourbar
cb = plt.colorbar(sc)

logticks = cb.get_ticks()[0:]
labels = [str('%0.1f' % 10**x) for x in logticks]
labels[0] = str('%0.2f' % 10**logticks[0])
labels[1] = str('%0.2f' % 10**logticks[1])
labels[2] = str('%0.2f' % 10**logticks[2])
labels[3] = str('%0.2f' % 10**logticks[3])
labels[4] = str('%0.2f' % 10**logticks[4])

#ticks = arange(0,len(logstressbins))
cb.set_ticks(logticks)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=270, fontsize=18, labelpad=22) # 


"""                        
##########################################################################################

# export Brune data
pklfile = open('brune_data.pkl', 'wb')
pickle.dump(events, pklfile, protocol=-1)
pklfile.close()

# write csv
txt = 'EVENT,LON,LAT,DEP,OMAG,OMAG_TYPE,BRUNE_MAG,STRESS_DROP,CORN_FREQ,FMIN,FMAX,QUALITY\n'

for ev in events_dict:
    txt += ','.join((str(ev['evdt']),str(ev['lon']),str(ev['lat']),str(ev['dep']),str(ev['omag']),ev['omag_type'], \
                     str(ev['brune_mw']),str(ev['brune_sd']),str(ev['brune_f0']),str(ev['minf']),str(ev['maxf']),str(ev['qual']))) + '\n'
"""
# write to file
f = open('ml_mw_stats.csv', 'w')
f.write(magcsv)
f.close()

# write to file
f = open('ml_sta_stats.csv', 'w')
f.write(stacsv)
f.close()


# now show figs 
plt.savefig('ml_vs_mw_brune.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

   