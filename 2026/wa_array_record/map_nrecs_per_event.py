import pickle
from sys import argv
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, \
                  argsort, percentile
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from mag_tools import m02mw
from mapping_tools import distance, drawshapepoly
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
#from ltsfit import lts_linefit
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors, colorbar
from obspy import UTCDateTime
from mpl_toolkits.basemap import Basemap
from misc_tools import remove_last_cmap_colour, listdir_extension
#from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon

mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

shpfile = 'WAarray_ply2025.shp'

sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))


    
def parse_filtering_data():
    
    filtdat = []
    # read parameter file
    #lines = open('brune_stats.csv').readlines()[1:]
    lines = open('../../2026/source_params_hazard_sensitivity/brune_stats.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'minf': float(dat[15]), \
        	      'maxf': float(dat[16]), 'qual':float(dat[17]), 'lat': float(dat[3]), 'lon':float(dat[2])}
    
        filtdat.append(filt)
    
    return filtdat
    
def correct_atten(rec, coeffs, kapdat):
    channel = rec['channels'][0]
    raw_fds = rec[channel]['swave_spec']#[:-8]
    sn_ratio = rec[channel]['sn_ratio']#[:-8]
    freqs = rec[channel]['freqs']#[:-8]
    
    sn_thresh = 4. # used 4 in model regression
    
    # loop thru freqs
    distterms = []
    regterms = []
    for c in coeffs:
        # get distance term
        distterm = get_distance_term(rec['rhyp'], c)
        
        distterms.append(distterm)
        
        regterm = get_regional_term(rec['rhyp'], c, rec['eqdom'])
        
        regterms.append(regterm)
        
    # set kappa
    kappa = kapdat[-1]['kappa0'] # default kappa
    
    # get site kappa
    for kap in kapdat:
        if kap['sta'] == rec['sta']:
            if not isnan(kap['kappa0']):
                kappa = kap['kappa0'] # + kap['kappa_r'] * rec['rhyp'] 
                #print(kap['kappa0'])
    
    k_term = log10(exp(-1 * pi * freqs * kappa))
    
    # correct to source
    cor_fds = 10**(log10(raw_fds) - distterms - regterms - k_term)
    
    # get data exceeding SN ratio
    idx = where(sn_ratio < sn_thresh)[0]
    cor_fds_nan = cor_fds.copy()
    cor_fds_nan[idx] = nan
    
    # ensure high SN for periods affected by 2ndary microseisms
    #idx = where((sn_ratio < 20) & (freqs >= 0.09) & (freqs <= 0.35))[0]
    #cor_fds_nan[idx] = nan
    
    # if short period only use f > 0.5
    if  channel.startswith('SH') or channel.startswith('EH'):
        if rec['pazfile'].endswith('s6000-2hz.paz'):
            idx = where(freqs < 0.9)[0]
        else:
            idx = where(freqs < 0.5)[0]
        cor_fds_nan[idx] = nan
        
    # ignore dodgy CMSA data
    if rec['sta'] == 'CMSA_FIXED' or rec['sta'] == 'PI207': 
        idx = where(freqs < 0.5)[0]
        cor_fds_nan[idx] = nan
        
    if rec['sta'] == 'NPS': 
        idx = where(freqs < 100)[0]
        cor_fds_nan[idx] = nan
        
    if rec['sta'] == 'STKA': 
        idx = where(freqs < 0.1)[0]
        cor_fds_nan[idx] = nan
        
    if rec['sta'] == 'AS17' and rec['evdt'].year > 2013 and rec['evdt'].year < 2018: 
        idx = where(freqs < 100)[0]
        cor_fds_nan[idx] = nan
        
    # if sta starts with AQT
    if rec['net'].startswith('1Q'):
        print(rec['sta'])
        idx = where(freqs > 5.0)[0]
        cor_fds_nan[idx] = nan
            
    return cor_fds, cor_fds_nan, freqs

# now fit Brune model
def fit_brune_model(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(c[0] / (1 + (f / (c[1]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve


    
def fit_brune_model_fixed_omega(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 4.2 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_petermann(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 0.69 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_murrayville(c, f):
    from numpy import array, log
    
    fixed_omega = 0.0095 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve

def fit_brune_model_fixed_omega_carnarvon(c, f):
    from numpy import array, log
    
    fixed_omega = 0.6 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_leongatha(c, f):
    from numpy import array, log

    fixed_omega = 0.01 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_marblebar(c, f):
    from numpy import array, log
    
    fixed_omega = 0.038 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_2020bowen(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 0.0081 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    
    return FittedCurve
   

###############################################################################
# load pickle 
###############################################################################

#recs = pickle.load(open('onger_fft_data.pkl', 'rb' ))
recs = pickle.load(open('../../2023/au_stress_drop/fft_data.pkl', 'rb' ))

# convert mags to MW
for i, rec in enumerate(recs):
    if rec['magType'].startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(datetimes[i]) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.05
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])

# load atten coeffs
pklfile = '../../2023/au_stress_drop/fft_data.pkl'

ignorePoorQual = True # True or False
if ignorePoorQual == 'True':
    ignorePoorQual = True
else:
    ignorePoorQual = False

coeffs = pickle.load(open(pklfile, 'rb' ))

###############################################################################
# set data to use
###############################################################################
print('run map_nrecs_per_event.py ../../2023/au_stress_drop/fft_data.pkl True')
# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'M8', 'DU', 'WG', '4N', \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX', 'AM', 'YW', '3B', '1K', \
                 '1Q', '3O', '7F', '6K', '5G', '5C'])
                 

# get stas to ignore
ignore_stas = open('../../2023/au_stress_drop/sta_ignore.txt').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

'''
SWN20, SWNNG

'''


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


# set M-R lookup
mdist_lookup_mags = arange(3.25,7.1,0.5)
mdist_lookup_dists = array([550, 1200, 1700, 2000, 2200, 2200, 2200, 2200])

# get event filters
filtdat = parse_filtering_data()

events_dict = []


ev_cnt = []
ev_mw = []
ev_qual = []
lats = []
lons = []

for e, event in enumerate(events): #[::-1]): #[-2:-1]:
    
    # get upper & lower f for filtering
    qual = 1
    for fdat in filtdat:
        if fdat['ev'] == event:
            minf = fdat['minf']
            maxf = fdat['maxf']
            qual = fdat['qual']
            mag = fdat['mw']
            lat = fdat['lat']
            lon = fdat['lon']
            

    #print(','.join((event,str(minf),str(maxf))))
    
    # check if we want to plot 
    plotTrue = False
    if ignorePoorQual == True:
        if qual == 1:
            plotTrue = True
            
    else:
        plotTrue = True
        
    if plotTrue == True:    
        print(event)
        
        i = 0
        pltboxes = True
        for rsi in rhyp_sort_idx:
            rec = recs[rsi]
            #for rec in recs:
            if rec['net'] in keep_nets:
                #if not rec['sta'] in ignore_stas:
                    if len(rec['channels']) > 0:
                        if rec['evdt'] == event: # will need to cahnge to rec['datetime']
                            print('   '+rec['sta'])
                            
                            # get distance cut-off
                            idx = where((rec['mag']) >= mdist_lookup_mags)[0]
                            if len(idx) == 0:
                                mag_dist = mdist_lookup_dists[0]
                            else:
                                mag_dist = mdist_lookup_dists[idx[-1]]
                            
                            # do skip sta checks
                            skip_sta = False
                            if rec['sta'] == 'QIS' and rec['evdt'].year <= 1999:
                                skip_sta = True
                            elif rec['sta'] == 'RMQ' and rec['evdt'].year <= 1998:
                                skip_sta = True
                                
                            chan = rec['channels'][0]
                            snr = rec[chan]['sn_ratio'][99] # 2 Hz
                            
                            if snr < 4.0 and rec['rhyp'] > mag_dist:
                                skip_sta = True
                                
                            #print(skip_sta, rec['rhyp'], rec['mag'], mag_dist)
                            if skip_sta == False:
                                i += 1
                                
    ev_cnt.append(i)        
    ev_mw.append(mag)
    ev_qual.append(qual)
    lats.append(lat)
    lons.append(lon)
     
##########################################################################################
# get sta data
##########################################################################################
sta_cnt = []
sta_net = []
sta_lat = []
sta_lon = []

for sta in stations:
    net = ''
    i = 0
    stla = nan
    stlo = nan
       
    for rec in recs:
       skip_sta = False
       
       
       if len(rec['channels']) > 0 and rec['sta'] == sta:
           # get distance cut-off
           idx = where((rec['mag']) >= mdist_lookup_mags)[0]
           if len(idx) == 0:
               mag_dist = mdist_lookup_dists[0]
           else:
               mag_dist = mdist_lookup_dists[idx[-1]]
           
           # do skip sta checks
           skip_sta = False
           if rec['sta'] == 'QIS' and rec['evdt'].year <= 1999:
               skip_sta = True
           elif rec['sta'] == 'RMQ' and rec['evdt'].year <= 1998:
               skip_sta = True
               
           chan = rec['channels'][0]
           snr = rec[chan]['sn_ratio'][99] # 2 Hz
           
           if snr < 4.0 and rec['rhyp'] > mag_dist:
               skip_sta = True
               
           if skip_sta == False:
               i += 1
               net = rec['net']
               stla = rec['stla']
               stlo = rec['stlo']
    
    sta_cnt.append(i)
    sta_net.append(net)
    sta_lat.append(stla)
    sta_lon.append(stlo)
    print(sta, i)
    
##########################################################################################
# add epicentres
##########################################################################################

bounds = array([1, 5, 10, 20, 30, 50, 70, 100, 1000])
bounds = array([10, 15, 20, 30, 50, 1000])
cmap = plt.get_cmap('plasma_r', len(bounds))
cs = (cmap(arange(len(bounds))))

urcrnrlat = -25.5
llcrnrlat = -35.75
urcrnrlon = 125.
llcrnrlon = 113.5
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(12,8))
plt.tick_params(labelsize=14)
ax = plt.subplot(121)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=200.)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,3.), labels=[1,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,4.), labels=[0,0,0,1], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.fillcontinents(color='w', lake_color='0.9') #, zorder=0)
m.drawmapboundary(fill_color='0.9')            

# get zorder for plotting
sortidx = argsort(argsort(ev_mw))
for i in range(0, len(ev_mw)):
    if ev_mw[i] >= 3. and ev_qual[i] == 1.0:
        #get colour idx
        
        cidx = where(bounds >= ev_cnt[i])[0][0] # take first idx
        col= tuple(cs[cidx][:-1])

        x, y = m(lons[i], lats[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=col, markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=-10 + ev_mw[i]*5, zorder=zo, alpha=0.8)
    
# make legend
legmag = [3.5, 4.5, 5.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-10 + lm*5), alpha=1., zorder=len(ev_mw)+1, lw=2)
    legh.append(h[0])

l = plt.legend(legh, ('$\mathregular{M_{W}}$ 3.5', '$\mathregular{M_{W}}$ 4.5', '$\mathregular{M_{W}}$ 5.5'), \
               loc=2, numpoints=1, fontsize=10)
l.set_zorder(len(ev_mw)+5)


##########################################################################################
# add statioms
##########################################################################################

bounds = array([1, 3, 5, 10, 15, 20, 30, 50, 1000])
bounds = array([1, 3, 5, 10, 20, 1000])
cmap = plt.get_cmap('plasma_r', len(bounds))
cs = (cmap(arange(len(bounds))))

ax = plt.subplot(122)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=200.)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,3.), labels=[1,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,4.), labels=[0,0,0,1], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
m.fillcontinents(color='w', lake_color='0.9') #, zorder=0)
m.drawmapboundary(fill_color='0.9')

# add phase polys
drawshapepoly(m, plt, sf, edgecolor='k', alpha=1, cmap=-99, zorder=1, lw=1)

syms = ['^','v','s','d','H','<']
nets = ['WG','2P','AU','S1','Other']

for i, sta in enumerate(stations):
    if sta_cnt[i] >= 1.: # and ev_qual[i] == 1.0:
        #get colour idx
        
        cidx = where(bounds >= sta_cnt[i])[0][0] # take first idx
        col= tuple(cs[cidx][:-1])

        x, y = m(sta_lon[i], sta_lat[i])
        zo = i + 20
        
        sym = syms[4]
        for j, net in enumerate(nets):
            if net == sta_net[i]:
                sym = syms[j]
            	
        plt.plot(x, y, marker=sym, mfc=col, markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=6, zorder=zo, alpha=0.8)
    
# make legend
legmag = [4., 5., 6.]
legh = []
for i, net in enumerate(nets):
    x, y = m(0, 0)
    h = m.plot(x, y, syms[i], mfc='0.8', mec='k', markersize=6, alpha=1., zorder=len(stations)+1, lw=0.25)
    legh.append(h[0])

l = plt.legend(legh, nets, loc=1, numpoints=1, fontsize=10)
l.set_zorder(len(stations)+5)

##########################################################################################
# add colourbar 1
##########################################################################################
'''
bounds = array([1, 3, 5, 10, 20, 30, 50, 70, 100, 150, 200, 1000])
ticks = arange(0.5,len(bounds))
cnt_rng = ['1', '2-3', '4-5', '6-10', '11-20', '21-30', '31-50', '51-70', '71-100', '101-150', '151-200', '200+']
'''

#plt.gcf().subplots_adjust(bottom=0.0)
cax = fig.add_axes([0.1,0.16,0.403,0.025]) # setup colorbar axes.
norm = mpl.colors.Normalize(vmin=0, vmax=len(bounds))#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', alpha=0.8)

# set cb labels
bounds = array([10, 15, 20, 30, 50, 1000])
ticks = arange(0.5,len(bounds))
cnt_rng = ['1-10', '11-15', '16-20', '21-30', '31-50', '51+']
cb.set_ticks(ticks)
cb.set_ticklabels(cnt_rng)
cb.ax.tick_params(labelsize=12)

titlestr = 'Number of Records per Event'
cb.set_label(titlestr, fontsize=14)

##########################################################################################
# add colourbar 2
##########################################################################################

#plt.gcf().subplots_adjust(bottom=0.0)
cax = fig.add_axes([0.525,0.16,0.403,0.025]) # setup colorbar axes.
norm = mpl.colors.Normalize(vmin=0, vmax=len(bounds))#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', alpha=0.8)

# set cb labels
bounds = array([1, 3, 5, 10, 20, 1000])

ticks = arange(0.5,len(bounds))
cnt_rng = ['1-2', '3-5', '6-10', '11-15', '16-20', '21+']
cb.set_ticks(ticks)
cb.set_ticklabels(cnt_rng)
cb.ax.tick_params(labelsize=12)

titlestr = 'Number of Records per Station'
cb.set_label(titlestr, fontsize=14)

plt.savefig('recs_per_event.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('recs_per_event.pdf', fmt='pdf', dpi=300, bbox_inches='tight')
plt.show()