# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 08:36:05 2017

@author: u56903
"""

def parse_usgs_events(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
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


# reads sa data files and returns period (T) and acceleration (SA) vectors
def read_sa(safile):
    from numpy import array
    from scipy.constants import g

    lines = open(safile).readlines()
    
    safile = lines[0].strip().split('\t')[-1]
    chan = lines[0].strip().split('.')[-2]
    datestr = lines[1].strip().split('\t')[-1].replace(' ','T')    
    sta = lines[2].strip().split('\t')[-1]
    stla = float(lines[3].strip().split('\t')[-1])
    stlo = float(lines[4].strip().split('\t')[-1])
    sr = float(lines[5].strip().split('\t')[-1].split()[0])
    eqla = float(lines[6].strip().split('\t')[-1])
    eqlo = float(lines[7].strip().split('\t')[-1])
    dep = float(lines[8].strip().split('\t')[-1].split()[0])
    mag = float(lines[9].strip().split('\t')[-1])
    rhyp = float(lines[10].strip().split('\t')[-1].split()[0])
    azim = float(lines[11].strip().split('\t')[-1].split()[0])
    pga = float(lines[12].strip().split('\t')[-1].split()[0]) / (1000. * g) # convert mm/s**2 to g
    pgv = float(lines[13].strip().split('\t')[-1].split()[0]) / 10. # convert mm/s to cm/s
    net = safile.split('.')[2]

    SA = []
    T = []
    
    for line in lines[24:]:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))
    
    T = array(T)
    SA = array(SA) / g

    rec = {'chan':chan, 'datestr':datestr, 'sta':sta, 'stla':stla, 'stlo':stlo, \
           'sr':sr, 'eqla':eqla, 'eqlo':eqlo, 'dep':dep, 'mag':mag, 'rhyp':rhyp, \
           'azim':azim, 'pga':pga, 'pgv':pgv, 'per':T, 'sa':SA, 'safile':safile, \
           'network':net}
    
    return rec

################################################################################
# import funcs
################################################################################
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors, colorbar
from numpy import where, zeros, exp, log, interp, array, logspace, log10, sqrt, arange, \
                  concatenate, around, mean, polyfit, delete, isnan, nan, isfinite, unique
from misc_tools import listdir_extension, dictlist2array, get_binned_stats, savitzky_golay
#from calc_oq_gmpes import inslab_gsims
from os import path
from sys import argv
from scipy.stats import linregress
import pickle
mpl.style.use('classic')
import warnings
warnings.filterwarnings("ignore")
maxDist = 1750.

# get zone
region = argv[1]

################################################################################
# loop through sa files and get data
################################################################################
#usgscsv = '20190625_merged_events.csv'
#evdict = parse_usgs_events(usgscsv)
"""
folder = 'psa'
extension = 'psa'
safiles = listdir_extension(folder, extension)
recs = []
print('Building dataset...')
for saf in safiles:
    rec = read_sa(path.join(folder, saf))
    recs.append(rec)
    '''
    if rec['network'] == 'AU' or rec['network'] == 'S':
        recs.append(rec)
    
    elif rec['rhyp'] < 400:
        recs.append(rec)
    '''

################################################################################
# loop stations
################################################################################

stdict = []
for saf in safiles:
    newsta = True    
    fsplit = saf.split('.')
    # if station info exists
    for sd in stdict:
        if sd['filedate'] == saf[0:16] and sd['sta'] == fsplit[3]:
            if fsplit[4].endswith('E') or fsplit[4].endswith('1'):
                sd['chans'][0] = 1
                sd['chstr'][0] = fsplit[-2]
            elif fsplit[4].endswith('N') or fsplit[4].endswith('2'):
                sd['chans'][1] = 1
                sd['chstr'][1] = fsplit[-2]
            elif fsplit[4].endswith('Z'):
                sd['chans'][2] = 1
                sd['chstr'][2] = fsplit[-2]
            sd['net'] = fsplit[2]
            newsta = False
    
    # if new event & new stn
    if newsta == True:
        chans = zeros(3)
        chstr = ['', '', '']
        if fsplit[4].endswith('E') or fsplit[4].endswith('1'):
            chans[0] = 1
            chstr[0] = fsplit[-2]
        elif fsplit[4].endswith('N') or fsplit[4].endswith('2'):
            chans[1] = 1
            chstr[1] = fsplit[-2]
        elif fsplit[4].endswith('Z'):
            chans[2] = 1
            chstr[2] = fsplit[-2]
        
        rec = read_sa(path.join(folder, saf))
        datestr = rec['datestr']
        sd = {'ev':datestr, 'sta':fsplit[3], 'chans':chans, 'chstr':chstr, 'net':fsplit[2], \
        	    'filedate':saf[0:16], 'azim':0.0, 'rhyp':9999}
        
        #if not sd['sta'] == 'AS31':
        stdict.append(sd)
            
################################################################################
# now get geometric mean for each stn
################################################################################

for i, st in enumerate(stdict):
    geom = nan
    # get geometric mean of horizontals
    if st['chans'][0] == int(1) and st['chans'][1] == int(1):
        for rec in recs:
            # get east channel
            if rec['datestr'] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][0]:
                
                edat = rec['sa']
                
                # add additional info
                stdict[i]['rhyp'] = rec['rhyp']
                stdict[i]['azim'] = rec['azim']
                stdict[i]['dep'] = rec['dep']
                stdict[i]['mag'] = rec['mag']
                stdict[i]['stla'] = rec['stla']
                stdict[i]['stlo'] = rec['stlo']
                stdict[i]['eqla'] = rec['eqla']
                stdict[i]['eqlo'] = rec['eqlo']
                stdict[i]['per'] = rec['per']
                stdict[i]['pga'] = rec['pga']
                stdict[i]['pgv'] = rec['pgv']
                stdict[i]['network'] = rec['network']
                stdict[i]['date'] = rec['datestr']
                
            # get north channel
            if rec['datestr'] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][1]:
                
                ndat = rec['sa']
        
        geom = exp((log(edat) + log(ndat)) / 2.)        
        stdict[i]['geom'] = geom
    
    # else, get on or other horizontal
    elif st['chans'][0] == 1 or st['chans'][1] == 1:
        for rec in recs:
            # get east channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][0]:
                
                stdict[i]['geom'] = rec['sa']
                # add additional info
                stdict[i]['rhyp'] = rec['rhyp']
                stdict[i]['azim'] = rec['azim']
                stdict[i]['dep'] = rec['dep']
                stdict[i]['mag'] = rec['mag']
                stdict[i]['stla'] = rec['stla']
                stdict[i]['stlo'] = rec['stlo']
                stdict[i]['eqla'] = rec['eqla']
                stdict[i]['eqlo'] = rec['eqlo']
                stdict[i]['per'] = rec['per']
                stdict[i]['pga'] = rec['pga']
                stdict[i]['pgv'] = rec['pgv']
                stdict[i]['network'] = rec['network']
                stdict[i]['date'] = rec['datestr']
                
            # get north channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][1]:
                
                stdict[i]['geom'] = rec['sa']
                # add additional info
                stdict[i]['rhyp'] = rec['rhyp']
                stdict[i]['azim'] = rec['azim']
                stdict[i]['dep'] = rec['dep']
                stdict[i]['mag'] = rec['mag']
                stdict[i]['stla'] = rec['stla']
                stdict[i]['stlo'] = rec['stlo']
                stdict[i]['eqla'] = rec['eqla']
                stdict[i]['eqlo'] = rec['eqlo']
                stdict[i]['per'] = rec['per']
                stdict[i]['pga'] = rec['pga']
                stdict[i]['pgv'] = rec['pgv']
                stdict[i]['network'] = rec['network']
                stdict[i]['date'] = rec['datestr']
    
    else:
        stdict[i]['azim'] = 0.0 # skip odd recs
        #print(st['sta'])
    
    '''            
    # else, just use vertical
    elif st['chans'][2] == 1:
        for rec in recs:
            # get east channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][2]:
                
                stdict[i]['geom'] = rec['sa']
                # add additional info
                stdict[i]['rhyp'] = rec['rhyp']
                stdict[i]['azim'] = rec['azim']
                stdict[i]['dep'] = rec['dep']
                stdict[i]['mag'] = rec['mag']
                stdict[i]['stla'] = rec['stla']
                stdict[i]['stlo'] = rec['stlo']
                stdict[i]['eqla'] = rec['eqla']
                stdict[i]['eqlo'] = rec['eqlo']
                stdict[i]['per'] = rec['per']
                stdict[i]['pga'] = rec['pga']
                stdict[i]['pgv'] = rec['pgv']
                stdict[i]['network'] = rec['network']
                stdict[i]['date'] = rec['datestr']
    '''
                   
didx = []
stdict = array(stdict)
for j, st in enumerate(stdict):
    if st['chstr'][0] == '' and st['chstr'][0] == '':
        didx.append(j)
    elif st['net'] == 'GE' and st['rhyp'] > 400:
        didx.append(j)
    elif st['net'] == 'GE' and st['azim'] > 210:
        didx.append(j)
    elif st['net'] == 'GE' and st['azim'] < 150:
        didx.append(j)
    elif st['sta'] == 'PMG' and st['rhyp'] > 400:
        didx.append(j)
    elif st['sta'] == 'PMG' and st['azim'] > 210:
        didx.append(j)
    elif st['sta'] == 'PMG' and st['azim'] < 150:
        didx.append(j)
    elif st['azim'] > 240 and st['net'] != 'AEES':
        didx.append(j)
    elif st['azim'] > 225 and region != 'NGH':
        didx.append(j)
    elif st['azim'] < 135 and st['net'] != 'AEES':
        didx.append(j)
    elif st['rhyp'] > maxDist:
        didx.append(j)

del recs
print(len(stdict))
stdict = delete(stdict, array(didx))
print(len(stdict))
print(len(didx))

################################################################################
# save/load pickle
################################################################################

print('Saving pkl file...')
pklfile = open("stdict.pkl", "wb" )
pickle.dump(stdict, pklfile) #, protocol=-1)
pklfile.close()
"""
print('Loading pkl file...')
stdict = pickle.load(open("stdict.pkl", "rb" ))

################################################################################
# cross-check with S/N data
################################################################################
print('Loading pkl file...')
cnt = 0
sndict = pickle.load(open("nac_fft_data.pkl", "rb" ))
cnt = 0
for i, sd in enumerate(stdict):
    hif_limit = 0.1
    lof_limit = 1000.
    
    recFound = False
    for sn in sndict:
        '''
        if sd['sta'].startswith('DRS') and sn['sta'].startswith('DRS') and sd['ev'].startswith('1995-12-25'):
            print(sd['sta'], sd['ev'], sn['ev'])
        
        '''
        if sd['ev'] == sn['ev'].replace('.',':') and sd['sta'] == sn['sta']:
            #print(sd['ev'], sd['sta'], 'blah')
            #print(sd['geom'])
            
            for chan in sn['channels']:
                if chan.endswith('E') or chan.endswith('N'):
                    if sn[chan]['hif_limit'] > hif_limit:
                       hif_limit = sn[chan]['hif_limit']
                       
                    if sn[chan]['lof_limit'] < lof_limit:
                       lof_limit = sn[chan]['lof_limit']
                       
                    if sd['chstr'][0].startswith('SH') and lof_limit < 0.25:
                        lof_limit = 0.25
                    
                    '''
                    # don't use wrab, but keep for site response
                    if sd['sta'] == 'WRAB':
                        lof_limit = 1000.
                        hif_limit = 0.01
                    
                    # don't use wrab, but keep for site response
                    if sd['sta'] == 'COEN':
                        lof_limit = 0.5
                    '''    
            recFound = True
            
    # remove noisy psa data
    idx = where(sd['per'] > (1./lof_limit))[0]
    stdict[i]['geom'][idx] = nan
    
    idx = where(sd['per'] < (1./hif_limit))[0]
    stdict[i]['geom'][idx] = nan
    
    
    if recFound == False:
        print('Not found: '+sd['ev'] + ' ' + sd['sta'])
        cnt += 1
   
del sndict 
################################################################################
# setup inversions
################################################################################
import scipy.odr.odrpack as odrpack

if region == 'BS':
    xref = 650.
    xref_hinge = 650.
    minDist = 100
elif region == 'NGH':
    xref = 1000 # NGH
    xref_hinge = 800.
    minDist = 10
elif region == 'OBE':
    xref = 1000 # OBE
    xref_hinge = 1000 # OBE
    minDist = 100

mrng = arange(5.3, 7.9, 0.1)
mpltrng = 0.05

def fit_atten(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] - c[1]*(x) - c[2]*log10(x)
    
    return ans
    
def fit_gr(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] + c[1]*log10(x)
    
    return ans

def polyfitmag(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] + c[1]*(x-6)**2
    
    return ans
    
def fit_quadratic_vertex(c, x):
    
    xx = x - 8.
    
    return c[0] * xx**2 + c[1]   

# bi-linear function for area and width!
def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x >= hx
    xmod[idx] = 1
    return xmod
    
def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    from numpy import zeros_like
    hx = 0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)
    
    return ans1 + ans2

hxfix = log10(xref_hinge) #4.0 # hinge distance
def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    #hxfix = log10(800) #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    modx_lo = lowside(x, hxfix)
    modx_hi = highside(x, hxfix)

    ans1 = modx_lo * (c[0] * x + c[1])
    yhinge = c[0] * hxfix + c[1]
    ans2 = modx_hi * (c[2] * (x-hxfix) + yhinge)

    return ans1 + ans2
    
# normalise data to get atten pattern
def normalise_data(stdict, T):
    
    print('Regressing T =', str(T))
    init_mw = []
    norm_amp_all = []
    norm_rhyp = []
    norm_dep = []
    norm_mag = []
    norm_stas = []
    norm_date = []
    
    for i, mplt in enumerate(mrng):
                
        # for each record, log interpolate across period
        amp_plt = []
        addData = True # just for PGA & PGV
        for st in stdict:
            if T == 0.01:
                amp_plt.append(st['pga'])
            elif T == -99:
                amp_plt.append(st['pgv'])
            else:
                # add logic to remove short T data for long T
                addData = False
                if T > 4.0:
                    if st['chstr'][-1].startswith('H') or st['chstr'][-1].startswith('B'):
                        addData = True
                else:
                    addData = True
        
            if addData == True:
                amp_plt.append(exp(interp(log(T), log(st['per']), log(st['geom']))))
                idx = where((st['geom'] > 1E5) | (st['geom'] < 1E-10))[0]
                if len(idx) > 0:
                    print(' '.join(('BAD DATA', st['sta'], st['ev'])))
            else:
                amp_plt.append(nan)
        
        amp_plt = array(amp_plt)
        
        # get data
        rhyp = dictlist2array(stdict, 'rhyp')
        mags = dictlist2array(stdict, 'mag')
        deps = dictlist2array(stdict, 'dep')
        azim = dictlist2array(stdict, 'azim')
        auth = dictlist2array(stdict, 'network')
        stas = dictlist2array(stdict, 'sta')
        date = dictlist2array(stdict, 'date')
        
        # get events within mag and T bin
        #midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (deps >= 30.))[0]
        midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)))[0]
        
        if len(midx) > 0:
    
            # get binned medians
            bins = arange(log10(minDist), log10(maxDist), 0.1)
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[midx]), log(amp_plt[midx]))

            # normalise at XREF km
            bidx = where(around(binstrp, 1) == around(log10(xref),1))[0]
            
            if len(bidx) == 1:
                norm_amp = amp_plt[midx] / exp(logmedamp[bidx])
                idx = where((norm_amp < 1E-6) | (norm_amp > 1E6))[0]
                if len(idx) > 0:
                    for ix in idx:
                        print(' '.join(('BAD DATA', stas[ix], date[ix])))
                norm_amp_all = concatenate((norm_amp_all, norm_amp))
                norm_rhyp = concatenate((norm_rhyp, rhyp[midx]))
                norm_dep = concatenate((norm_dep, deps[midx]))
                norm_mag = concatenate((norm_mag, mags[midx]))
                norm_stas = concatenate((norm_stas, stas[midx]))
                norm_date = concatenate((norm_date, date[midx]))
                
    return norm_rhyp, norm_amp_all, norm_dep, norm_stas, norm_mag, norm_date

################################################################################
# now get geometric atten for all mags
################################################################################
def regress_zone(stdict, zgroup):
    
    Tplt = array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.]) # secs; PGA = 0.01; PGV = -99
    Tplt = array([0.067, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 9.5, 10.]) #, 10.]) # secs; PGA = 0.01; PGV = -99
    bins = arange(log10(minDist), log10(maxDist), 0.1)
    
    # compute inslab gmpes
    rrup = logspace(log10(minDist), log10(maxDist))

    #fig = plt.figure(1, figsize=(22,18))
    ii = 0    
    fig = plt.figure(1, figsize=(18,10))
    init_c0 = []
    init_c2 = []
    init_c3 = []
    bl_init_c0 = []
    bl_init_c1 = []
    bl_init_c2 = []
    bl_init_c3 = []
    blf_init_c0 = []
    blf_init_c1 = []
    blf_init_c2 = []

    for ii, T in enumerate(Tplt):
    
        norm_rhyp, norm_amp_all, norm_dep, norm_stas, norm_mag, norm_date = normalise_data(stdict, T)
        
        ################################################################################
        # now get atten for normalised data
        ################################################################################
         
        ax = plt.subplot(4,5,ii+1)
        idx = where(norm_amp_all < 1000)[0]
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyp[idx]), log(norm_amp_all[idx]))
        
        norm = mpl.colors.Normalize(vmin=0, vmax=500)
        
        #if ii < 20:
        #print(norm_amp_all)
        sc = plt.scatter(norm_rhyp, norm_amp_all, c=norm_dep, marker='o', s=25, \
                         cmap='Spectral_r', norm=norm, alpha=0.6)
                         
        # find suspect data
        norm_amp_all=array(norm_amp_all)
        index = where((norm_amp_all > 10) & (isnan(norm_amp_all) == False))[0]
        
        if len(index) > 0:
            for idx in index:
                print('Suspect:', norm_rhyp[idx], norm_mag[idx], norm_dep[idx], norm_stas[idx])
        
        #plt.plot(norm_rhyp, norm_amp_all,'bo')
        plt.loglog(10**medx, exp(logmedamp), 'rs', ms=6.5)
        
        ################################################################################
        # now fit normalised data - use binned data > hinge for far field
        ################################################################################
        # use binned data > hinge for far field
        hxfix = log10(xref_hinge)
        ridx = where(10**medx >= 10**hxfix)[0]
        data = odrpack.RealData(10**medx[ridx], logmedamp[ridx])
         
        ''' GR only '''
        afit = odrpack.Model(fit_gr)
        odr = odrpack.ODR(data, afit, beta0=[1.0, 1.0])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        c = out.beta
        
        # now plt
        attenfit = c[0] + c[1]*log10(rrup)
        #plt.loglog(rrup, exp(attenfit), 'k-', lw=2)
        init_c0.append(c[0])
        init_c2.append(c[1]) # this is now far-field slope!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    ################################################################################
    # fit bi-linear fixed hinge & slope
    ################################################################################
    # smooth bi-linear fixed coef
    sg_window = 7
    sg_poly = 2
    smooth_c2 = savitzky_golay(array(init_c2), sg_window, sg_poly) # slope
    
    for ii, T in enumerate(Tplt):
        ax = plt.subplot(4,5,ii+1)
        
        # prep data again
        norm_rhyp, norm_amp_all, norm_dep, norm_stas, norm_mag, norm_date = normalise_data(stdict, T)
        idx = where(norm_amp_all > 1E-6)[0]
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyp[idx]), log(norm_amp_all[idx]))
        data = odrpack.RealData(medx, logmedamp)
        
        # fix far-field slope
        sff_gr = smooth_c2[ii] # using smoothed far-field slope
        
        hxfix = log10(xref_hinge)
        def bilinear_reg_fix_hinge_slope(c, x):
            from numpy import zeros_like
            
            ans2 = zeros_like(x)
            ans1 = zeros_like(x)
        
            modx_lo = lowside(x, hxfix)
            modx_hi = highside(x, hxfix)
        
            ans1 = modx_lo * (c[0] * x + c[1])
            yhinge = c[0] * hxfix + c[1]
            ans2 = modx_hi * (sff_gr * (x-hxfix) + yhinge)
        
            return ans1 + ans2
            
        bilin_fix_hinge_slope = odrpack.Model(bilinear_reg_fix_hinge_slope)
        odr = odrpack.ODR(data, bilin_fix_hinge_slope, beta0=[-5, 1.])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
        out = odr.run()
        
        af = out.beta[0]
        bf = out.beta[1]
        
        if af > 0:
            af = 0
            bf = 1. 
        
        blf_init_c1.append(af)
        blf_init_c0.append(bf)
        blf_init_c2.append(sff_gr)
        
        attenfit_blf = bf + af * log10(rrup)
        yhinge = bf + af * hxfix
        idx = log10(rrup) > hxfix
        attenfit_blf[idx] = sff_gr * (log10(rrup[idx])-hxfix) + yhinge
        plt.loglog(rrup, exp(attenfit_blf), 'k-', lw=2)  
        
        ################################################################################
        # make pretty
        ################################################################################
        
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.xlim([50, maxDist])
        plt.ylim([1E-3, 1E2])
        plt.grid(which='both')
        #plt.title('Normalised Amplitudes')
        
        plt.ylabel('SA('+str(T)+') in g')
        if ii > 15:
            plt.xlabel('Hypocental Distance (km)')
        
    '''
    if zgroup == 'NGH':
        bl_init_c0 = array(init_c0)
        bl_init_c1 = array(init_c2) # c2 from above
        bl_init_c2 = array(init_c2)*0.
        hxfix = 4 # no hinge
    '''
    if zgroup == 'BS' or zgroup == 'NGH' or zgroup == 'OBE':
        bl_init_c0 = array(blf_init_c0)
        bl_init_c1 = array(blf_init_c1)
        bl_init_c2 = array(blf_init_c2)
    else:
        bl_init_c0 = array(bl_init_c0)
        bl_init_c1 = array(bl_init_c1)
        bl_init_c2 = array(bl_init_c2)
    
    # smooth bi-linear - near-field slope
    sg_window = 7
    sg_poly = 2
    #smooth_c1 = savitzky_golay(array(bl_init_c1), sg_window, sg_poly) # far-field slope
            
    plt.suptitle('G(R) and Q')
    plt.savefig('.'.join(('gr',zgroup,'png')), fmt='png', bbox_inches='tight')
    plt.show()
            
    ################################################################################
    # plot G(R) coeffs
    ################################################################################
    
    fig = plt.figure(2, figsize=(12, 9))
    
    # plt far-field slope
    plt.subplot(311)
    #plt.semilogx(Tplt, init_c2, 'bo')
    
    # add BL coeffs
    plt.semilogx(array(Tplt), -1*array(init_c2), 'go', label='raw')
    
    plt.semilogx(Tplt, -1*smooth_c2, 'ro', label='smoothed')
    plt.ylabel('far-field c2 (GR)')
    plt.legend(loc=3)
    
    # plt near-field slope
    plt.subplot(312)
    plt.semilogx(array(Tplt), -1*array(bl_init_c1), 'bo', label='raw')
    
    # smooth NF-slope
    smooth_c1 = savitzky_golay(bl_init_c1, sg_window, sg_poly)
    plt.semilogx(Tplt, -1*array(smooth_c1), 'ro', label='smoothed')
    
    plt.ylabel('near-field slope c1')
    plt.xlabel("T (sec)")
    plt.legend()
    
    plt.show()

    ################################################################################
    # loop through T and M to get mag scaling
    ################################################################################
    print('Fitting magnitude dependence...')
    fig = plt.figure(3, figsize=(12,6))
    
    tidx = 2
    #init_mw = zeros((len(Tplt), len(mrng)))
    ii = 0
    mag_per_fit = []
    mag_c0_fit = []
    m0_array = []
    m1_array = []
    m2_array = []
    
    t_dept_c0 = []
    t_dept_mw = []
    
    t_amplitudes = []
    
    # get data
    rhyp = dictlist2array(stdict, 'rhyp')
    mags = dictlist2array(stdict, 'mag')
    deps = dictlist2array(stdict, 'dep')
    azim = dictlist2array(stdict, 'azim')
    stas = dictlist2array(stdict, 'sta')
    date = dictlist2array(stdict, 'date')
    '''
    print('\nMW 6.3')
    wmag = 6.3
    midx = where((mags >= (wmag-mpltrng)) & (mags < (wmag+mpltrng)))[0]
    print(stas[midx])
    print(date[midx])
    
    print('\nMW 7.5')
    wmag = 7.5
    midx = where((mags >= (wmag-mpltrng)) & (mags < (wmag+mpltrng)))[0]
    print(stas[midx])
    print(date[midx])
    '''
    for tt, T in enumerate(Tplt):
        init_mw = []
        init_c0 = []
        
        # for each record, log interpolate across period
        amp_plt = []
        for st in stdict:
            if T == 0.01:
                amp_plt.append(st['pga'])
            elif T == -99:
                amp_plt.append(st['pgv'])
            else:
                amp_plt.append(exp(interp(log(T), log(st['per']), log(st['geom']))))
        
        amp_plt = array(amp_plt)
        t_amplitudes.append(amp_plt)
        
        for i, mplt in enumerate(mrng):
            
            # for banda sea - just use h > 40 km to calibrate mag scaling
            if zgroup == 'BS':
                midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (isnan(amp_plt) == False) & (deps > 40.))[0]
            else:
                midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (isnan(amp_plt) == False) & (deps < 300.))[0]
            
            if len(midx) > 0:
                # get binned medians
                bins = arange(log10(minDist), log10(maxDist), 0.1)
        
                logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[midx]), log(amp_plt[midx]))
                
                ii += 1
                
                # fit BL atten intercept for mag scaling
                c1 = smooth_c1[tt]
                c2 = smooth_c2[tt]
                def fit_bl_mag_scaling(c, x):
                    from numpy import log10, zeros_like
                    
                    ans2 = zeros_like(x)
                    ans1 = zeros_like(x)
                    
                    modx_lo = lowside(log10(x), hxfix)
                    modx_hi = highside(log10(x), hxfix)
                    
                    ans1 = modx_lo * (c1 * log10(x) + c[0])
                    yhinge = c1 * hxfix + c[0]
                    ans2 = modx_hi * (c2 * (log10(x)-hxfix) + yhinge)
                    
                    return ans1 + ans2
                    
                data = odrpack.RealData(rhyp[midx], log(amp_plt[midx]))
                
                afit = odrpack.Model(fit_bl_mag_scaling)
                odr = odrpack.ODR(data, afit, beta0=[1.0])
                odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
                out = odr.run()
                c = out.beta
                
                # add to array
                init_c0.append(c[0])
                init_mw.append(mplt)
        
        # fit straight quadratic
        nn = where(isnan(init_c0)==False)[0]
        m1, m2, m0 = polyfit((array(init_mw)[nn]-6), array(init_c0)[nn], 2)
        m0_array.append(m0)
        m1_array.append(m1)
        m2_array.append(m2)
        
        t_dept_c0.append(array(init_c0))
        t_dept_mw.append(array(init_mw))
        
    m0_array = array(m0_array)
    m1_array = array(m1_array)
    m2_array = array(m2_array)
    
    # plot polynomial coefs
    plt.subplot(223)
    plt.semilogx(Tplt, m0_array, 'bo',label='m0')
    plt.semilogx(Tplt, m1_array, 'go',label='m1')
    plt.semilogx(Tplt, m2_array, 'ro',label='m2')
    plt.legend()
                
    ######################################################################################
    # plot M1
    plt.subplot(221)
    plt.semilogx(Tplt, array(m0_array), 'bo', label='raw')
    smooth_m0 = savitzky_golay(array(m0_array), sg_window, sg_poly)
    plt.semilogx(Tplt, smooth_m0, 'ro', label='smoothed')
    plt.title('m0')
    plt.legend(loc=0)
    
    # smoothing of M1 is terrible, so use orig vals!
    smooth_m0 = m0_array
    smooth_m1 = m1_array
    smooth_m2 = m2_array
    """
    ######################################################################################
    # refit M1 with smoothed M!
    refit_m1 = []
    refit_m2 = []

    for c0_dat, mw_dat, sm0 in zip(t_dept_c0, t_dept_mw, smooth_m0):
        # check nan values
        nn = where(isnan(c0_dat)==False)[0]  
        
        def fit_fixed_polynomial1(c, x):
            return sm0 + c[0]*(mw_dat[nn]-6.)**2 + c[1]*(mw_dat[nn]-6.)
        
        data = odrpack.RealData(mw_dat[nn], c0_dat[nn])
    
        magfit = odrpack.Model(fit_fixed_polynomial1)
        odr = odrpack.ODR(data, magfit, beta0=[-0.1, 1.])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        m = out.beta
        
        refit_m1.append(m[0])
        refit_m2.append(m[1])
    
    smooth_m1 = savitzky_golay(array(refit_m1), sg_window, sg_poly)
    """
    # plot M2            
    plt.subplot(222)
    plt.semilogx(Tplt, array(m1_array), 'bo', label='orig')
    #plt.semilogx(Tplt, array(refit_m1), 'go', label='refit')
    plt.semilogx(Tplt, smooth_m1, 'ro', label='smooth')
    plt.title('m1')
    plt.legend(loc=0)
    
    ######################################################################################
    # refit M2 with smoothed M!
    refit_m2 = []
    """
    for c0_dat, mw_dat, sm0, sm1 in zip(t_dept_c0, t_dept_mw, smooth_m0, smooth_m1):
        # check nan values
        nn = where(isnan(c0_dat)==False)[0]  
        
        def fit_fixed_polynomial2(c, x):
            
            return sm0 + sm1*(mw_dat[nn]-6.)**2 + c[0]*(mw_dat[nn]-6.)
    
        data = odrpack.RealData(array(mw_dat[nn]), array(c0_dat[nn]))
    
        magfit = odrpack.Model(fit_fixed_polynomial2)
        odr = odrpack.ODR(data, magfit, beta0=[1.])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        m = out.beta
        
        refit_m2.append(m[0])
    
    #smooth_m2 = savitzky_golay(array(refit_m2), sg_window, sg_poly)
    smooth_m2 = refit_m2
    """
    # plot M2
    plt.subplot(224)
    plt.semilogx(Tplt, array(m2_array), 'bo', label='orig')
    #plt.semilogx(Tplt, refit_m2, 'go', label='refit')
    plt.semilogx(Tplt, smooth_m2, 'ro', label='smoothed')
    plt.title('m2')
    plt.legend(loc=0)
    
    plt.suptitle('Mag Coefs')
    plt.show()
    
    ######################################################################################
    # loop throught periods and plot mag scaling
    ######################################################################################
    fig = plt.figure(10, figsize=(18,10))
    
    i = 0
    for T, ym, xm in zip(Tplt, t_dept_c0, t_dept_mw):
        plt.subplot(4,5,i+1)
        plt.plot(xm, ym, 'bo')
        plt.ylabel(str(T))
        
        #mfit = m0_array[i] + m1_array[i]*(mrng-6)**2 + m2_array[i]*(mrng-6)
        mfit = smooth_m0[i] + smooth_m1[i]*(mrng-6)**2 + smooth_m2[i]*(mrng-6)
        plt.plot(mrng, mfit, 'g-', lw=2.)
        #plt.ylim([4.5, 7])
        
        if i >= 16:
            plt.xlabel('MW')
        
        i += 1
    
    plt.suptitle('mag scaling')
    plt.savefig('.'.join(('mag_scaling', zgroup, 'png')), fmt='png', bbox_inches='tight')
    plt.show()
    #blah=blah
    ######################################################################################
    # loop through periods and plot 1st distance residual
    ######################################################################################
    fig = plt.figure(30, figsize=(18,10))
    
    for i, T in enumerate(Tplt):
        #attenCor = (- smooth_c2[i]*log10(rhyp) - smooth_c3[i]*(rhyp)) # old
        atten_cor = smooth_c1[i] * log10(rhyp)
  
        yhinge = smooth_c1[i] * hxfix
        idx = log10(rhyp) > hxfix
        atten_cor[idx] = smooth_c2[i] * (log10(rhyp[idx])-hxfix) + yhinge
                    
        logY = smooth_m0[i] + smooth_m1[i]*(mags-6)**2 + smooth_m2[i]*(mags-6) \
               + atten_cor
              
        resY = log(t_amplitudes[i]) - logY
    
        # plot
        ax = plt.subplot(4,5,i+1)
        
        norm = mpl.colors.Normalize(vmin=0, vmax=700)
        
        sc = plt.scatter(log10(rhyp), resY, c=deps, marker='o', s=20, \
                         cmap='viridis_r', norm=norm, alpha=1.0)
        plt.plot([1,3.4], [0,0], 'r--', lw=1.5)
        
        # get binned medians
        bins = arange(2.0, log10(maxDist), 0.1)
        nn = where((isnan(resY) == False) & (isnan(log10(rhyp)) == False) & (isfinite(resY) == True))[0]
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[nn]), resY[nn])
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
        
        plt.xlim([2,3.4])
        plt.ylim([-6, 6])
        plt.ylabel(str(T))
        
        if i >= 16:
            plt.xlabel('log Rhyp (km)')
    
    plt.suptitle('First Dist Residuals')
    plt.savefig('init_dist_res.png', fmt='png', bbox_inches='tight')
    plt.show()
    
    ######################################################################################
    # calculate model residual and plot with depth
    ######################################################################################
    print('Fitting depth dependence...')
    #log Y = c0 + c1*(M-6)**2 + c2*(M-6) - c3*log Rhyp - c4*Rhyp
    fig = plt.figure(11, figsize=(18,10))
    mod_res = []
    log_drng = arange(1, 2.9, 0.1)
    
    # set coef array
    d0_array = []
    d1_array = []
    d2_array = []
    d3_array = []
    t_dept_res = []
    t_dept_logdep = []
    
    for i, T in enumerate(Tplt):
        #attenCor = (- smooth_c2[i]*log10(rhyp) - smooth_c3[i]*(rhyp)) # old
        atten_cor = smooth_c1[i] * log10(rhyp)
  
        yhinge = smooth_c1[i] * hxfix
        idx = log10(rhyp) > hxfix
        atten_cor[idx] = smooth_c2[i] * (log10(rhyp[idx])-hxfix) + yhinge
                    
        logY = smooth_m0[i] + smooth_m1[i]*(mags-6)**2 + smooth_m2[i]*(mags-6) \
               + atten_cor
              
        resY = log(t_amplitudes[i]) - logY
         
        # plot
        ax = plt.subplot(4,5,i+1)
         
        norm = mpl.colors.Normalize(vmin=0, vmax=2000)
        
        sc = plt.scatter(log10(deps), resY, c=rhyp, marker='o', s=20, \
                         cmap='viridis_r', norm=norm, alpha=1.0)
        plt.plot([0.4,3], [0,0], 'r--', lw=1.5)
        
        plt.ylabel(str(T))
        if i >= 16:
            plt.xlabel('log Depth (km)')
        plt.xlim([0.4, 3])
        plt.ylim([-6, 6])
        
        # get binned medians
        bins = arange(0.5, 2.8, 0.1)
        
        if zgroup == 'OBE':
            nn = where((isnan(resY) == False) & (isnan(log10(deps)) == False) & (isfinite(resY) == True) & (deps < 300.))[0]
        else:
            nn = where((isnan(resY) == False) & (isnan(log10(deps)) == False) & (isfinite(resY) == True))[0]
        
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(deps[nn]), resY[nn])
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
        
        if zgroup == 'NGH' or zgroup == 'OBE': # fit linear
            d3, d0 = polyfit(log10(deps[nn]), resY[nn], 1) # full data
            #d3, d0 = polyfit(array(medx), array(logmedamp), 1) # binned data
            d0_array.append(d0)
            d1_array.append(0.0)
            d2_array.append(0.0)
            d3_array.append(d3)
            d1 = 0.
            d2 = 0.
            
        else:
            # fit straight cubic
            #d1, d2, d3, d0 = polyfit(array(medx), array(logmedamp), 3) # binned data
            d1, d2, d3, d0 = polyfit(log10(deps[nn]), resY[nn], 3) # full data
            d0_array.append(d0)
            d1_array.append(d1)
            d2_array.append(d2)
            d3_array.append(d3)
        
        '''
        # binned data
        t_dept_res.append(logmedamp)
        t_dept_logdep.append(medx)
        '''
        # full data
        t_dept_res.append(resY[nn])
        t_dept_logdep.append(log10(deps[nn]))
        
        # plot quadratics
        dx = arange(0.5, 2.9, 0.01)
        dy = d0 + d1*dx**3 + d2*dx**2 + d3*dx
        #print(dx, dy)
        plt.plot(dx, dy, 'g-', lw=2.)
        
    d0_array = array(d0_array)
    d1_array = array(d1_array)
    d2_array = array(d2_array)
    d3_array = array(d3_array)
    
    plt.suptitle('Depth Residuals')
    plt.savefig('.'.join(('depth_residuals',zgroup,'png')), fmt='png', bbox_inches='tight')
    plt.show()
    
    # smooth coefs
    smooth_d0 = savitzky_golay(d0_array, sg_window, sg_poly)
    
    # plot depth coeffs
    fig = plt.figure(13, figsize=(18,10))
    # plot D0
    plt.subplot(221)
    plt.semilogx(Tplt, array(d0_array), 'bo')
    plt.semilogx(Tplt, smooth_d0, 'ro')
    plt.title('d0')
    #print(smooth_d0)
    
    plt.subplot(222)
    plt.semilogx(Tplt, array(d1_array), 'bo')
    #plt.semilogx(dx, smooth_d0, 'ro')
    plt.title('d1')
    
    plt.subplot(223)
    plt.semilogx(Tplt, array(d2_array), 'bo')
    #plt.semilogx(dx, smooth_d0, 'ro')
    plt.title('d2')
    
    plt.subplot(224)
    plt.semilogx(Tplt, array(d3_array), 'bo')
    #plt.semilogx(dx, smooth_d0, 'ro')
    plt.title('d3')
    
    ######################################################################################
    # refit D1 with smoothed DO
    ######################################################################################
    
    if zgroup == 'BS':
        refit_d1 = []
        refit_d2 = []
        refit_d3 = []
        # m0_array[i] + d1_array[i]*(mrng-6)**2 + d2_array[i]*(mrng-6)
        #for c0_dat, sd1 in zip(mag_c0_fit, smooth_m0): # t_dept_c0, t_dept_mw
        for t_res, dep_dat, sd0 in zip(t_dept_res, t_dept_logdep, smooth_d0):
             
            # check nan values
            nn = where(isnan(t_res)==False)[0]  
            
            def fit_fixed_polynomial1(c, x):
                return sd0 + c[0]*dep_dat[nn]**3 + c[1]*dep_dat[nn]**2 + c[2]*dep_dat[nn]
            
            data = odrpack.RealData(dep_dat[nn], t_res[nn])
            #data = odrpack.RealData(dep_dat, t_res)
        
            depfit = odrpack.Model(fit_fixed_polynomial1)
            odr = odrpack.ODR(data, depfit, beta0=[0.1, 1., 1.])
            
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            m = out.beta
            
            refit_d1.append(m[0])
            refit_d2.append(m[1])
            refit_d3.append(m[2])
        
        smooth_d1 = savitzky_golay(array(refit_d1), sg_window, sg_poly)
        
        ######################################################################################
        # refit d2 with smoothed d1!
        refit_d2 = []
        refit_d3 = []
        # m0_array[i] + d1_array[i]*(mrng-6)**2 + d2_array[i]*(mrng-6)
        #for c0_dat, sd1 in zip(mag_c0_fit, smooth_m0): # t_dept_c0, t_dept_mw
        for t_res, dep_dat, sd0, sd1 in zip(t_dept_res, t_dept_logdep, smooth_d0, smooth_d1):
             
            # check nan values
            nn = where(isnan(t_res)==False)[0]  
            
            def fit_fixed_polynomial1(c, x):
                return sd0 + sd1*dep_dat[nn]**3 + c[0]*dep_dat[nn]**2 + c[1]*dep_dat[nn]
            
            data = odrpack.RealData(dep_dat[nn], t_res[nn])
        
            depfit = odrpack.Model(fit_fixed_polynomial1)
            odr = odrpack.ODR(data, depfit, beta0=[1., 1.])
            
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            m = out.beta
            
            refit_d2.append(m[0])
            refit_d3.append(m[1])
        
        smooth_d2 = savitzky_golay(array(refit_d2), sg_window, sg_poly)
        
        ######################################################################################
        # refit d3 with smoothed M!
        refit_d3 = []
        # m0_array[i] + d1_array[i]*(mrng-6)**2 + d2_array[i]*(mrng-6)
        #for c0_dat, sd1 in zip(mag_c0_fit, smooth_m0): # t_dept_c0, t_dept_mw
        for t_res, dep_dat, sd0, sd1, sd2 in zip(t_dept_res, t_dept_logdep, smooth_d0, smooth_d1, smooth_d2):
             
            # check nan values
            nn = where(isnan(t_res)==False)[0]  
            
            def fit_fixed_polynomial1(c, x):
                return sd0 + sd1*dep_dat[nn]**3 + sd2*dep_dat[nn]**2 + c[0]*dep_dat[nn]
            
            data = odrpack.RealData(dep_dat[nn], t_res[nn])
            #data = odrpack.RealData(dep_dat, t_res)
        
            depfit = odrpack.Model(fit_fixed_polynomial1)
            odr = odrpack.ODR(data, depfit, beta0=[1.])
            
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            m = out.beta
            
            refit_d3.append(m[0])
        
        #smooth_d3 = savitzky_golay(array(refit_d3), sg_window, sg_poly)
        smooth_d3 = refit_d3
    
    else:
        smooth_d1 = d0_array * 0.
        smooth_d2 = d0_array * 0.
        
        ######################################################################################
        # refit d3 with smoothed M!
        refit_d3 = []
        for t_res, dep_dat, sd0 in zip(t_dept_res, t_dept_logdep, smooth_d0):
             
            # check nan values
            nn = where(isnan(t_res)==False)[0]  
            
            def fit_fixed_polynomial1(c, x):
                return sd0 + c[0]*dep_dat[nn]
            
            data = odrpack.RealData(dep_dat[nn], t_res[nn])
            #data = odrpack.RealData(dep_dat, t_res)
        
            depfit = odrpack.Model(fit_fixed_polynomial1)
            odr = odrpack.ODR(data, depfit, beta0=[1.])
            
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            m = out.beta
            
            refit_d3.append(m[0])
        
        smooth_d3 = savitzky_golay(array(refit_d3), sg_window, sg_poly)
    
    # over ride smoothing
    smooth_d0 = d0_array
    smooth_d1 = d1_array
    smooth_d2 = d2_array
    smooth_d3 = d3_array
    
    
    # plt smoothed params
    plt.subplot(222)
    plt.semilogx(Tplt, smooth_d1, 'ro')
    
    plt.subplot(223)
    plt.semilogx(Tplt, smooth_d2, 'ro')
    
    plt.subplot(224)
    plt.semilogx(Tplt, refit_d3, 'go')
    plt.semilogx(Tplt, smooth_d3, 'ro')
    
    plt.suptitle('Depth Coefs')
    plt.show()
    
    ######################################################################################
    # calculate model residual and plot with distance again
    ######################################################################################
    print('Fitting 2nd distance dependence...')
    #log Y = c0 + c1*(M-6)**2 + c2*(M-6) - c3*log Rhyp - c4*Rhyp
    fig = plt.figure(14, figsize=(18,10))
    mod_res = []
    log_drng = arange(1, 2.9, 0.1)
    
    for i, T in enumerate(Tplt):
        #atten_cor = (- smooth_c2[i]*log10(rhyp) - smooth_c3[i]*(rhyp)) # old
        atten_cor = smooth_c1[i] * log10(rhyp)
        yhinge = smooth_c1[i] * hxfix
        idx = log10(rhyp) > hxfix
        atten_cor[idx] = smooth_c2[i] * (log10(rhyp[idx])-hxfix) + yhinge
        
        dep_cor = smooth_d0[i] + smooth_d1[i]*log10(deps)**3 + smooth_d2[i]*log10(deps)**2 + smooth_d3[i]*log10(deps)
        
        logY = smooth_m0[i] + smooth_m1[i]*(mags-6)**2 + smooth_m2[i]*(mags-6) \
               + atten_cor + dep_cor
              
        resY = log(t_amplitudes[i]) - logY
         
        # plot
        ax = plt.subplot(4,5,i+1)
        
        norm = mpl.colors.Normalize(vmin=0, vmax=700)
        
        sc = plt.scatter(log10(rhyp), resY, c=deps, marker='o', s=20, \
                         cmap='viridis_r', norm=norm, alpha=1.0)
        
        plt.plot([1,3.4], [0,0], 'r--', lw=1.5)
        
        # get binned medians
        bins = arange(log10(minDist), log10(maxDist), 0.1)
        nn = where((isnan(resY) == False) & (isnan(log10(rhyp)) == False) & (isfinite(resY) == True))[0]
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[nn]), resY[nn])
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
        
        plt.xlim([2,3.4])
        plt.ylim([-4, 4])
        plt.ylabel(str(T))
        
        if i >= 16:
            plt.xlabel('log Rhyp (km)')
    
        sidx = where((resY < -10) & (log10(rhyp)<2.6))[0]
        #sidx = where(resY .)[0]
        print(stas[sidx])
        print(unique(array(stas)))
        #print(azim[sidx])
        print(date[sidx])
    plt.suptitle('Dist Residuals')
    plt.savefig('final_res.'+region+'.png', fmt='png', bbox_inches='tight')
    plt.show()
    
    ######################################################################################
    # calculate model residual and plot with mag
    ######################################################################################
    print('Plotting mag residuals...')
    #log Y = c0 + c1*(M-6)**2 + c2*(M-6) - c3*log Rhyp - c4*Rhyp
    fig = plt.figure(15, figsize=(18,10))
    mod_res = []
    bins = arange(5.4, 8.1, 0.2)
    
    for i, T in enumerate(Tplt):
        #atten_cor = (- smooth_c2[i]*log10(rhyp) - smooth_c3[i]*(rhyp)) # old
        atten_cor = smooth_c1[i] * log10(rhyp)
        yhinge = smooth_c1[i] * hxfix
        idx = log10(rhyp) > hxfix
        atten_cor[idx] = smooth_c2[i] * (log10(rhyp[idx])-hxfix) + yhinge
        
        dep_cor = smooth_d0[i] + smooth_d1[i]*log10(deps)**3 + smooth_d2[i]*log10(deps)**2 + smooth_d3[i]*log10(deps)
        
        logY = smooth_m0[i] + smooth_m1[i]*(mags-6)**2 + smooth_m2[i]*(mags-6) \
               + atten_cor + dep_cor
              
        resY = log(t_amplitudes[i]) - logY
         
        # plot
        ax = plt.subplot(4,5,i+1)
        
        norm = mpl.colors.Normalize(vmin=0, vmax=700)
        
        sc = plt.scatter(mags, resY, c=deps, marker='o', s=20, \
                         cmap='viridis_r', norm=norm, alpha=1.0)
        
        plt.plot([5,8], [0,0], 'r--', lw=1.5)
        
        # get binned medians
        #bins = arange(log10(minDist), log10(maxDist), 0.1)
        nn = where((isnan(resY) == False) & (isnan(mags) == False) & (isfinite(resY) == True))[0]
        logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, mags[nn], resY[nn])
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
        
        plt.xlim([5.,8])
        plt.ylim([-4, 4])
        plt.ylabel(str(T))
        
        if i >= 16:
            plt.xlabel('MW')
    
        sidx = where((resY < -3) & (log10(rhyp)<2.6))[0]
        #sidx = where(resY .)[0]
        print(medx)
        #print(azim[sidx])
        print(nperbin)
    plt.suptitle('Mag Residuals')
    plt.savefig('mag_res.'+region+'.png', fmt='png', bbox_inches='tight')
    plt.show()
    #################################################################################
    # export coeffs
    #################################################################################
    
    #ctxt = '#ln Y = c0 + c1*(M-6)**2 + c2*(M-6) - c3*log10(Rhyp) - c4*(Rhyp) + (d0 + d1*log10(h)**3 + d2*log10(h)**2 + d3*log10(h))\nT, c0, c1, c2, c3, c4, d0, d1, d2, d3\n'
    #ctxt = '#Rhyp <= hx: ln Y = c0 + c1*(M-6)**2 + c2*(M-6) + (c3*log10(Rhyp)) + (d0 + d1*log10(h)**3 + d2*log10(h)**2 + d3*log10(h))\n'
    ctxt = '#Rhyp <= hx: ln Y = c0 + c1*(M-6)**2 + c2*(M-6) + (c3*hx +  c4*(log10(Rhyp)-hx)) + (d0 + d1*log10(h)**3 + d2*log10(h)**2 + d3*log10(h))\nT, c0, c1, c2, c3, c4, d0, d1, d2, d3, hx\n'
    for i, t in enumerate(Tplt):
        '''
        ctxt += ','.join((str(t), str('%0.5f' % smooth_m0[i]), str('%0.5f' % smooth_m1[i]), str('%0.5f' % smooth_m2[i]), \
                          str('%0.5f' % smooth_c2[i]), str('%0.5e' % smooth_c3[i]), \
                          str('%0.5f' % smooth_d0[i]), str('%0.5f' % smooth_d1[i]), \
                          str('%0.5f' % smooth_d2[i]), str('%0.5f' % smooth_d3[i]))) + '\n'
        '''
        ctxt += ','.join((str(t), str('%0.5f' % smooth_m0[i]), str('%0.5f' % smooth_m1[i]), str('%0.5f' % smooth_m2[i]), \
                          str('%0.5f' % smooth_c1[i]), str('%0.5e' % smooth_c2[i]), \
                          str('%0.5f' % smooth_d0[i]), str('%0.5f' % smooth_d1[i]), \
                          str('%0.5f' % smooth_d2[i]), str('%0.5f' % smooth_d3[i]), \
                          str('%0.5f' % hxfix))) + '\n'
                          
        '''
        atten_cor = bl_init_c1[i] * log10(rhyp)
        yhinge = bl_init_c1[i] * hxfix
        idx = log10(rhyp) > hxfix
        atten_cor[idx] = bl_init_c2[i] * (log10(rhyp[idx])-hxfix) + yhinge
        '''
    f = open('.'.join(('ncc_gmm_coeffs', zgroup, 'csv')), 'w')
    f.write(ctxt)
    f.close()
    
    return resY, deps, norm_stas

################################################################################
# parse shapefile and filter stdict, and regress
################################################################################
'''
BS = Banda Sea
NGH = New Guinea Highlands
OBW = Oceanic Basin-Extended Margin West
OBE = Oceanic Basin-Extended Margin East
'''
# load shape
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')
zone_group = get_field_data(sf, 'ZONE_GROUP', 'str')

# loop thru zones
i = 0
reg_stdict = []

if region == 'NGH':

    zgroup1 = 'NGH'
    zgroup2 = 'NGH'

elif region == 'OBE':
    zgroup1 = 'OBE'
    zgroup2 = 'OBW'

elif region == 'BS':
    zgroup1 = 'BS'
    zgroup2 = 'BS'

if zgroup1 == 'BS':
    mmin = 5.25
elif zgroup1 == 'OBE':
    mmin = 5.85
else:
    mmin = 5.75
print('Starting inversion...')
for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
    for sd in stdict:
        pt = Point(sd['eqlo'], sd['eqla'])
        if pt.within(poly) and sd['mag'] >= mmin:
            if zgroup == zgroup1 or zgroup == zgroup1:
                reg_stdict.append(sd)
            
# do regression for each polygon
resY, deps, norm_stas = regress_zone(reg_stdict, zgroup1)
