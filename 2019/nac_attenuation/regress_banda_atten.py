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
    datestr = lines[1].strip().split('\t')[-1]    
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
                  concatenate, around, mean, polyfit, delete, isnan
from misc_tools import listdir_extension, dictlist2array, get_binned_stats, savitzky_golay
#from calc_oq_gmpes import inslab_gsims
from os import path
from scipy.stats import linregress

################################################################################
# loop through sa files and get data
################################################################################
#usgscsv = '20190625_merged_events.csv'
#evdict = parse_usgs_events(usgscsv)

folder = 'psa'
extension = 'psa'
safiles = listdir_extension(folder, extension)
recs = []
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
        if sd['ev'] == saf[0:16] and sd['sta'] == fsplit[3]:
            if fsplit[4].endswith('E'):
                sd['chans'][0] = 1
                sd['chstr'][0] = fsplit[-2]
            elif fsplit[4].endswith('N'):
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
        if fsplit[4].endswith('E'):
            chans[0] = 1
            chstr[0] = fsplit[-2]
        elif fsplit[4].endswith('N'):
            chans[1] = 1
            chstr[1] = fsplit[-2]
        elif fsplit[4].endswith('Z'):
            chans[2] = 1
            chstr[2] = fsplit[-2]
        
        sd = {'ev':saf[0:16], 'sta':fsplit[3], 'chans':chans, 'chstr':chstr, 'net':fsplit[2]}
        
        stdict.append(sd)
            
################################################################################
# now get geometric mean for each stn
################################################################################

for i, st in enumerate(stdict):
    # get geometric mean of horizontals
    if st['chans'][0] == 1 and st['chans'][1] == 1:
        for rec in recs:
            # get east channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
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
                
            # get north channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][1]:
                
                ndat = rec['sa']
                
        stdict[i]['geom'] = exp((log(edat) + log(ndat)) / 2.)
    
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
                
didx = []
'''
for j, st in enumerate(stdict):
    if 'rhyp' in st:
        if st['net'] == 'GE' and st['rhyp'] > 400:
             didx.append(j)
    else:
        didx.append(j)
'''
stdict = delete(stdict, array(didx))
         	

################################################################################
# setup inversions
################################################################################
import scipy.odr.odrpack as odrpack

def fit_atten(c, x):
    from numpy import sqrt, log10
    
    xref = 100.
    ans = c[0] - c[1]*(x) - c[2]*log10(x)
    
    return ans
    
def fit_gr(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] - c[1]*log10(x)
    
    return ans
    


'''
pwr = 2.
def fit_atten(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] - c[1] * (x**pwr + c[2]**pwr)**(1./pwr)
    
    return ans
'''

################################################################################
# parse shapefile and filter stdict
################################################################################
'''
BS = Banda Sea
NGH = New Guinea Highlands
OB = Oceanic Basin-Extended Margin
'''
# load shape
import shapefile
from shapely import Point, Polygon
shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

# loop thru zones
i = 0

reg_stdict = []
for sd in stdict:
    for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
        pt = Point(sd['eqlo'], sd['eqla'])
        if zgroup == 'BS' and pt.within(poly):
            reg_stdict.append(sd)
    

################################################################################
# now get geometric atten for all mags
################################################################################
    
mrng = arange(5.3, 7.7, 0.1)
mpltrng = 0.05
Tplt = array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]) # secs; PGA = 0.01; PGV = -99

# compute inslab gmpes
rrup = logspace(2, log10(3000))
xref = 100.
#fig = plt.figure(1, figsize=(22,18))
ii = 0    
fig = plt.figure(1, figsize=(18,10))
init_c0 = []
init_c2 = []
init_c3 = []
for ii, T in enumerate(Tplt):
    
    init_mw = []
    norm_amp_all = []
    norm_rhyp = []
    norm_dep = []
    norm_mag = []
    
    for i, mplt in enumerate(mrng):
                
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
        
        # get data
        rhyp = dictlist2array(stdict, 'rhyp')
        mags = dictlist2array(stdict, 'mag')
        deps = dictlist2array(stdict, 'dep')
        azim = dictlist2array(stdict, 'azim')
        auth = dictlist2array(stdict, 'network')
        
        # get events within mag and T bin
        #midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (deps >= 30.))[0]
        midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)))[0]
        
        if len(midx) > 0:
    
            # get binned medians
            bins = arange(2.0, log10(3000), 0.1)
    
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[midx]), log10(amp_plt[midx]))
            
            # normalise at 1000 km
            bidx = where(around(binstrp, 1) == around(2.8,1))[0]
            if len(bidx) == 1:
                norm_amp = amp_plt[midx] / 10**logmedamp[bidx]
                norm_amp_all = concatenate((norm_amp_all, norm_amp))
                norm_rhyp = concatenate((norm_rhyp, rhyp[midx]))
                norm_dep = concatenate((norm_dep, deps[midx]))
                norm_mag = concatenate((norm_mag, mags[midx]))
    
    ################################################################################
    # now get atten for normalised data
    ################################################################################
    
    ax = plt.subplot(4,5,ii+1)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyp), log10(norm_amp_all))
    
    norm = mpl.colors.Normalize(vmin=0, vmax=500)
    
    sc = plt.scatter(norm_rhyp, norm_amp_all, c=norm_dep, marker='o', s=50, \
                     cmap='Spectral_r', norm=norm, alpha=0.6)
    
    #plt.plot(norm_rhyp, norm_amp_all,'bo')
    plt.loglog(10**medx, 10**logmedamp, 'rs', ms=6.5)
    
    ################################################################################
    # now fit normalised data
    ################################################################################
    #ridx = where(norm_rhyp <= 1000.)[0]
    #data = odrpack.RealData(norm_rhyp[ridx], log10(norm_amp_all[ridx]))
    
    # use binned data
    ridx = where(10**medx <= 1000.)[0]
    data = odrpack.RealData(10**medx[ridx], logmedamp[ridx])
    
    ''' GR + Q '''
    #afit = odrpack.Model(fit_atten)
    #odr = odrpack.ODR(data, afit, beta0=[1.0, 1.0, 1.0])
    
    ''' GR only '''
    afit = odrpack.Model(fit_gr)
    odr = odrpack.ODR(data, afit, beta0=[1.0, 1.0])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    c = out.beta
    
    # now plt
    #attenfit = c[0] - c[1]*rrup - c[2]*log10(rrup)
    attenfit = c[0] - c[1]*log10(rrup)
    plt.loglog(rrup, 10**attenfit, 'k-', lw=2)
    init_c0.append(c[0])
    init_c2.append(c[1])
    
################################################################################
# smooth GR and calc Q
################################################################################    
smooth_c2 = savitzky_golay(init_c2, 7, 3)

for ii, T in enumerate(Tplt):
	
    ax = plt.subplot(4,5,ii+1)
    
    init_mw = []
    norm_amp_all = []
    norm_rhyp = []
    norm_dep = []
    norm_mag = []
    
    for i, mplt in enumerate(mrng):
                
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
        
        # get data
        rhyp = dictlist2array(stdict, 'rhyp')
        mags = dictlist2array(stdict, 'mag')
        deps = dictlist2array(stdict, 'dep')
        azim = dictlist2array(stdict, 'azim')
        auth = dictlist2array(stdict, 'network')
        
        # get events within mag and T bin
        midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (deps >= 30.))[0]
        #midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)))[0]
        
        if len(midx) > 0:
    
            # get binned medians
            bins = arange(2.0, log10(3000), 0.1)
    
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[midx]), log10(amp_plt[midx]))
            
            # normalise at 1000 km
            bidx = where(around(binstrp, 1) == around(2.8,1))[0]
            if len(bidx) == 1:
                norm_amp = amp_plt[midx] / 10**logmedamp[bidx]
                norm_amp_all = concatenate((norm_amp_all, norm_amp))
                norm_rhyp = concatenate((norm_rhyp, rhyp[midx]))
                norm_dep = concatenate((norm_dep, deps[midx]))
                norm_mag = concatenate((norm_mag, mags[midx]))
	
    ################################################################################
    # make pretty
    ################################################################################
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.xlim([100, 3000])
    plt.ylim([1E-3, 1E3])
    plt.grid(which='both')
    #plt.title('Normalised Amplitudes')
    
    plt.ylabel('SA('+str(T)+') in g')
    if ii > 5:
        plt.xlabel('Hypocental Distance (km)')
        
    # use mean c2
    #c2_mean = mean(array(init_c2))
    
    # use period dependent c2
    c2T = smooth_c2[ii]
    
    def fit_q(c, x):
        from numpy import sqrt, log10
        
        ans = c[0] - c[1]*x - c2T*log10(x)
        
        return ans

    
    ################################################################################
    # fit GR + Q 
    ################################################################################
    #ridx = where((norm_rhyp > 1500.) & (norm_mag > 6.5))[0] # test not using M < 6.5 for R > 1500
    #data = odrpack.RealData(delete(norm_rhyp, ridx), log10(delete(norm_amp_all, ridx)))
    
    # get binned data
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyp), log10(norm_amp_all))
    
    # use binned data
    ridx = where(10**medx < 1750.)[0]
    data = odrpack.RealData(10**medx[ridx], logmedamp[ridx])
    
    #ridx = where(norm_rhyp < 2000.)[0]
    #data = odrpack.RealData(norm_rhyp[ridx], log10(norm_amp_all[ridx]))
    
    #data = odrpack.RealData(10**medx, logmedamp)
    afit = odrpack.Model(fit_q)
    odr = odrpack.ODR(data, afit, beta0=[1.0, 1.0])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    c = out.beta
    '''
    if c[1] < 0.:
        c[1] = 0.
    '''
    # now plt
    attenfit = c[0] - c[1]*rrup - c2T*log10(rrup)
    plt.loglog(rrup, 10**attenfit, 'r-', lw=2)
    init_c3.append(c[1]) # Q
    
    # keep intercept

plt.show()
        
################################################################################
# fit period dependent Q
################################################################################

fig = plt.figure(2, figsize=(12, 9))

plt.subplot(311)
plt.semilogx(Tplt, init_c2, 'bo')
plt.semilogx(Tplt, smooth_c2, 'ro')
plt.ylabel('c2 (GR)')

plt.subplot(312)
plt.semilogx(Tplt, init_c3, 'bo')

# smooth Q
smooth_c3 = savitzky_golay(init_c3, 5, 3)
plt.semilogx(Tplt, smooth_c3, 'ro')

plt.ylabel('c3 (Q)')
#plt.semilogx(Tplt, c2_bilin, 'r-', lw=2)

plt.subplot(313)
plt.semilogx(Tplt, init_c0, 'bo')
plt.ylabel('c0 init mag') # i think
plt.xlabel("T (sec)")

plt.show()

#slope, intercept, r_value, p_value, std_err = linregress(Tplt, smooth_c2)

# fit porabola
def inv_porabola(c, x):
    from numpy import sqrt, log10
    
    ans = -1*x**2 - c[0]
    
    return ans
    
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
    
################################################################################
# loop through T and M to get mag scaling
################################################################################
fig = plt.figure(3, figsize=(18,10))
fig = plt.figure(3, figsize=(12,6))

# Tplt = array([0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0]) # secs; PGA = 0.01; PGV = -99

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

for tt, T in enumerate(Tplt):
    init_mw = []
    init_c0 = []

    for i, mplt in enumerate(mrng):
                
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
        
        # get data
        rhyp = dictlist2array(stdict, 'rhyp')
        mags = dictlist2array(stdict, 'mag')
        deps = dictlist2array(stdict, 'dep')
        azim = dictlist2array(stdict, 'azim')
        
        midx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)) & (deps >= 30.))[0]
        
        if len(midx) > 0:
    
            # get binned medians
            bins = arange(2.0, 3.51, 0.1)
    
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[midx]), log10(amp_plt[midx]))
            
            ii += 1
            '''
            ax = plt.subplot(4,4, ii)
            norm = mpl.colors.Normalize(vmin=0, vmax=500)
            sc = plt.scatter(rhyp[midx], amp_plt[midx], c=deps[midx], marker='o', s=50, \
                        cmap='Spectral_r', norm=norm, alpha=0.8)
            
            plt.plot(10**medx, 10**logmedamp, 'rs', ms=6.5)
            ax.set_xscale("log")
            ax.set_yscale("log")
            plt.xlim([100, 3000])
            plt.ylim([1E-7, 1E-1])
            plt.grid(which='both')
            plt.title('MW '+str(mplt))
            
            if ii == 0:
                plt.legend(fontsize=10)
                
            if ii == 1 or ii == 5 or ii == 9 or ii == 13:
                plt.ylabel('SA('+str(Tplt[tidx])+') in g')
            if ii > 10:
                plt.xlabel('Hypocental Distance (km)')
            '''
            # fit m-scaling
            #c3 = init_c3[tidx]
            c3 = smooth_c3[tt]
            c2 = smooth_c2[tt] #_mean
            def fit_atten(c, x):
                from numpy import sqrt, log10
                
                #xref = 100.
                ans = c[0] - c2*log10(x) - c3*x 
                
                return ans
                
            data = odrpack.RealData(rhyp[midx], log10(amp_plt[midx]))
            #data = odrpack.RealData(10**medx, logmedamp)
            
            afit = odrpack.Model(fit_atten)
            odr = odrpack.ODR(data, afit, beta0=[1.0])
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            c = out.beta
            
            '''
            # now plt
            attenfit = c[0] - c3*rrup - c2*log10(rrup)
            plt.loglog(rrup, 10**attenfit, 'k-', lw=1.5)
            '''
            # add to array
            init_c0.append(c[0])
            init_mw.append(mplt)   
            
            #plt.savefig('banda_gmm_fit_'+str(Tplt[tidx])+'.png', fmt='png', bbox_inches='tight')
            #plt.show()
            
    
    # fit straight quadratic
    nn = where(isnan(init_c0)==False)[0]
    m1, m2, m0 = polyfit((array(init_mw)[nn]-6), array(init_c0)[nn], 2)
    m0_array.append(m0)
    m1_array.append(m1)
    m2_array.append(m2)
    
    t_dept_c0.append(array(init_c0))
    t_dept_mw.append(array(init_mw))
    
    '''
    # fit mag scaling for period
    def fit_quadratic_vertex(c, x):
        
        xx = x - 8.
        
        return c[0] * xx**2 + c[1]

    
    data = odrpack.RealData(array(init_mw), array(init_c0))

    magfit = odrpack.Model(fit_quadratic_vertex)
    odr = odrpack.ODR(data, magfit, beta0=[1., 1.])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    m = out.beta
    
    
    mag_per_fit.append(m)
    mag_c0_fit.append(init_c0)
    '''
    
    
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
plt.semilogx(Tplt, array(m0_array), 'bo')
smooth_m0 = savitzky_golay(array(m0_array), 7, 3)
plt.semilogx(Tplt, smooth_m0, 'ro')
plt.title('m0')
# plot M2            
plt.subplot(222)
plt.semilogx(Tplt, array(m1_array), 'bo')
smooth_m1 = savitzky_golay(array(m1_array), 7, 3)
plt.semilogx(Tplt, smooth_m1, 'ro')
plt.title('m1')

# plot M2            
plt.subplot(224)
m_scaling = m0_array + m1_array*(7.0-6)**2 + m2_array*(7.0-6)
plt.semilogx(Tplt, m_scaling, 'bo')
plt.title('m7')



######################################################################################
# refit M2 with smoothed M!
refit_m2 = []
for c0_dat, sm1 in zip(mag_c0_fit, smooth_m0):
	
    
    # fit mag scaling for period
    def fit_m2_vertex(c, x):
        
        xx = x - 8.
        
        return sm1 * xx**2 + c[0]

    data = odrpack.RealData(array(init_mw), array(c0_dat))

    magfit = odrpack.Model(fit_m2_vertex)
    odr = odrpack.ODR(data, magfit, beta0=[1.])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    m = out.beta
    
    refit_m2.append(m[0])

#plt.semilogx(Tplt, array(refit_m2), 'go-')

#smooth_m2 = savitzky_golay(array(mag_per_fit)[:,1], 5, 3)

plt.show()

######################################################################################
# loop throught periods and plot mag scaling
######################################################################################
fig = plt.figure(10, figsize=(18,10))

i = 0
for T, ym, xm in zip(Tplt,t_dept_c0, t_dept_mw):
    plt.subplot(4,5,i+1)
    plt.plot(xm, ym, 'bo')
    plt.ylabel(str(T))
    
    mfit = m0_array[i] + m1_array[i]*(mrng-6)**2 + m2_array[i]*(mrng-6)
    plt.plot(mrng, mfit, 'g-', lw=2.)
    #plt.ylim([4.5, 7])
    
    if i >= 16:
        plt.xlabel('MW')
    
    i += 1

plt.show()


def polyfitmag(c, x):
    from numpy import sqrt, log10
    
    ans = c[0] + c[1]*(x-6)**2
    
    return ans
    
def fit_quadratic_vertex(c, x):
    
    xx = x - 8.
    
    return c[0] * xx**2 + c[1]

'''
fig = plt.figure(4)
plt.plot(init_mw, init_c0, 'bo') # ignore last data point

#data = odrpack.RealData(array(init_mw[0:-1]), array(init_c0[0:-1]))
data = odrpack.RealData(array(init_mw), array(init_c0))

magfit = odrpack.Model(fit_quadratic_vertex)
odr = odrpack.ODR(data, magfit, beta0=[1., 1.])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
m = out.beta
print m

mfit = m[0]*(mrng-8)**2 + m[1]
plt.plot(mrng, mfit, 'g-', lw=2.)


m0, m1, m2 = polyfit(array(init_mw)-6, init_c0, 2)
mfit = m0*(mrng-6)**2 + m1*(mrng-6) + m2
plt.plot(mrng, mfit, 'r-', lw=2.)
print 'DO POLYNOMAL FUNCTION!!!!'
# 

plt.xlabel("MW")
plt.ylabel('c0')
plt.show()
'''

#################################################################################
# export coeffs
#################################################################################

ctxt = '#log Y = c0 + c1*(M-6)**2 + c2*(M-6) - c3*log Rhyp - c4*Rhyp\nT, c0, c1, c2, c3, c4\n'
for i, t in enumerate(Tplt):
    ctxt += ','.join((str(t), str('%0.5f' % m0_array[i]), str('%0.5f' % m1_array[i]), str('%0.5f' % m2_array[i]), \
                      str('%0.5f' % smooth_c2[i]), str('%0.5e' % init_c3[i]))) + '\n'
                      
f = open('ncc_gmm_coeffs.csv', 'wb')
f.write(ctxt)
f.close()

'''

data = odrpack.RealData(log10(Tplt), array(init_c2))
bilin = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data, bilin, beta0=[.5, 3., -0.5])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()

a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = 0. #out.beta[3]

c2_bilin = b + a * log10(Tplt)
hy = b + a * hx
idx = log10(Tplt) > hx
c2_bilin[idx] = c * (log10(Tplt)[idx]-hx) + hy

fig = plt.figure(2)
plt.semilogx(Tplt, init_c2, 'bo')
#c2_fit = slope * Tplt + intercept
#plt.semilogx(Tplt, c2_fit, 'k-', lw=2)

#c2_fit = -1 * log10(Tplt)**2 - c
plt.semilogx(Tplt, c2_bilin, 'r-', lw=2)
plt.show()
'''

"""

        ii += 1
    
        ax = plt.subplot(4,4, ii)
        norm = mpl.colors.Normalize(vmin=0, vmax=500)
        sc = plt.scatter(rhyp[idx], amp_plt[midx], c=deps[midx], marker='o', s=50, \
                    cmap='Spectral_r', norm=norm, alpha=0.8)
        
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.xlim([100, 3000])
        plt.ylim([1E-7, 1E-1])
        plt.grid(which='both')
        plt.title('MW '+str(mplt))
        
        if ii == 0:
            plt.legend(fontsize=10)
            
        if ii == 1 or ii == 5 or ii == 9 or ii == 13:
            plt.ylabel('SA('+str(Tplt)+') in g')
        if ii > 10:
            plt.xlabel('Hypocental Distance (km)')

        ################################################################################
        # now fit binned GR
        ################################################################################
        
        #data = odrpack.RealData(10**medx, logmedamp)#, sx=sx, sy=sy)
        #data = odrpack.RealData(rhyp[idx], log10(amp_plt[idx]))#, sx=sx, sy=sy)
        
        ridx = where(rhyp <= 1000.)[0]
        logmedamp_all, stdbin, medx_all, binstrp, nperbin = get_binned_stats(bins, log10(rhyp[ridx]), log10(amp_plt[ridx]))
        data = odrpack.RealData(10**medx_all, logmedamp_all)
        #data = odrpack.RealData(rhyp, log10(amp_plt))
        
        afit = odrpack.Model(fit_atten)
        odr = odrpack.ODR(data, afit, beta0=[1.0, 1.0, 1.0])
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        c = out.beta
        print c
        
        # now plt
        #attenfit = c[0] - c[1]*(rrup) - c[2]*log10(rrup)
        attenfit = c[0] - c[1]*rrup - c[2]*log10(rrup)
        #attenfit = c[0] - c[1] * (rrup**pwr + c[2]**pwr)**(1./pwr)
        plt.loglog(rrup, 10**attenfit, 'k-', lw=1.5)
        
        ################################################################################
        # now fit binned data for each mag
        ################################################################################
        
        geospread = c[2]
        qfactor = c[1]
        
        def fit_atten_fix_gr_q(c, x):
            from numpy import sqrt, log10
            
            xref = 100.
            ans = c[0] - qfactor*(x) - geospread*log10(x)
            
            return ans
            
        # set mag specific reg
        data = odrpack.RealData(10**medx_all, logmedamp_all)#, sx=sx, sy=sy)
        data = odrpack.RealData(rhyp[idx], log10(amp_plt[idx]))#, sx=sx, sy=sy)
        afit = odrpack.Model(fit_atten_fix_gr_q)
        odr = odrpack.ODR(data, afit, beta0=[1.0])
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        c = out.beta
        print c
        
        # now plt
        attenfit = c[0] - qfactor*(rrup) - geospread*log10(rrup)
        plt.loglog(rrup, 10**attenfit, 'r-', lw=1.5)
        
        if min(rhyp[idx]) < 1000. and mplt < 7.7:
            init_c0.append(c[0])
            #init_c1.append(c[1])
            init_mw.append(mplt)


# make colorbar
plt.gcf().subplots_adjust(right=0.93)
cax = fig.add_axes([0.94,0.33,0.02,0.33]) # setup colorbar axes.
#norm = colors.Normalize(vmin=0, vmax=500)
cb = colorbar.ColorbarBase(cax, cmap='Spectral_r', norm=norm, orientation='vertical', alpha=0.8)
cb.set_label("Hypocentral Depth (km)", fontsize=12)
    
plt.savefig('banda_atten.'+str(Tplt)+'.png', fmt='png', bbox_inches='tight')
plt.show()


################################################################################
# plot mag scaling
################################################################################

fig = plt.figure(2)
plt.plot(init_mw, init_c0, 'bo')
plt.xlabel('MW')
plt.ylabel('c0')

slope, intercept, r_value, p_value, std_err = linregress(init_mw, init_c0)
mscale = slope * mrng +intercept
plt.plot(mrng, mscale, 'k-', lw=2.)
plt.xlim([5.9, 7.7])
plt.show()
"""








