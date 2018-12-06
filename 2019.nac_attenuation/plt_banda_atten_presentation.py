# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 08:36:05 2017

@author: u56903
"""
def calc_nac_gmm_spectra(mag, rhyp):
    from numpy import loadtxt, log10, log
    
    coeffile = '//Users//tallen//Documents//Earthquake_Data//2017_NAC_GMM//ncc_gmm_coeffs.csv'
    
    coeffs = loadtxt(coeffile, delimiter=',', skiprows=2)  
    
    T  = coeffs[:,0]
    c0 = coeffs[:,1]
    c1 = coeffs[:,2]
    c2 = coeffs[:,3]
    c3 = coeffs[:,4]
    c4 = coeffs[:,5]
    
    lnsa = log(10**(c0 + c1*(mag-6)**2 + c2*(mag-6) - c3*log10(rhyp) - c4*rhyp))
    
    A17imt = {'per':T, 'sa':lnsa}

    return A17imt

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
    from numpy import array, nan
    from scipy.constants import g

    lines = open(safile).readlines()
    #print safile
    
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
    SA = []
    T = []
    #print lines[14]
    if lines[14].strip() == '':
        azim = float(lines[11].strip().split('\t')[-1].split()[0])
        pga = float(lines[12].strip().split('\t')[-1].split()[0]) / (1000. * g) # convert mm/s**2 to g
        pgv = float(lines[13].strip().split('\t')[-1].split()[0]) / 10. # convert mm/s to cm/s
        
        for line in lines[24:]:
            dat = line.strip().split('\t')
            T.append(float(dat[0]))
            SA.append(float(dat[1]))
    else:
        pga = float(lines[11].strip().split('\t')[-1].split()[0]) / (1000. * g) # convert mm/s**2 to g
        pgv = float(lines[12].strip().split('\t')[-1].split()[0]) / 10. # convert mm/s to cm/s
        azim = nan
        
        for line in lines[23:]:
            dat = line.strip().split('\t')
            print dat
            T.append(float(dat[0]))
            SA.append(float(dat[1]))
        
    T = array(T)
    SA = array(SA) / g

    rec = {'chan':chan, 'datestr':datestr, 'sta':sta, 'stla':stla, 'stlo':stlo, \
           'sr':sr, 'eqla':eqla, 'eqlo':eqlo, 'dep':dep, 'mag':mag, 'rhyp':rhyp, \
           'azim':azim, 'pga':pga, 'pgv':pgv, 'per':T, 'sa':SA, 'safile':safile}
    
    return rec

################################################################################
# import funcs
################################################################################
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors, colorbar
from numpy import where, zeros, exp, log, interp, array, logspace, log10, sqrt, arange
from misc_tools import listdir_extension, dictlist2array
from calc_oq_gmpes import inslab_gsims, scr_gsims
from os import path
mpl.style.use('classic')


################################################################################
# loop through sa files and get data
################################################################################
#usgscsv = 'neic_banda_png_gt_6.csv'
#evdict = parse_usgs_events(usgscsv)

folder = 'psa'
extension = 'psa'
safiles = listdir_extension(folder, extension)
recs = []
for saf in safiles:
    rec = read_sa(path.join(folder, saf))
    
    recs.append(rec)

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
            '''
            if st['ev'].startswith('1995-12-25T04.43') and rec['safile'].startswith('1995'):
                print st['ev'], st['sta'], rec['safile'][0:16] , st['ev'], rec['safile'][0:16] == st['ev']
                print st['ev'], st['sta'], rec['sta'], st['sta'] , rec['sta'] == st['sta'] 
                print st['ev'], st['sta'], rec['chan'], st['chstr'][0], rec['chan'] == st['chstr'][0]
            '''    
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
                
            # get north channel
            if rec['safile'][0:16] == st['ev'] and rec['sta'] == st['sta'] \
                and rec['chan'] == st['chstr'][1]:
                
                ndat = rec['sa']
                
        stdict[i]['geom'] = exp((log(edat) + log(ndat)) / 2.)
    
    # else, get one or other horizontal
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

################################################################################
# now get data in given mag range
################################################################################

mrng = arange(6.0, 8.95, 0.1)
mrng = arange(6.0, 8., 0.25)
mpltrng = 0.25 / 2.
Tplt = 0.5 # secs

# compute inslab gmpes
rrup = logspace(2, log10(3000))
dep = 350.
fig = plt.figure(1, figsize=(10,13))
ii = 0    
for i, mplt in enumerate(mrng):
            
    # for each record, log interpolate across period
    aplt = []
    for j, st in enumerate(stdict):
        aplt.append(exp(interp(log(Tplt), log(st['per']), log(st['geom']))))
    aplt = array(aplt)
    
    # get data
    rhyp = dictlist2array(stdict, 'rhyp')
    mags = dictlist2array(stdict, 'mag')
    deps = dictlist2array(stdict, 'dep')
    azim = dictlist2array(stdict, 'azim')
    
    idx = where((mags >= (mplt-mpltrng)) & (mags < (mplt+mpltrng)))[0]
    
    
    if len(idx) > 0:
        ii += 1
    
        ax = plt.subplot(3,2, ii)
        norm = mpl.colors.Normalize(vmin=0, vmax=500)
        sc = plt.scatter(rhyp[idx], aplt[idx], c=deps[idx], marker='o', s=50, \
                    cmap='Spectral_r', norm=norm, alpha=0.8)
        
        ab03 = []
        gea05 = []
        mp10 = []
        a15 = []
        ab06 = []
        bea14 = []
        sea09 = []
        a12 = []
        a17 = []
        for rr in rrup:
            rjb = sqrt(rr**2 - 50**2)
            YYea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, Aea15imt \
                = inslab_gsims(mplt, dep, 40, 30., 90., rr, rjb, 760)
                
            Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt, YA15imt, SP16imt \
                = scr_gsims(mplt, dep, 40, 30., 90., rr, rjb, 760)
            
            # subduction gmms
            ab03.append(interp(Tplt, AB03imt['per'], AB03imt['sa']))
            gea05.append(interp(Tplt, Gea05imt['per'], Gea05imt['sa']))
            mp10.append(interp(Tplt, MP10imt['per'], MP10imt['sa']))
            a15.append(interp(Tplt, Aea15imt['per'], Aea15imt['sa']))
            
            # scr gmms
            ab06.append(interp(Tplt, AB06imt['per'], AB06imt['sa']))
            sea09.append(interp(Tplt, Sea09imt['per'], Sea09imt['sa']))
            a12.append(interp(Tplt, A12imt['per'], A12imt['sa']))
            
            # plt nga-w2
            bea14.append(interp(Tplt, Bea14imt['per'], Bea14imt['sa']))
            
            A17imt = calc_nac_gmm_spectra(mplt, rr)
            a17.append(interp(Tplt, A17imt['per'], A17imt['sa']))

        # plt AB03
        plt.loglog(rrup, exp(ab03), '-', c='darkred', lw=1.5, label='AB03')
        
        # plt Gea05
        plt.loglog(rrup, exp(gea05), '-', c='orangered', lw=1.5, label='Gea05')
        
        # plt MP10
        plt.loglog(rrup, exp(mp10), '-', c='peru', lw=1.5, label='MP10')
        
        # plt Aea15
        plt.loglog(rrup, exp(a15), '-', c='gold', lw=1.5, label='Aea15')
        
        # plt AB06
        plt.loglog(rrup, exp(ab06), '-', c='seagreen', lw=1.5, label='AB06')
        
        # plt Sea09
        plt.loglog(rrup, exp(sea09), '-', c='cyan', lw=1.5, label='Sea09')
        
        # plt A12
        plt.loglog(rrup, exp(a12), '-', c='dodgerblue', lw=1.5, label='A12')
        
        # plt Bea14
        plt.loglog(rrup, exp(bea14), '-', c='magenta', lw=1.5, label='Bea14')
        
        # plt A17
        plt.loglog(rrup, exp(a17), '-', c='k', lw=1.5, label='A18')
        
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.xlim([100, 3000])
        plt.ylim([1E-7, 1E-1])
        plt.grid(which='both')
        plt.title('MW '+str(mplt))
        
        if ii == 1:
            plt.legend(fontsize=10, loc=3)
            
        if ii == 1 or ii == 5 or ii == 9 or ii == 13:
            plt.ylabel('SA('+str(Tplt)+') in g')
        if ii > 10:
            plt.xlabel('Hypocental Distance (km)')

# make colorbar
plt.gcf().subplots_adjust(right=0.93)
cax = fig.add_axes([0.94,0.33,0.02,0.33]) # setup colorbar axes.
#norm = colors.Normalize(vmin=0, vmax=500)
cb = colorbar.ColorbarBase(cax, cmap='Spectral_r', norm=norm, orientation='vertical', alpha=0.8)
cb.set_label("Hypocentral Depth (km)", fontsize=12)
    
plt.savefig('pres_banda_atten.'+str(Tplt)+'.png', fmt='png', bbox_inches='tight')
plt.show()





















