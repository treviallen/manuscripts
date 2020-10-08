from numpy import array, arange, sqrt, exp, log, log10, logspace, where, interp, unique, vstack, \
                  argwhere, ones_like, zeros_like, pi, concatenate, mean, nanmean, std, nanstd, isfinite
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from mag_tools import m02mw
from misc_tools import get_log_xy_locs, listdir_extension, savitzky_golay, get_mpl2_colourlist
from os import path
from scipy import interpolate
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
import matplotlib.pyplot as plt
import pickle
import matplotlib as mpl
mpl.style.use('classic')

folder = argv[1] # fds folder
fdsfiles = listdir_extension(folder, 'fds')

interpfreqs = logspace(-1, log10(50), 100)
maxr = 430.
minr = 50.

# get stn from file list
fdsdict = []
for files in fdsfiles:
    stn = files.split('.')[1]
    fdspath = path.join(folder,files)
    lines = open(fdspath).readlines()
    stla = float(lines[3].split('\t')[1])
    stlo = float(lines[4].split('\t')[1])
    sps = float(lines[5].strip().split('\t')[1].strip(' Hz'))
    mag = float(lines[9].split('\t')[1])
    rhyp = float(lines[10].strip().split('\t')[1].strip(' km'))
    
    if rhyp < maxr and rhyp >= minr:
        print(files)
    
        # re-read and remover headers
        freqs = []
        fds = []
        if lines[21].startswith('---'):
            sidx = 22
        else:
            sidx = 21
            
        lines = open(fdspath).readlines()[sidx:]
        for line in lines:
            dat = line.strip().split('\t')
            freqs.append(float(dat[0]))
            fds.append(float(dat[1]))
            
        fdict = {'stn':stn, 'stla':stla, 'stlo':stlo, 'sps':sps, 'mag':mag, \
        	       'rhyp':rhyp, 'freqs':array(freqs), 'fds':array(fds)}
        	       	
        # use spline interpolation smooth and interp to interpfreqs
        n = len(fdict['freqs'])
        #index = argwhere(interpfreqs < abs(fdict['freqs'][n])).flatten()
        #w = ones_like(fdict['freqs'])
        #tck = interpolate.splrep(abs(fdict['freqs']),log(abs(fdict['fds'])),w=w, s=n/5.,k=3)
        #tck = interpolate.splrep(abs(fdict['freqs']),log(abs(fdict['fds'])),k=2, s=3)
        #smoothfft = exp(interpolate.splev(interpfreqs,tck,der=0))
        
        if fdict['stn'] == 'MILA':
           sw = 5
        else:
            sw = 51
        smoothfft = exp(savitzky_golay(log(abs(fdict['fds'])), sw, 3))
        	       	
        fdict['smfreqs'] = interpfreqs
        fdict['smfds'] = exp(interp(log(interpfreqs), log(fdict['freqs']), log(smoothfft)))
        
        fdsdict.append(fdict)
        
        # plotting test
        '''
        plt.loglog(fdict['freqs'],abs(fdict['fds']), 'b-')
        plt.loglog(interpfreqs,smoothfft, 'r-', lw=2)
        plt.title(stn)
        plt.show()
        '''

print('Saving pickle files...')
pklpath = path.join(folder, 'fds.pkl')
pklfile = open(pklpath, 'wb')        
pickle.dump(fdsdict, pklfile, -1)
pklfile.close()


####################################################################################
# now atten correct
####################################################################################

def correct_2012_GR_model(freqs, logfds, r, kappa, region):
    # r = rhyp  
    if region == 'sea':
        # get G(R)
        r1 = 90.
        r2 = 150.
        b1 = -1.33
        b2 = 0.32
        b3 = -1.66
        
        if r <= r1:
            GRfact = b1*log10(r)
        elif r > r1 and r <= r2:
            GRfact = b1*log10(r1) + b2*log10(r/r1)
        elif r > r2:
            GRfact = b1*log10(r1) + b2*log10(r2/r1) + b3*log10(r/r2)
        
        GRfact = abs(GRfact)
        
        # get Q(f)
        c1 = 5.85E-3
        c2 = -0.015
        Qfact = array(c1 + c2*log10(freqs))
        findex = where(freqs <= 2.0)[0]
        Qfact[findex] = 0.
        
        corfds = 10**(GRfact + logfds - log10(exp(Qfact*(r**3 + r1**3)**(1./3.))) - log10(exp(-1 * pi * freqs * kappa))) # in  m-s 
                 
        #corfds = (10**GRfact * 10**logfds) # in  m-s 
            
    return corfds
    
def correct_2007_GR_model(freqs, logfds, r, kappa, region):
    # r = rhyp  
    if region == 'sea':
        # get G(R)
        r1 = 90.
        r2 = 160.
        b1 = -1.3
        b2 = 0.1
        b3 = -1.6
        vs = 3.6
        
        if r <= r1:
            GRfact = b1*log10(r)
        elif r > r1 and r <= r2:
            GRfact = b1*log10(r1) + b2*log10(r/r1)
        elif r > r2:
            GRfact = b1*log10(r1) + b2*log10(r2/r1) + b3*log10(r/r2)
        
        GRfact = abs(GRfact)
        
        # get Q(f)
        c1 = 5.85E-3
        c2 = -0.015
        freqs = array(freqs)
        Qfact = zeros_like(freqs)
        
        fidx = where(freqs <= 3.92)[0]
        Qfact[fidx] = 10**(3.66 - 1.05 * log10(freqs[fidx]))
        fidx = where((freqs > 3.92) & (freqs <= 9.83))[0]
        Qfact[fidx] = 10**(3.01 + 0.03 * log10(freqs[fidx]))
        fidx = where(freqs > 9.83)[0]
        Qfact[fidx] = 10**(2.56 + 0.48 * log10(freqs[fidx]))
        
        Qterm = exp(-1 * pi * freqs * r / (vs * Qfact))
        
        corfds = 10**(GRfact + logfds - log10(Qterm) - log10(exp(-1 * pi * freqs * kappa))) # in m-s 
                 
        #corfds = (10**GRfact * 10**logfds) # in  m-s 
            
    return corfds


####################################################################################
# set def params
####################################################################################

vs = 3.6 # km/s
vsm = vs*1000.
rho = 2800 # kg/m^3
C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.707)

RP=0.55
VHC=0.7071068
FSE=2.0
shear_vel=3.6
density=2.8
R0=1.0
C_USGS = 1/ (RP * VHC * FSE / (4 * pi * density * shear_vel**3 * R0) * 1e-20)

print(C, C_USGS)

nf = 12 # 5 Hz
mindist = minr
maxdist = maxr
minmag = 3.8
minf = 12
maxf = 22

plt.rcParams['pdf.fonttype'] = 42

pklpath = path.join(folder, 'fds.pkl')
pklfile = open(pklpath, 'rb')        
fdsdict = pickle.load(pklfile)

####################################################################################
# get kappas
####################################################################################

kapstn = []
kapval = []
lines = open('kappa_list.dat').readlines()
for line in lines:
    dat = line.strip().split('\t')
    kapstn.append(dat[0].strip())
    kapval.append(float(dat[1].strip()))
        
####################################################################################
# start main
####################################################################################
fig = plt.figure(1, figsize=(7,7))
ax = plt.subplot(111)
stack_logfds = []

# loop through & plot station  data
cmap = plt.get_cmap('hsv_r', len(fdsdict)+1)
#cs = (cmap(arange(len(fdsdict)+1)))
cs = get_mpl2_colourlist()

handles1 = []
labels1 = []
for i, fds in enumerate(fdsdict):
    # correct G(R)
    if fds['rhyp'] < maxdist:
        # correct kappa
        kappa = 0.006 # default kappa
        # get site specific kappa
        for k, stn in enumerate(kapstn):
            if stn == fds['stn']:
                kappa = kapval[k]
        
        print(fds['stn'], fds['rhyp'], kappa, len(fds['fds']))
        
        #if fds['stn'] == 'FSHM':
        corfds = correct_2007_GR_model(fds['smfreqs'], log10(fds['smfds']), fds['rhyp'], kappa, 'sea')
        fdsdict[i]['corfds']  = corfds
        if i <= 9:
            h1, = plt.loglog(fds['smfreqs'],fds['corfds'],'-', c=cs[i], lw=1, label=fds['stn'])
        elif i <= 19:
            h1, = plt.loglog(fds['smfreqs'],fds['corfds'],'--', c=cs[i-10], lw=1, label=fds['stn'])
        else:
            h1, = plt.loglog(fds['smfreqs'],fds['corfds'],'-.', c=cs[i-20], lw=1, label=fds['stn'])
        handles1.append(h1)
        labels1.append(fds['stn'])
        
        # label curves
        minf = 0.4
        maxsf = 0.4
        sfidx = where((interpfreqs >= minf) & (interpfreqs <= maxsf))[0]
        meanamp = exp(nanmean(log(fds['corfds'][sfidx])))
                
        # get data for average mean
        if stack_logfds == []:
            stack_logfds = log(corfds)
        else:
            stack_logfds = vstack((stack_logfds, log(corfds)))

leg1 = plt.legend(handles=handles1, loc=3, fontsize=13, ncol=2)

# get mean of logfds
mean_fds = exp(stack_logfds.mean(axis=0))

# nanmean does not work above
mean_fds = []
for f in range(0, len(stack_logfds[0])):
    idx = where(isfinite(stack_logfds[:, f]))[0]
    mean_fds.append(exp(nanmean(stack_logfds[:, f][idx])))
mean_fds = array(mean_fds)

#mean_fds = exp(stack_logfds)
h2, = plt.loglog(fds['smfreqs'],mean_fds,'--', color='0.2', lw=2, label='Mean Source Spectrum')
plt.grid(which='both', color='0.75')

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
    FittedCurve = log(c[0] / (1 + (exp(f) / (c[1]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve

# call model fit and plot
if folder.startswith('fds/20120720'):
    minf = .5
else:
    minf = .4

maxsf = 0.45
maxf = 15.
fidx = where((interpfreqs >= minf) & (interpfreqs <= maxf))[0]
#sfidx = where((interpfreqs >= minf) & (interpfreqs <= maxsf))[0] # for labelling curves

data = odrpack.RealData(log(interpfreqs[fidx]), log(mean_fds[fidx]))
fitted_brune = odrpack.Model(fit_brune_model)
odr = odrpack.ODR(data, fitted_brune, beta0=[1.,1.])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()
#out.pprint()

omega0 = out.beta[0]
f0 = out.beta[1]
print('f0', f0)

m0 = C * omega0
print('\nMw', m02mw(m0))

# calc stress drop
r0 = 2.34 * vsm / (2 * pi * f0)

sd = 7. * m0 / (16. * r0**3) / 10**6 # in MPa
print('\nSD', sd, 'MPa')
print('\nf0', f0, 'Hz' )

# plot fitted curve
fitted_curve = omega0 / (1 + (interpfreqs / f0)**2)
h3, = plt.loglog(interpfreqs, fitted_curve, 'k-', lw=2., label='Fitted Brune Model')
plt.legend(handles=[h2, h3], loc=1, fontsize=16)
plt.gca().add_artist(leg1)

plt.xlim([0.3, 35])
if folder.startswith('fds/20120720'):
    plt.ylim([1E-6, .03])
else:
    plt.ylim([1E-5, .3])
plt.xlabel('Frequency (Hz)', fontsize=18)
plt.ylabel('Fourier Displacement Spectra (m-s)', fontsize=18)
plt.savefig(folder.split('/')[1]+'.Moe.brunefit.png', format='png', dpi=150, bbox_inches='tight')
plt.show() 

# now fit Brune model with predefined omega0
def fit_brune_model_fixed_m0(c, f):
    from numpy import array, log

    # fit curve
    FittedCurve = log(omega0 / (1 + (exp(f) / (c[0]))**2))
    
    return FittedCurve


# given mag, fit individual recs fro stress drop uncert
f0i = []
sdi = [] 
smw = []
for lfds in stack_logfds:
    
    data = odrpack.RealData(log(interpfreqs[fidx]), lfds[fidx])
    fitted_brune = odrpack.Model(fit_brune_model_fixed_m0)
    odr = odrpack.ODR(data, fitted_brune, beta0=[1.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    
    f0i.append(out.beta[0])
    r0 = 2.34 * vsm / (2 * pi * f0i[-1])
    sdi.append(7. * m0 / (16. * r0**3) / 10**6) # in MPa
    
    # get mag uncert
    fitted_brune = odrpack.Model(fit_brune_model)
    odr = odrpack.ODR(data, fitted_brune, beta0=[1.,1.])
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    #out.pprint()

    omega0 = out.beta[0]
    f0 = out.beta[1]

    m0 = C * omega0
    #print('Mw', m02mw(m0))  
    smw.append(m02mw(m0))

print('\nMW uncert', nanstd(smw))
print('SD uncert', exp(nanstd(log(sdi))), 'MPa')
print('f0 uncert', nanstd(f0i), 'Hz')