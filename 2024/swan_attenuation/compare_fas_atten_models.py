import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, logspace, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, polyfit, isfinite
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl
from get_mag_dist_terms_swan import get_distance_term, get_magnitude_term, get_kappa_term
from misc_tools import get_ga_master_colours_2022, dictlist2array, get_log_xy_locs
from fas_tools import get_dist_atten_allen_2006, get_dist_atten_allen_2007, get_dist_atten_atkinson_boore_2014, get_dist_atten_atkinson_2004
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

###############################################################################
# load coeffs
###############################################################################

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

mcoeffs = pickle.load(open('mag_coeffs.pkl', 'rb' ))

cols = get_ga_master_colours_2022()

###############################################################################
# get kappa dat
###############################################################################

def parse_kappa_data():
    from numpy import array, loadtxt
    
    kapdat = []
    # read parameter file
    lines = open('site_kappa.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        #kap = {'sta':dat[0], 'kappa_0': float(dat[1]), 'kappa_r': float(dat[2])}
        kap = {'sta':dat[0], 'kappa0': float(dat[1]), 'cnt': float(dat[2])}
    
        kapdat.append(kap)
    
    return kapdat
kapdat = parse_kappa_data()

###############################################################################
# get new atten
###############################################################################
rhyps = logspace(0, log10(500), 100)

# set freqs    
fidx = [51, 99] # 1 & 5 Hz

mw = 5

fig = plt.figure(1, figsize=(18,6))
xplt = arange(1,500,1)
for p, fi in enumerate(fidx):
    
    ax = plt.subplot(1,2,p+1)
    
    '''
    get Allen etal (2006) for swwa
    '''
    a06_fds = get_dist_atten_allen_2006(mw, rhyps, freqs[fi]) / 10**0.18 # H-V quick correction 
    plt.loglog(rhyps, a06_fds, '-', c=cols[2], lw=2, label='Allen et al (2006; Burakin)')
    
    '''
    get Allen etal (2006) for sea
    '''
    a07_fds = get_dist_atten_allen_2007(mw, rhyps, freqs[fi]) 
    plt.loglog(rhyps, a07_fds, '-', c=cols[3], lw=2, label='Allen et al (2007; SE Aust)')
    
    '''
    get AB (2014) for ENA
    '''
    #ab14_fds = get_dist_atten_atkinson_boore_2014(mw, rhyps, freqs[fi]) 
    #plt.loglog(rhyps, ab14_fds, '-', c=cols[4], lw=2, label='Atkinson & Boore (2014; ENA)')
    
    '''
    get A04 for ENA
    '''
    a04_fds = get_dist_atten_atkinson_2004(mw, rhyps, freqs[fi]) 
    plt.loglog(rhyps, a04_fds, '-', c=cols[5], lw=2, label='Atkinson (2004; ENA)')
    
    
    '''
    get present study
    '''
    magterm = get_magnitude_term(mw, mcoeffs[fi])
    
    distterms = []
    for r in rhyps:
        distterm = get_distance_term(r, coeffs[fi])
        distterms.append(distterm)
    
    #	get distance independent kappa
    kappa = kapdat[-1]['kappa0'] # default kappa
    
    k_term = log10(exp(-1 * pi * freqs[fi] * kappa))
    
    # correct to source
    atten_dat = (10**(magterm + distterms + k_term)) * 1000 # convert from m-s to mm-s
    
    plt.loglog(rhyps, atten_dat, '-', c=cols[0], lw=2, label='This Study')
    
    plt.xlabel('Hypocentral Distance (km)', fontsize=17)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    if p == 0:
        plt.ylabel('Fourier Displacement Amplitude (mm-s)', fontsize=17)
    plt.xlim([1,500])
    
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    xloc = get_log_xy_locs(xlims, 0.04)
    yloc = get_log_xy_locs(ylims, 0.04)
    props = dict(boxstyle='round', facecolor='w', alpha=1)
    if p == 0:
       fstr = str('%0.2f' % freqs[fi])
    else:
        fstr = str('%0.1f' % freqs[fi])
    plt.text(xloc, yloc, 'f = '+fstr+' Hz', va='bottom', ha ='left', fontsize=16, bbox=props)
    plt.grid(which='both')
       
ax = plt.subplot(1,2,1)
plt.legend(loc=1, fontsize=14)
plt.subplots_adjust(wspace=0.1)
plt.savefig('fas_comparison.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
