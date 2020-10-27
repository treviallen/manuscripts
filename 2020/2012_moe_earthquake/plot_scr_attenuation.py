# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder, t):
    from fnmatch import filter
    from os import path, walk, system
    from numpy import amax
    from scipy.constants import g

    zfile = ''
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn) >= 0 and filename.find('HE') >= 0:
                efile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('NE') >= 0:
                efile = path.join(root, filename)

            if filename.find(stn) >= 0 and filename.find('HN') >= 0:
                nfile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('NN') >= 0:
                nfile = path.join(root, filename)
            
            if filename.find(stn) >= 0 and filename.find('HZ') >= 0:
                zfile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('NZ') >= 0:
                zfile = path.join(root, filename)
                
            
    # get geomean
    try:
        T, SAe, PGAe = read_psa(efile)
        
        try:
            T, SAn, PGAn = read_psa(nfile)
            
            # get geometric mean and convert to g
            print(efile, nfile)
            geomean = exp((log(SAe) + log(SAn)) / 2.) / g # convert from m/s**2 to g
            PGA = max([PGAe, PGAn])
            
        
        # just use E comp
        except:
            print(efile)
            T, geomean, PGA = read_psa(efile)
            geomean /= g

    except:
        print(zfile)
        T, geomean, PGA = read_psa(zfile)
        geomean /= g
        
    return T, geomean, PGA

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array
    from scipy.constants import g

    lines = open(psafile).readlines()
    
    pga = (float(lines[11].split('\t')[-1].split()[0]) / 1000) / g
    #print('pga', pga
    
    # remove header lines
    if lines[22].startswith('---'):
        stidx = 23
    else:
        stidx = 24

    SA = []
    T = []
    for line in lines[stidx:]:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA), pga
    
# get Rhyp
def get_rhyp_from_psa(psafile):
    line = open(psafile).readlines()[10] # read only line 11
    rstr = line.split('\t')[1]
    return float(rstr.split()[0])

# get event id from psa files
def get_evid(stns):
    for stn in stns:
        for root, dirnames, filenames in walk(folder):
            for filename in filenames:
                if filename.find(stn) >= 0:
                    evid = filename.split('.')[0]
    return evid
'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp, unique, vstack, nan
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from calc_oq_gmpes import inslab_gsims, scr_gsims, tang2019_cam_gsim, \
                          nga_east_mean, get_station_vs30, adjust_gmm_with_SS14, \
                          adjust_gmm_with_nga_east
from misc_tools import get_log_xy_locs, remove_last_cmap_colour
from mapping_tools import distance
import matplotlib as mpl
mpl.style.use('classic')

'''
usage: run plot_scr_attenuation.py Moe_5.4/psa True
'''

plt.rcParams['pdf.fonttype'] = 42
import matplotlib 
matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 

folder = argv[1]
prefix = argv[2]
colTrue = 'True'

# set event details
if prefix.startswith('201206'):
    mag  = 5.14
    dep = 18.0
    eqlat = -38.259
    eqlon = 146.290

elif prefix.startswith('201207'):
    mag  = 4.32
    dep = 12.41
    eqlat = -38.231
    eqlon = 146.215

ztor = dep-1 # guess
rake = 90. # USGS CMT
dip  = 30.

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

# get stn from file list
stns = []
stla = []
stlo = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            lines = open(path.join(root,filename)).readlines()
            stla.append(float(lines[3].split('\t')[1]))
            stlo.append(float(lines[4].split('\t')[1]))

# get unique stns and locs
tstn = [stns[0]]
tla = [stla[0]]
tlo = [stlo[0]]
for i, s1 in enumerate(stns):
    addstn = True
    for s2 in tstn:
        if s1 == s2:
            addstn = False
    
    if addstn == True:
        tstn.append(s1)
        tla.append(stla[i])
        tlo.append(stlo[i])
        
# recast stn locs
stns = tstn
stla = tla
stlo = tlo
        
#stns = unique(array(stns))

# now get distances from psa files
dist = []
for stn in stns:
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn) >= 0:
                d = get_rhyp_from_psa(path.join(folder, filename))
    dist.append(d)
    
# get evid
evid = get_evid(stns)
    
# set site details
vs30 = 760.

rjb = logspace(log10(8.),log10(600.), 50)
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb

fig = plt.figure(1, figsize=(11, 13))

ncols = 8
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

titles = ['Sa(0.01)','Sa(0.1)', 'Sa(0.2)','Sa(0.5)','Sa(1.0)','Sa(2.0)']
Tplot = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0]
props = dict(boxstyle='round', facecolor='w', alpha=1)

'''
titles = ['Sa(0.01)','Sa(0.2)','Sa(1.0)','Sa(2.0)']
Tplot = [0.01, 0.2, 1.0, 2.0]

titles = ['PGA','Sa(0.2)','Sa(1.0)','Sa(2.0)']
Tplot = [0.0, 0.2, 1.0, 2.0]
'''
#Tplot = [0.0]
# loop thru periods
for j, t in enumerate(Tplot):
    ax = plt.subplot(3, 2, j+1)
    Tea02r = []
    C03r = []
    AB06r = []
    Sea09r = []
    Pea11r = []
    A12r = []
    AA13r = []
    Bea14r = [] 
    YA15r = []
    G90r = []
    Tea19r = []
    NGAEr = []
    for i,r in enumerate(rrup):

        # get ground motion estimates from GMPEs
        Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt \
                 = scr_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)
        
        Tea19imt = tang2019_cam_gsim(mag, dep, rrup[i], vs30)
        
        # adjust some GMMs using Seyhan & Stewart 2014
        A12imt = adjust_gmm_with_SS14(A12imt, 820., vs30)
        Sea09imt = adjust_gmm_with_SS14(Sea09imt, 865., vs30)
        
        # get mean USGS NGA-E
        nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup[i], vs30)
        
        # adjust NGA-E from 3000 -> target
        nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
        
        syms = ['s', 'x', '+', '1', '2', 'o']
        
        if t == 0.0:
            #Tea02r.append(Tea02imt['pga'][0])
            Tea02r.append(Tea02imt['sa'][0])
            #C03r.append(C03imt['pga'][0])
            C03r.append(C03imt['sa'][0])
            AB06r.append(AB06imt['pga'][0])
            Sea09r.append(Sea09imt['pga'][0])
            #Pea11r.append(Pea11imt['pga'][0])
            #Pea11r.append(nan)
            A12r.append(A12imt['sa'][0])
            #AA13r.append(AA13imt['pga'][0])
            Bea14r.append(Bea14imt['pga'][0])
            NGAEr.append(nga_e_imt['pga'][0])
            Tea19r.append(Tea19imt['pga'][0])
            #G90r.append(log((12.2 * exp(1.04*5.35) * r**-1.18)/750))
        else:
            #ti = get_T_index(Tea02imt, t)
            # interpolate log values to correct period
            Tea02r.append(interp(t, Tea02imt['per'], Tea02imt['sa']))
            C03r.append(interp(t, C03imt['per'], C03imt['sa']))
            AB06r.append(interp(t, AB06imt['per'], AB06imt['sa']))
            Sea09r.append(interp(t, Sea09imt['per'], Sea09imt['sa']))
            #Pea11r.append(interp(t, Pea11imt['per'], Pea11imt['sa']))
            A12r.append(interp(t, A12imt['per'], A12imt['sa']))
            #AA13r.append(interp(t, AA13imt['per'], AA13imt['sa']))
            Bea14r.append(interp(t, Bea14imt['per'], Bea14imt['sa']))
            NGAEr.append(interp(t, nga_e_imt['per'], nga_e_imt['sa']))
            Tea19r.append(interp(t, Tea19imt['per'], Tea19imt['sa']))
            #G90r.append(log((12.2 * exp(1.04*5.35) * r**-1.18)/750))

    if colTrue == 'True':
        '''
        h1 = plt.loglog(rjb, exp(Tea02r), '.-' , lw=1., color=cs[0])
        h2 = plt.loglog(rjb, exp(C03r),   '--', lw=1., color=cs[1])
        h3 = plt.loglog(rjb, exp(AB06r),  '-.', lw=1., color=cs[2])
        h4 = plt.loglog(rjb, exp(Sea09r), '-' , lw=1., color=cs[3])
        h5 = plt.loglog(rjb, exp(Pea11r), '-' , lw=2.5, color=cs[4])
        h6 = plt.loglog(rjb, exp(A12r),   '-.', lw=2.5, color=cs[5])
        h7 = plt.loglog(rjb, exp(Bea14r), '--', lw=2.5, color=cs[6])
        h8 = plt.loglog(rjb, exp(YA15r),  ':' , lw=2.5, color=cs[7])
        '''

        #h1 = plt.loglog(rjb, exp(Tea02r), '-', lw=1.5, color=cs[0])
        #h2 = plt.loglog(rjb, exp(C03r),   '-', lw=1.5, color=cs[1])
        h1 = plt.loglog(rjb, exp(AB06r), syms[0], ls='-', lw=1., color=cs[0], \
                        ms=5, mec=cs[0], mfc='none',  mew=1., markevery=8)
        h2 = plt.loglog(rjb, exp(Sea09r), syms[1], ls='-', lw=1., color=cs[1], \
                        ms=5, mec=cs[1], mfc='none',  mew=1., markevery=8)
        #h5 = plt.loglog(rjb, exp(Pea11r), '-', lw=1.5, color=cs[4])
        h3 = plt.loglog(rjb, exp(A12r), syms[2], ls='-', lw=1., color=cs[2], \
                        ms=5, mec=cs[2], mfc='none',  mew=1., markevery=8)
        h4 = plt.loglog(rjb, exp(Bea14r), syms[3], ls='-', lw=1., color=cs[3], \
                        ms=5, mec=cs[3], mfc='none',  mew=1., markevery=8)
        h5 = plt.loglog(rjb, exp(NGAEr), syms[4], ls='-', lw=1., color=cs[5], \
                        ms=5, mec=cs[5], mfc='none',  mew=1., markevery=8)
        h6 = plt.loglog(rjb, exp(Tea19r), syms[5], ls='-', lw=1., color=cs[6], \
                        ms=5, mec=cs[6], mfc='none',  mew=1., markevery=8)
        
        #h11 = plt.loglog(rjb, exp(AA13r),  '-', lw=2.5, color='k')
        #if t <= 0.01:
        #    h11 = plt.loglog(rjb, exp(G90r),  '-', lw=2.5, color='k')
        
        
    # get recorded SeisComP3 data and plot
    for k, s in enumerate(stns):
        #try:
        # get station details
        vs30, isproxy, usgsvs, asscmvs, kvs, stla, stlo = get_station_vs30(s)
        
        # now calculate distance on the fly
        starepi = distance(eqlat, eqlon, stla, stlo)[0]
        starhyp = sqrt(starepi**2 + dep**2)
        print(starepi)
        
        T, geomean, pga = get_site_geomean(s, folder, t)
        # get interpolated value
        if t == 0.0:
            specval = pga
        else:
            specval = interp(t, T, geomean)
            
        if stla > -40.:
            h7 = plt.loglog(starepi, specval, 'k+', markersize=10, markeredgewidth=1.75)
        else: # for tassie sites
            h8 = plt.loglog(starepi, specval, 'x', color='0.5', markersize=8, markeredgewidth=1.75)

        #except:
        #    print('Cannot find data for site:', s)

    if j >= 4:
        plt.xlabel('$\mathregular{R_{JB}}$ (km)', fontsize=16)
    
    if j == 0 or j == 2 or j == 4:
        plt.ylabel('Spectral Acceleration (g)', fontsize=16)
        
    plt.xlim([8, 600])
    if evid.startswith('20120720'):
        if t < 0.5:
            plt.ylim([2E-5, 0.2])
        elif t >= 0.5 and t < 1.0:
            plt.ylim([1E-5, 0.1])
        elif t >= 1. and t < 2.0:
            plt.ylim([1E-6, 0.01])
        else:
            plt.ylim([2E-7, 0.002]) 
    else:
        if t < 0.5:               
            plt.ylim([1E-4, 1])   
        elif t >= 0.5 and t <= 1.0:             
            plt.ylim([1E-5, 0.1])
        else:                     
            plt.ylim([1E-6, 0.01])

    #plt.title(titles[j])
    xtxt = get_log_xy_locs(ax.get_xlim(), 0.95)
    ytxt = get_log_xy_locs(ax.get_ylim(), 0.95)
    plt.text(xtxt, ytxt, titles[j], size=17, horizontalalignment='right', verticalalignment='top', weight='normal', bbox=props)
    plt.grid(which='both', color='0.5')
    ylims = ax.get_ylim()
    plt.text(7., ylims[1]*1.25, letters[j], va='bottom', ha ='right', fontsize=16)

    if j == 1:
        #plt.legend((h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0], h8[0], h11[0], h9[0], h10[0]), \
        #           ['T02', 'C03','AB06','Sea09','Pea11','A12','Bea14','YA15', 'Gea90', 'Data', 'Data (TAS)'],loc=3,numpoints=1,fontsize=9.)
        plt.legend((h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0], h8[0]), \
                   ['AB06','Sea09', 'A12','Bea14', 'NGA-E', 'Tea20', 'Data', 'Data (TAS)'],loc=3,numpoints=1,fontsize=9.)
        #plt.legend(['T02', 'C03','AB06','Sea09','Pea11','A12','Data'],loc=3,numpoints=1,fontsize=10)

plt.savefig(evid+'_atten.png', format='png', dpi=300, bbox_inches='tight')
plt.show()