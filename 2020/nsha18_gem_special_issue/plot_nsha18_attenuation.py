
'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp, unique, vstack, nan
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from calc_oq_gmpes import scr_gsims, crustal_gsims, gaull1990_gsim
from misc_tools import get_log_xy_locs
import matplotlib as mpl
mpl.style.use('classic')

plt.rcParams['pdf.fonttype'] = 42
import matplotlib 
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 

mag  = 6.0
dep = 7
ztor = 7. # guess
rake = 90. # USGS CMT
dip  = 30.

# set site details
vs30 = 760.
rjb = logspace(0,2.7,100)
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb

fig = plt.figure(1, figsize=(9, 8))
ncols = 9
cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:3], cs[4:]))

titles = ['PGA','Sa(0.5)','Sa(1.0)','Sa(2.0)']
Tplot = [0.0]

#Tplot = [0.0]
# loop thru periods
for j, t in enumerate(Tplot):
    ax = plt.subplot(1, 1, j+1)
    Tea02r = []
    C03r = []
    AB06r = []
    Sea09r = []
    Sea09YCr = []
    Pea11r = []
    A12r = []
    AA13r = []
    Bea14r = [] 
    YA15r = []
    G90WAr = []
    G90SEAr = []
    G90WArPGV = []
    G90SEArPGV = []
    
    for i,r in enumerate(rrup):

        # get ground motion estimates from GMPEs
        Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt , YA15imt, SP16imt \
            = scr_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)
        # get one of the NGA-W2 GMMs
        # get ground motion estimates from GMPEs
        G90WAimt, G90SEAimt, G90INDimt, G90WA_PGVimt, G90SEA_PGVimt, G90IND_PGVimt = gaull1990_gsim(mag, dep, rrup[i])
        if t == 0.0:
            #Tea02r.append(Tea02imt['pga'][0])
            Tea02r.append(Tea02imt['sa'][0])
            #C03r.append(C03imt['pga'][0])
            C03r.append(C03imt['sa'][0])
            AB06r.append(AB06imt['pga'][0])
            Sea09r.append(Sea09imt['pga'][0])
            Sea09YCr.append(Sea09YCimt['pga'][0])
            Pea11r.append(nan)
            A12r.append(A12imt['sa'][0])
            Bea14r.append(Bea14imt['pga'][0])
            YA15r.append(YA15imt['pga'][0])
            G90WArPGV.append(G90WA_PGVimt['pga'][0])
            G90SEArPGV.append(G90SEA_PGVimt['pga'][0])
            G90WAr.append(G90WAimt['pga'][0])
            G90SEAr.append(G90SEAimt['pga'][0])
        else:
            #ti = get_T_index(Tea02imt, t)
            # interpolate log values to correct period
            Tea02r.append(interp(t, Tea02imt['per'], Tea02imt['sa']))
            C03r.append(interp(t, C03imt['per'], C03imt['sa']))
            AB06r.append(interp(t, AB06imt['per'], AB06imt['sa']))
            Sea09r.append(interp(t, Sea09imt['per'], Sea09imt['sa']))
            Pea11r.append(interp(t, Pea11imt['per'], Pea11imt['sa']))
            A12r.append(interp(t, A12imt['per'], A12imt['sa']))
            AA13r.append(interp(t, AA13imt['per'], AA13imt['sa']))
            Bea14r.append(interp(t, Bea14imt['per'], Bea14imt['sa']))
            YA15r.append(interp(t, YA15imt['per'], YA15imt['sa']))
            G90r.append(log((12.2 * exp(1.04*5.35) * r**-1.18)/750))

    plt.loglog(rjb, exp(G90WAr), '-', lw=2.0, color='k', label='Gaull et al. (1990; WA)')
    plt.loglog(rjb, exp(G90WArPGV), '--', lw=2.0, color='k', label='Gaull et al. (1990; WA[PGV])')
    plt.loglog(rjb, exp(G90SEAr),  '-', lw=2.0, color=cs[0], label='Gaull et al. (1990; SEA)')
    plt.loglog(rjb, exp(G90SEArPGV), '--', lw=2.0, color=cs[0], label='Gaull et al. (1990; SEA[PGV])')
    plt.loglog(rjb, exp(AB06r), '-', lw=2.0, color=cs[2], label='Atkinson & Boore (2006)')
    plt.loglog(rjb, exp(Sea09r),   '-', lw=2.0, color=cs[3], label='Somerville et al. (2009; Non-craton)')
    plt.loglog(rjb, exp(Sea09YCr), '-', lw=2.0, color='gold', label='Somerville et al. (2009; Yilgarn Craton)')
    plt.loglog(rjb, exp(A12r),  '-', lw=2.0, color=cs[5], label='Allen (2012)')
    plt.loglog(rjb, exp(Bea14r),  '-', lw=2.0, color=cs[6], label='Boore et al. (2014)')

    plt.xlabel(r'$\mathregular{R_{JB}}$ (km)', fontsize=18)
    plt.ylabel('PGA (g)', fontsize=18)
    plt.xlim([1, 500])
    if t < 0.5:               
        plt.ylim([1E-3, 2])   
    elif t > 1.0:             
        plt.ylim([1E-5, 0.01])
    else:                     
        plt.ylim([1E-4, 0.1]) 

    #plt.title(titles[j])
    xtxt = get_log_xy_locs(ax.get_xlim(), 0.95)
    ytxt = get_log_xy_locs(ax.get_ylim(), 0.95)
    #plt.text(xtxt, ytxt, titles[j], size=18, horizontalalignment='right', verticalalignment='top', weight='normal')
    plt.grid(which='both', color='0.5')

    plt.legend(loc=3,numpoints=1,fontsize=13.25)
        
plt.savefig('nsha18_atten.png', format='png', dpi=800, bbox_inches='tight')
plt.show()