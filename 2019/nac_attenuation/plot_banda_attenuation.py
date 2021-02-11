def calc_nac_gmm_spectra(mag, rhyp, dep, vs30, region):
    from numpy import loadtxt, log10, log
    
    if region == 'BS':
        coeffile = 'ncc_gmm_coeffs.BS.csv'
    elif region == 'NGH':
        coeffile = 'ncc_gmm_coeffs.NGH.csv'
    elif region == 'OBE':
        coeffile = 'ncc_gmm_coeffs.OBE.csv'
    
    coeffs = loadtxt(coeffile, delimiter=',', skiprows=2)  
    
    T  = coeffs[:,0] #[2:-2]
    c0 = coeffs[:,1] #[2:-2]
    c1 = coeffs[:,2] #[2:-2]
    c2 = coeffs[:,3] #[2:-2]
    c3 = coeffs[:,4] #[2:-2]
    c4 = coeffs[:,5] #[2:-2]
    d0 = coeffs[:,6] #[2:-2]
    d1 = coeffs[:,7] #[2:-2]
    d2 = coeffs[:,8] #[2:-2]
    d3 = coeffs[:,9] #[2:-2]
    n0 = coeffs[:,10]
    hx = coeffs[:,11] #[2:-2]
    
    logdep = log10(dep)
#    lnsa = c0 + c1*(mag-6)**2 + c2*(mag-6) - c3*log10(rhyp) - c4*rhyp # \
#               + (d0 + d1*logdep**3 + d2*logdep**2 + d3*logdep)
    '''
    Rhyp <= hx: ln Y = c0 + c1*(M-6)**2 + c2*(M-6) + \
    (c3*hx +  c4*(log10(Rhyp)-hx)) + (d0 + d1*log10(h)**3 + d2*log10(h)**2 + d3*log10(h))
    '''
    mag_term = c0 + c1*(mag-6)**2 + c2*(mag-6)
    
    if log10(rhyp) < hx[0]:
        atten_term = c3 * log10(rhyp)
    else:
        hy = c3 * hx[0]
        atten_term = c4 * (log10(rhyp)-hx[0]) + hy
    
    # get near source correction
    if log10(rhyp) < hx[0]:
        near_field_term = n0 * (log10(rhyp)-hx[0])
    else:
        near_field_term = 0.0
    
    dep_term = d0 + d1*logdep**3 + d2*logdep**2 + d3*logdep
    
    # get site coefs
    sitefile = 'nac_site_amp_coeffs.csv'
    coeffs = loadtxt(sitefile, delimiter=',', skiprows=1)  
    
    T  = coeffs[:,0]
    s0 = coeffs[:,1]
    s1 = coeffs[:,2]
    	
    site_term = s0 + s1 / (log10(vs30) - log10(150))
    
    lnsa = mag_term + atten_term + dep_term + near_field_term + site_term
           
    A19imt = {'per':T, 'sa':lnsa}

    return A19imt

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp, unique, vstack, nan
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from calc_oq_gmpes import scr_gsims, inslab_gsims, nga_east_mean, adjust_gmm_with_SS14, \
                          adjust_gmm_with_nga_east
from misc_tools import get_log_xy_locs, get_mpl2_colourlist
import matplotlib as mpl
mpl.style.use('classic')

plt.rcParams['pdf.fonttype'] = 42
import matplotlib 
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 

mag = 7.0
dep = 150
ztor = 180. # guess
rake = -90. # USGS CMT
dip  = 30.

# set site details
vs30 = 760.
rjb = logspace(2,log10(1500),75)
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
#rjb = rrup

fig = plt.figure(1, figsize=(18, 8))
'''
ncols = 9
cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'

cmap, zvals = cpt2colormap(cptfile, ncols)

cs = (cmap(arange(ncols)))
'''
cs = get_mpl2_colourlist()
#cs = vstack((cs[0:3], cs[4:]))

titles = ['Sa(0.2)','Sa(2.0)','Sa(2.0)']
Tplot = [0.2, 2.0]
letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
syms = ['o', '^', 's', 'd', 'v', '<', 'h']
props = dict(boxstyle='round', facecolor='w', alpha=1)

#Tplot = [0.0]
# loop thru periods
for j, t in enumerate(Tplot):
    ax = plt.subplot(1, 2, j+1)
    Tea02r = []
    C03r = []
    AB06r = []
    AB03r = []
    Sea09r = []
    Sea09YCr = []
    Pea11r = []
    A12r = []
    AA13r = []
    Bea14r = [] 
    MP10r = []
    NGAEr = []
    Aea16sr = []
    A20_BSr = []
    
    for i, r in enumerate(rrup):

        # get ground motion estimates from GMPEs
        Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt , SP16imt \
            = scr_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)
        
        A12imt = adjust_gmm_with_SS14(A12imt, 820., vs30)
        
        nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup[i], vs30)
        
        # adjust NGA-E from 3000 -> target
        nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
        
        # get in-slab models
        Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06imt, MP10imt, Aea16imt, Zea16imt \
            = inslab_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)
        
        # get BS model
        A20imt_BS = calc_nac_gmm_spectra(mag, rrup[i], dep, vs30, 'BS') # use rrup
        
        if t == 0.0:
            #Tea02r.append(Tea02imt['pga'][0])
            Tea02r.append(Tea02imt['sa'][0])
            #C03r.append(C03imt['pga'][0])
            C03r.append(C03imt['sa'][0])
            AB06r.append(AB06imt['pga'][0])
        else:
            #ti = get_T_index(Tea02imt, t)
            # interpolate log values to correct period
            A12r.append(interp(t, Tea02imt['per'], Tea02imt['sa']))
            if i < len(rjb)-1:
                NGAEr.append(interp(t, nga_e_imt['per'], nga_e_imt['sa']))
            A20_BSr.append(interp(t, A20imt_BS['per'], A20imt_BS['sa']))
            Aea16sr.append(interp(t, Aea16imt['per'], Aea16imt['sa']))
            MP10r.append(interp(t, MP10imt['per'], MP10imt['sa']))
            Bea14r.append(interp(t, Bea14imt['per'], Bea14imt['sa']))
            AB06r.append(interp(t, AB06imt['per'], AB06imt['sa']))
            AB03r.append(interp(t, AB03imt['per'], AB03imt['sa']))
            

    plt.loglog(rjb, exp(Bea14r), syms[0], ls='-', c=cs[0], lw=2, \
               ms=8, mec=cs[0], mfc='none', mew=2, markevery=8, label='Boore et al. (2014)')
    plt.loglog(rjb, exp(AB06r), syms[1], ls='-', c=cs[1], lw=2, \
               ms=8, mec=cs[1], mfc='none', mew=2, markevery=8, label='Atkinson & Boore (2006)')
    plt.loglog(rjb, exp(A12r), syms[2], ls='-', c=cs[2], lw=2, \
               ms=8, mec=cs[2], mfc='none', mew=2, markevery=8, label='Allen (2012)')
    plt.loglog(rjb[:-1], exp(NGAEr), syms[3], ls='-', c=cs[3], lw=2, \
               ms=8, mec=cs[3], mfc='none', mew=2, markevery=8, label='Goulet et al. (2017)')
    '''
    plt.loglog(rjb, exp(MP10r), syms[2], ls='-', c=cs[2], lw=2, \
    	         ms=8, mec=cs[2], mfc='none', mew=2, markevery=8, label='Megawati & Pan (2010; interface)')
    '''
    plt.loglog(rjb, exp(AB03r), syms[4], ls='-', c=cs[4], lw=2, \
               ms=8, mec=cs[4], mfc='none', mew=2, markevery=8, label='Atkinson & Boore (2003; slab)')
    plt.loglog(rjb, exp(Aea16sr), syms[5], ls='-', c=cs[5], lw=2, \
               ms=8, mec=cs[5], mfc='none', mew=2, markevery=8, label='Abrahamson et al. (2016; slab)')
    plt.loglog(rjb, exp(A20_BSr), syms[6], ls='-', c=cs[6], lw=2, \
               ms=8, mec=cs[6], mfc='none', mew=2, markevery=8, label='Banda Sea Model')
    
    plt.xlabel(r'$\mathregular{R_{JB}}$ (km)', fontsize=18)
    plt.xlim([200, 1500])
    if t <= 0.5:               
        plt.ylim([1E-5, 0.1])   
    elif t > 1.0:             
        plt.ylim([1E-5, 0.1])
    else:                     
        plt.ylim([1E-4, 0.1]) 

    xtxt = get_log_xy_locs(ax.get_xlim(), 0.95)
    ytxt = get_log_xy_locs(ax.get_ylim(), 0.95)
    plt.text(xtxt, ytxt, titles[j], size=17, horizontalalignment='right', verticalalignment='top', weight='normal', bbox=props)
    plt.grid(which='both', color='0.5')
    ylims = ax.get_ylim()
    plt.text(175., ylims[1]*1.25, letters[j], va='bottom', ha ='right', fontsize=20)
    
    # set xticks
    xticks = [200, 300, 500, 1000, 1500]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks])
    
    if j == 0:
        plt.ylabel('Spectral Acceleration (g)', fontsize=18)
        plt.legend(loc=3,numpoints=1,fontsize=13.5)
        
        
plt.savefig('figures/banda_atten.png', format='png', dpi=300, bbox_inches='tight')
plt.show()