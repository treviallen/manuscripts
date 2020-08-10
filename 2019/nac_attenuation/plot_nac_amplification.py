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
        
    near_field_term = n0 * (log10(rhyp)-hx[0])
    
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

###############################################################################
# start main
###############################################################################

import matplotlib.pyplot as plt
from numpy import exp
from misc_tools import get_mpl2_colourlist
import matplotlib as mpl
mpl.style.use('classic')

fig = plt.figure(1, figsize=(7,7))
ax = plt.subplot(111)
cols = get_mpl2_colourlist()
syms = ['o', '^', 's', 'd']

# set vals
mag = 7.0
rhyp = 600
dep = 200
vs30s = [220, 270, 360, 760] #, 1100]

# calc 760
A19imt_760 = calc_nac_gmm_spectra(mag, rhyp, dep, 760., 'BS')

for i, vs30 in enumerate(vs30s):
    A19imt = calc_nac_gmm_spectra(mag, rhyp, dep, vs30, 'BS')
    
    rat_amp = A19imt['sa'] - A19imt_760['sa']
    
    plt.loglog(A19imt['per'][1:-2], exp(rat_amp)[1:-2], syms[i], ls='-', c=cols[i], lw=2, \
    	         ms=8, mec=cols[i], mfc='none', mew=2, markevery=4, label=str(vs30)+' m/s')
    
plt.legend(loc=1, numpoints=1)
plt.grid(which='both')
plt.xlim([0.05, 10])
plt.ylim([0.9, 10])
yticks = [1, 2, 5, 10]
ax.tick_params(axis='both', labelsize=14)

ax.set_yticks(yticks)
ax.set_yticklabels([str(x) for x in yticks])

plt.xlabel('Spectral Period (s)', fontsize=16)
plt.ylabel('Amplification Factor', fontsize=16)
    

plt.savefig('figures/nac_amplification.png', format='png', dpi=150, bbox_inches='tight')
plt.show()

