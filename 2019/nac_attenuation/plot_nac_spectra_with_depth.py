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
    '''
    sitefile = 'nac_site_amp_coeffs.csv'
    coeffs = loadtxt(sitefile, delimiter=',', skiprows=1)  
    
    T  = coeffs[:,0]
    s0 = coeffs[:,1]
    s1 = coeffs[:,2]
    	
    site_term = s0 + s1 / (log10(vs30) - log10(150))
    '''
    site_term = 0
    
    lnsa = mag_term + atten_term + dep_term + near_field_term + site_term
           
    A19imt = {'per':T, 'sa':lnsa}

    return A19imt

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort, isnan
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

from calc_oq_gmpes import inslab_gsims, scr_gsims, get_station_vs30
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from misc_tools import get_mpl2_colourlist

fig = plt.figure(1, figsize=(8,8))
ax = plt.subplot(1, 1, 1)

depths = [20, 50, 100, 150, 200, 300]
#cmap = plt.cm.get_cmap('jet',len(depths))
#cols = (cmap(arange(len(depths))))
cols = get_mpl2_colourlist()
syms = ['o', '^', 's', 'd', 'p', 'v']
mag = 7.3
rhyp = 700.
vs30 = 760.

for i, dep in enumerate(depths):
    A19imt_BS = calc_nac_gmm_spectra(mag, rhyp, dep, vs30, 'BS') # use rrup
    label = '$\mathregular{h_z}$ = '+str(int(round(dep)))+' km'
    if i== 4:
        plt.loglog(A19imt_BS['per'], exp(A19imt_BS['sa']), syms[i], ls='-', lw=1.5, color=cols[i], \
                   ms=8, mec=cols[i], mfc='none', mew=2, markevery=5, label=label)
    else:
        plt.loglog(A19imt_BS['per'], exp(A19imt_BS['sa']), syms[i], ls='-', lw=1.5, color=cols[i], \
                   ms=8, mec=cols[i], mfc='none', mew=2, markevery=4, label=label)
    
plt.xlabel('Period (s)', fontsize=20)
plt.ylabel('Spectral Acceleration (g)', fontsize=20)
    
plt.xlim([0.05, 10])
plt.grid(which='both', color='0.75')

plt.legend(loc=3, fontsize=16, numpoints=1)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('figures/spectra_with_depth.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

