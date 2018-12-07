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
    
################################################################################
# import funcs
################################################################################
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors, colorbar
from numpy import where, zeros, exp, log, interp, array, logspace, log10, sqrt, arange
from misc_tools import listdir_extension, dictlist2array, remove_last_cmap_colour
from calc_oq_gmpes import inslab_gsims, scr_gsims
from gmt_tools import cpt2colormap
from os import path
mpl.style.use('classic')

################################################################################
# nget cmap
################################################################################
ncolours=12
cptfile='/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

################################################################################
# now get data in given mag range
################################################################################

deprng = range(50, 600, 50)

Tplt = 0.2 # secs

# compute inslab gmpes
rrup = logspace(2, log10(3000))
mplt = 7.0
fig = plt.figure(1, figsize=(9,9))
ax = plt.subplot(111)

ii = 0    
for i, dep in enumerate(deprng):
        
    gea05 = []
    a17 = []
    for rr in rrup:
        rjb = sqrt(rr**2 - 50**2)
        YYea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, Aea15imt \
            = inslab_gsims(mplt, dep, 40, 30., 90., rr, rjb, 760)
            
        # subduction gmms
        gea05.append(interp(Tplt, Gea05imt['per'], Gea05imt['sa']))
        
        A17imt = calc_nac_gmm_spectra(mplt, rr)
        a17.append(interp(Tplt, A17imt['per'], A17imt['sa']))

    # plt Gea05
    plt.loglog(rrup, exp(gea05), '-', c=cs[i], lw=1.5, label=' '.join(('h =',str(dep),'km')))
    
    # plt A17
    plt.loglog(rrup, exp(a17), '--', c=cs[i], lw=1.0)
    
ax.set_xscale("log")
ax.set_yscale("log")
plt.xlim([100, 3000])
xticks = [100, 200, 500, 1000, 2000]
ax.set_xticks(xticks)
labels = [str('%0.0f' % x) for x in xticks]
ax.set_xticklabels(labels)

plt.ylim([1E-6, 1E0])
plt.grid(which='both')
#plt.title('MW '+str(mplt))

plt.legend(loc=3, fontsize=14)
    
plt.ylabel('SA('+str(Tplt)+') in g', fontsize=16)
plt.xlabel('Hypocental Distance (km)', fontsize=16)
    
plt.savefig('garcia_depth_var.'+str(Tplt)+'.png', fmt='png', bbox_inches='tight')
plt.show()





















