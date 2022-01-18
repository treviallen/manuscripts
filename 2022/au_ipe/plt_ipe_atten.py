import shapefile
from os import path
from mapping_tools import get_field_data, distance
from shakemap_tools import parse_faultdat, make_fault_mesh
from misc_tools import get_binned_stats
import matplotlib.pyplot as plt
import datetime as dt
from scipy.special import erf
from numpy import arange, array, delete, ones_like, nan, where, isnan, sqrt, mean, median, std, \
                  loadtxt, log, exp, interp, unique, logspace, log10
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, atkinson_worden_wald14_ceus_ipe, \
                      atkinson_worden_wald14_cal_ipe
import matplotlib as mpl
mpl.style.use('classic')
import warnings
warnings.filterwarnings("ignore")

# set rhyp
depd = 10.
deps = 5.
repi = logspace(0, log10(500.), 100)
rhypd = sqrt(repi**2 + depd**2)
rhyps = sqrt(repi**2 + deps**2)
rrupd = rhypd
rrups = rhyps

AW07ceus = []
AW07cal = []
Aea12rup = []
Aea12hypo = []
AWW14ceus = []
AWW14cal = []
L15 = []
A19_shallow = []
A19_deep = []
A22_shallow = []
A22_deep = []

# leonard params
lc0 = 3.5
lc1 = 1.05
lc2 = -1.09
lc3 = 1.1

####################################################################################
# parse pre-internet AU ipe
####################################################################################
lines = open('au_ipe_coefs_pre-internet.dat').readlines()

vert = 7.
c0 = float(lines[1].strip().split()[-1]) 
c1 = float(lines[2].strip().split()[-1]) 
c2 = float(lines[3].strip().split()[-1]) 
c3 = float(lines[4].strip().split()[-1]) 
rref = float(lines[5].strip().split()[-1]) 
xh = float(lines[6].strip().split()[-1]) 
h1 = float(lines[7].strip().split()[-1]) 
h2 = float(lines[8].strip().split()[-1]) 
h3 = float(lines[9].strip().split()[-1]) 

####################################################################################
# parse original (2020) AU ipe
####################################################################################
lines = open('au_ipe_coefs_2020.dat').readlines()

vert = 7.
oc0 = float(lines[1].strip().split()[-1]) 
oc1 = float(lines[2].strip().split()[-1]) 
oc2 = float(lines[3].strip().split()[-1]) 
oc3 = float(lines[4].strip().split()[-1]) 
orref = float(lines[5].strip().split()[-1]) 
oxh = float(lines[6].strip().split()[-1]) 
oh1 = float(lines[7].strip().split()[-1]) 
oh2 = float(lines[8].strip().split()[-1]) 
oh3 = float(lines[9].strip().split()[-1]) 

####################################################################################
# calc IPEs
####################################################################################

m = 5.4

for rrd, rrs, rh, re in zip(rrupd, rrups, rhypd, repi):
    if m > 0:
        # do AW07 CEUS
        mmipred, sig = atkinson_wald_ceus_ipe(m, [rrd])
        AW07ceus.append(mmipred[0])
        
        # do AW07 CA
        mmipred, sig = atkinson_wald_cal_ipe(m, [rrd])
        AW07cal.append(mmipred[0])
        
        # do Allen et al 2012 rrdup
        mmipred, Aea12sig = allen_etal_2012_rrup_ipe(m, rrd, 10.)
        Aea12rup.append(mmipred)
        
        # do Allen et al 2012 Rhypo
        mmipred, Aea12sig = allen_etal_2012_rhypo_ipe(m, [rh], 10.)
        Aea12hypo.append(mmipred[0])
        
        # do AWW14 CEUS
        mmipred, sig = atkinson_worden_wald14_ceus_ipe(m, [rh], [re])
        AWW14ceus.append(mmipred[0])
        
        # do AWW14 CA
        mmipred, sig = atkinson_worden_wald14_cal_ipe(m, [rh])
        AWW14cal.append(mmipred)
        
        # Leonard 15
        L15mmi = lc0 + lc1 * m + lc2 * log(sqrt(rrd**2 + (1+lc3*exp(m-5))**2))
        L15.append(L15mmi)
        
        # Allen 2021 shallow
        c3 = 0. # data biased at large distances
        if rrs > xh:
            A22_shallow.append(c0 * m + c1 + c2 * log10(sqrt(rrs**2 + rref**2)) + (c3 * (rrs - xh)) + (h1*erf((deps-vert)/(h2*sqrt(2))) + h3))
        else:
            A22_shallow.append(c0 * m + c1 + c2 * log10(sqrt(rrs**2 + rref**2)) + (h1*erf((deps-vert)/(h2*sqrt(2))) + h3))

        # Allen 2021 deep
        if rrd > xh:
            A22_deep.append(c0 * m + c1 + c2 * log10(sqrt(rrd**2 + rref**2)) + (c3 * (rrd - xh)) + (h1*erf((depd-vert)/(h2*sqrt(2))) + h3))
        else:
            A22_deep.append(c0 * m + c1 + c2 * log10(sqrt(rrd**2 + rref**2)) + (h1*erf((depd-vert)/(h2*sqrt(2))) + h3))
        
        # Allen 2019 shallow
        if rrs > oxh:
            A19_shallow.append(oc0 * m + oc1 + oc2 * log10(sqrt(rrs**2 + orref**2)) + (oc3 * (rrs - oxh)) + (oh1*erf((deps-vert)/(oh2*sqrt(2))) + oh3))
        else:
            A19_shallow.append(oc0 * m + oc1 + oc2 * log10(sqrt(rrs**2 + rref**2)) + (oh1*erf((deps-vert)/(oh2*sqrt(2))) + oh3))

        # Allen 2019 deep
        if rrd > oxh:
            A19_deep.append(oc0 * m + oc1 + oc2 * log10(sqrt(rrd**2 + orref**2)) + (oc3 * (rrd - oxh)) + (oh1*erf((depd-vert)/(oh2*sqrt(2))) + oh3))
        else:
            A19_deep.append(oc0 * m + oc1 + oc2 * log10(sqrt(rrd**2 + rref**2)) + (oh1*erf((depd-vert)/(oh2*sqrt(2))) + oh3))

####################################################################################
# plot IPEs
####################################################################################

from gmt_tools import cpt2colormap
ncols = 9
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'
#cptfile = 'U:\\DATA\\GMT\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)
cs = (cmap(arange(ncols)))

titles = ['AW07ceus', 'AW07cal', 'Aea12rup', 'L15', 'AWW14ceus', 'AWW14cal', ]

# plt residuals with distance
fig = plt.figure(1, figsize=(12, 6.5))
xmax = 250

#plt.semilogx(repi, AW07ceus, lw=1.5, c=cs[0], label='AW07ceus')
#plt.semilogx(repi, AW07cal, lw=1.5, c=cs[1], label='AW07cal')
plt.semilogx(repi, Aea12rup, lw=1.5, c=cs[0], label='Aea12rup')
plt.semilogx(repi, AWW14ceus, lw=1.5, c=cs[1], label='AWW14ceus')
plt.semilogx(repi, AWW14cal, lw=1.5, c=cs[2], label='AWW14cal')
plt.semilogx(repi, L15, lw=1.5, c=cs[3], label='L15')
plt.semilogx(repi, A19_shallow, ls='-.', lw=1.5, c=cs[4], label='A19shallow')
plt.semilogx(repi, A19_deep, ls='-.', lw=1.5, c=cs[5], label='A19deep')

ridxd = where(rrupd < xh)[0]
ridxs = where(rrups < xh)[0]
plt.semilogx(repi, A22_shallow, ls='--', lw=1.5, c=cs[6])
plt.semilogx(repi, A22_deep, ls='--', lw=1.5, c=cs[-1])
plt.semilogx(repi[ridxs], array(A22_shallow)[ridxs], ls='-', lw=1.5, c=cs[6], label='A22shallow')
plt.semilogx(repi[ridxd], array(A22_deep)[ridxd], ls='-', lw=1.5, c=cs[-1], label='A22deep')

plt.grid(which='both', color='0.5')
plt.legend(loc=3)
plt.xlim([5, xmax])
plt.ylim([2, 7.5])

plt.xlabel('Epicentral Distance (km)', fontsize=14)
plt.ylabel('Predicted MMI', fontsize=14)
plt.title('MMI Attenuation for MW '+str(m), fontsize=16)

plt.savefig('mmi_model_cmp_atten.png', fmt='png', bbox_inches='tight')    
plt.show()
    
