#moe_mmi_atten.py
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, atkinson_worden_wald14_ceus_ipe, \
                      atkinson_worden_wald14_cal_ipe
from calc_oq_gmpes import scr_gsims
from mapping_tools import distance
from misc_tools import get_binned_stats
from numpy import array, arange, log10, exp, logspace, sqrt, where, zeros_like, isnan, nan
import matplotlib.pyplot as plt

# event details
mag = 5.15
eqdep = 9.4
eqlat = -38.252
eqlon = 146.234
ztor = 9.4
rrup = logspace(1., log10(300), 50)
rhypo = rrup
repi = sqrt(rhypo**2 - ztor**2)
rake = 90. # USGS CMT
dip  = 30.
vs30 = 360.

# do AW07 CEUS
AW07ceus, sig = atkinson_wald_ceus_ipe(mag, rrup)

# do AW07 CA
AW07cal, sig = atkinson_wald_cal_ipe(mag, rrup)

# do Allen et al 2012 Rrup
Aea12rup, Aea12sig = allen_etal_2012_rrup_ipe(mag, rrup, eqdep)

# do Allen et al 2012 Rhypo
Aea12hypo, Aea12sig = allen_etal_2012_rhypo_ipe(mag, rrup, eqdep)

# do AWW14 CEUS
AWW14ceus, sig = atkinson_worden_wald14_ceus_ipe(mag, rrup, repi)

# do Bea14 - W12 PGA
rjb = sqrt(rrup**2 - eqdep**2)
logB14pga = [] # tried pgv, but looked really bad!
for i in range(0, len(rjb)):
    Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Pea11imt, A12imt, Bea14imt, YA15imt \
    = scr_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    logB14pga.append(Bea14imt['pga'][0][0]) 

# convert ln g to cm/s**2
from scipy.constants import g
B14pga = exp(logB14pga) * 100. * g
# now get MMI from GMICE
#B14W12mmi, sig = pgm2mmi_worden12(10**array(logB14pgv), 'pga', mag, rrup)  
B14W12mmiPGA, sig = pgm2mmi_worden12(B14pga, 'pga', mag, rrup) 

# do Bea14 - W12 PGV
logB14pgv = [] # tried pgv, but looked really bad!
for i in range(0, len(rjb)):
    Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Pea11imt, A12imt, Bea14imt, YA15imt \
    = scr_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    logB14pgv.append(Bea14imt['pgv'][0][0]) 

B14pgv = exp(array(logB14pgv))
# now get MMI from GMICE
B14W12mmiPGV, sig = pgm2mmi_worden12(B14pgv, 'pgv', mag, rrup)  
B14W12mmiPGVnan, sig = pgm2mmi_worden12(B14pgv, 'pgv', nan, nan)
 
# now plot models
fig = plt.figure(figsize=(10, 6))
plt.tick_params(labelsize=12)

plt.semilogx(rhypo, AW07ceus, '-', color='c', lw=2)
plt.semilogx(rhypo, AW07cal, '-', color='g', lw=2)
plt.semilogx(rhypo, Aea12rup, 'b-', lw=2)
#plt.semilogx(rhypo, Aea12hypo, '-', color='orange', lw=2)
plt.semilogx(rhypo, B14W12mmiPGA, '-', color='orange', lw=2)
plt.semilogx(rhypo, B14W12mmiPGV, '-', color='r', lw=2)
#plt.semilogx(rhypo, B14W12mmiPGVnan, '--', color='r', lw=2)

# now read and plot DYFI data
dyfifile = 'usgs_geocoded_cdi.txt'
dyfidict = parse_usgs_dyfi_geocoded(dyfifile)
#dyfifile = 'usgs_zip_cdi.txt'
#dyfidict = parse_usgs_dyfi_zip(dyfifile)
dyfimmi  = []
dyfilat = []
dyfilon = []
for rec in dyfidict:
    if rec['nresp'] >= 1:
       dyfimmi.append(rec['cdi'])
       dyfilat.append(rec['lat'])
       dyfilon.append(rec['lon'])
       
# calc DYFI hypo dist
dyfilat = array(dyfilat)
dyfilon = array(dyfilon)
dyfirepi = []
for i in range(0, len(dyfilat)):
    dyfirepi.append(distance(eqlat, eqlon, dyfilat[i], dyfilon[i])[0])
dyfirhyp = sqrt(array(dyfirepi)**2 + eqdep**2)

# now plot
plt.semilogx(dyfirhyp, dyfimmi, '+', lw=2., color='0.25', ms=9)

# get binned stats
bins = arange(1.1, log10(300)+0.05, 0.1)
medres, stdres = get_binned_stats(bins, log10(dyfirhyp), dyfimmi)
idx=~isnan(array(medres))
plt.errorbar(10**bins[idx], medres[idx], yerr=stdres[idx], fmt='ks', ms=7) 

plt.legend(['AW07 CEUS', 'AW07 CA', 'Aea12 Rrup', 'Bea14-Wea12 (PGA)', 'Bea14-Wea12 (PGV)', 'USGS DYFI (GC)', 'Median MMI'], \
            fontsize=11, loc=1, numpoints=1)

plt.grid(which='both', color='0.5')
plt.xlim([10, 300])
plt.ylim([1, 7])
plt.xlabel('Hypocentral Distance (km)', fontsize=14)
plt.ylabel('Macroseismic Intensity', fontsize=14)


plt.savefig('Moe_MMI_atten.png', format='png', dpi=300, bbox_inches='tight')

plt.show()
