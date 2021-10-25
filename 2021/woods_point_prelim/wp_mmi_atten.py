#moe_mmi_atten.py
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, \
                      atkinson_worden_wald14_cal_ipe, atkinson_worden_wald14_ceus_ipe, \
                      atkinson_worden_wald14_ceus_oq, leonard15_ipe
from misc_tools import get_mpl2_colourlist
from calc_oq_gmpes import scr_gsims, crustal_gsims
from mapping_tools import distance
from shakemap_tools import write_mmi_obs_raw
from misc_tools import get_binned_stats_mean
from numpy import array, arange, log10, exp, logspace, sqrt, isnan, nan, hstack
from scipy.special import erf
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
mpl.style.use('classic')

plt.rcParams['pdf.fonttype'] = 42

fig = plt.figure(figsize=(10, 6))
plt.tick_params(labelsize=12)

# event details
mag = 5.9
eqdep = 12 
eqlat = -37.5063
eqlon = 146.4022
ztor = 4.
maxrrup = 510.
rjb = logspace(0, log10(maxrrup), 60)
rrup = sqrt(rjb**2 + eqdep**2)
rhypo = rrup
repi = rjb
rake = 0. # USGS CMT
dip  = 90.
vs30 = 360.

####################################################################################
# now read and plot DYFI data
dyfifile = 'felt_reports/felt_reports_1km.geojson'

with open(dyfifile) as f:
    data = json.load(f)

dyfimmi  = []
dyfilat = []
dyfilon = []
for feature in data['features']:
    rec = {'geomerty':feature['geometry']['coordinates'][0],
           'centroid':feature['properties']['center']['coordinates'],
           'intensity':feature['properties']['intensityFine'],
           'nresp':feature['properties']['nresp']}
    
    if rec['nresp'] >= 2:
       dyfimmi.append(rec['intensity'])
       dyfilat.append(rec['centroid'][1])
       dyfilon.append(rec['centroid'][0])
       
# calc DYFI hypo dist
dyfilat = array(dyfilat)
dyfilon = array(dyfilon)
dyfirepi = []
for i in range(0, len(dyfilat)):
    dyfirepi.append(distance(eqlat, eqlon, dyfilat[i], dyfilon[i])[0])
dyfirhyp = sqrt(array(dyfirepi)**2 + eqdep**2)

# now plot
d3 = plt.semilogx(dyfirepi, dyfimmi, '+', lw=2., color='0.6', ms=6)
#d3 = plt.semilogx(dyfirepi[0], dyfimmi[0], '+', lw=2., color='0.6', ms=6)

# write DYFI MMI4SM
#write_mmi_obs_raw('DYFI', dyfidict)

####################################################################################
# parse AU ipe
####################################################################################
lines = open('../au_ipe/au_ipe_coefs.dat').readlines()

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

# Allen 2019 deep
A19_deep = c0 * mag + c1 + c2 * log10(sqrt(rrup**2 + rref**2)) + (h1*erf((eqdep-vert)/(h2*sqrt(2))) + h3)

####################################################################################
# plot models   

# do AW07 CEUS
AW07ceus, sig = atkinson_wald_ceus_ipe(mag, rrup)

# do biased AW07 CEUS
#AW07ceus_bias, sig = atkinson_wald_ceus_ipe(mag+SMbias, rrup)

# do AW07 CA
AW07cal, sig = atkinson_wald_cal_ipe(mag, rrup)

# do Allen et al 2012 Rrup
Aea12rup, Aea12sig = allen_etal_2012_rrup_ipe(mag, rrup, eqdep)

# do Allen et al 2012 Rhypo
Aea12hypo, Aea12sig = allen_etal_2012_rhypo_ipe(mag, rrup, eqdep)

# do AWW14 CEUS
SMbias = -0.233
AWW14ceus, sig = atkinson_worden_wald14_ceus_ipe(mag, rhypo, repi)
AWW14ceus_bias = AWW14ceus + SMbias # from Moe

# do AWW14 CA
AWW14cal, sig = atkinson_worden_wald14_cal_ipe(mag, rhypo)

# do Leonard 2015
L15 = leonard15_ipe(mag, rrup)

#AWW14ceus_oq, sig_oq = atkinson_worden_wald14_ceus_oq(mag, rhypo, eqdep) # this works, but not needed

# do Bea14 - W12 PGA
logB14pga = [] # tried pgv, but looked really bad!
"""
for i in range(0, len(rjb)):
    #Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CB14imt, CY14imt \
    #= crustal_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Pea11imt, A12imt, Bea14imt, YA15imt \
    = scr_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    logB14pga.append(Bea14imt['pga'][0][0]) 
"""
# convert ln g to cm/s**2
from scipy.constants import g
#B14pga = exp(logB14pga) * 100. * g
# now get MMI from GMICE
#B14W12mmi, sig = pgm2mmi_worden12(10**array(logB14pgv), 'pga', mag, rrup)  
#B14W12mmiPGA, sig = pgm2mmi_worden12(B14pga, 'pga', mag, rrup) 

# do Bea14 - W12 PGV
logB14pgv = [] # tried pgv, but looked really bad!
"""
for i in range(0, len(rjb)):
    #Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CB14imt, CY14imt \
    #= crustal_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Pea11imt, A12imt, Bea14imt, YA15imt \
    = scr_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    if isnan(Bea14imt['pgv']):
        logB14pgv.append(nan) 
    else:
        logB14pgv.append(Bea14imt['pgv'][0][0]) 
"""
#B14pgv = exp(array(logB14pgv))
## now get MMI from GMICE
#B14W12mmiPGV, sig = pgm2mmi_worden12(B14pgv, 'pgv', mag, rrup)  
#B14W12mmiPGVnan, sig = pgm2mmi_worden12(B14pgv, 'pgv', nan, nan)
 
# now plot models
#SMbias = -0.16
cl = get_mpl2_colourlist()
syms = ['o', '^', 's', 'd', 'v', '<', 'h', '>', 'p']

# make secondary plots to get around color issues
#h1 = plt.semilogx(rjb, AW07ceus, syms[0], color=cl[0], ls='-', ms=7, mec=cl[0], markevery=5)
#h2 = plt.semilogx(rjb, AW07cal,  syms[1], color=cl[1], ls='-', ms=7, mec=cl[1], markevery=5)
h1 = plt.semilogx(rjb, Aea12rup, syms[0], color=cl[0], ls='-', ms=7, mec=cl[0], markevery=5)
h2 = plt.semilogx(rjb, AWW14ceus,syms[1], color=cl[1], ls='-', ms=7, mec=cl[1], markevery=5)
h3 = plt.semilogx(rjb, AWW14cal, syms[2], color=cl[2], ls='-', ms=7, mec=cl[2], markevery=5)
h4 = plt.semilogx(rjb, L15,      syms[3], color=cl[3], ls='-', ms=7, mec=cl[3], markevery=5)
h5 = plt.semilogx(rjb, A19_deep, syms[4], color=cl[4], ls='-', ms=7, mec=cl[4], markevery=5)

#h4 = plt.semilogx(rjb, B14W12mmiPGA, '-s', color=cl[3], lw=1.5, ms=7, mec=cl[3], markevery=5)

##################################################################################

# get binned stats
bins = arange(0.05, log10(maxrrup)+0.05, 0.1)
meanres, stdres, medx, outbins, nperbin = get_binned_stats_mean(bins, log10(dyfirepi), dyfimmi)

d4 = plt.errorbar(10**medx, meanres, yerr=stdres, fmt='ks', ms=7, elinewidth=1.5, capsize=5, zorder=30000) 

leg1 = plt.legend([h1[0], h2[0], h3[0], h4[0], h5[0], d3[0], d4[0]], \
           ['AWW12 ATR', 'AWW14 CEUS', 'AWW14 CA', 'L15 AU', 'A21 AU', 'Raw DYFI', 'Binned MMI'], fontsize=12, loc=3, numpoints=1)

plt.grid(which='both', color='0.5')
plt.xlim([8, maxrrup])
plt.ylim([1, 8])
plt.xlabel('Epicentral Distance (km)', fontsize=14)
plt.ylabel('Macroseismic Intensity', fontsize=14)

xtic = [10, 20, 50, 100, 200]
xlab = ['10', '20', '50', '100', '200']
plt.xticks(xtic, xlab)
ylab = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII']
ytic = range(1, 9)
plt.yticks(ytic, ylab)

plt.savefig('wp_mmi_atten.png', format='png', dpi=300, bbox_inches='tight')
#plt.savefig('figures/moe_mmi_atten.svg', format='svg', dpi=300, bbox_inches='tight')

plt.show()
