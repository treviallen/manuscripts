#moe_mmi_atten.py
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, \
                      atkinson_worden_wald14_cal_ipe, atkinson_worden_wald14_ceus_ipe
from calc_oq_gmpes import scr_gsims, crustal_gsims
from mapping_tools import distance
from shakemap_tools import write_mmi_obs_raw
from misc_tools import get_binned_stats_mean
from numpy import array, arange, log10, exp, logspace, sqrt, isnan, nan, hstack
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

plt.rcParams['pdf.fonttype'] = 42

fig = plt.figure(figsize=(10, 6))
plt.tick_params(labelsize=12)

# event details
mag = 5.0
eqdep = 18.0
eqlat = -38.259
eqlon = 146.290
ztor = eqdep
maxrrup = 410.
rjb = logspace(0, log10(maxrrup), 50)
rrup = sqrt(rjb**2 + eqdep**2)
rhypo = rrup
rake = 90. # USGS CMT
dip  = 30.
vs30 = 360.

####################################################################################
# now read GA MMI data at 0.05d grid
gammifile = 'mmi_data/GA_MMI_0.05d.csv'
lines = open(gammifile).readlines()
gammi = []
galat = []
galon = []
for line in lines:
    dat = line.strip().split(',')
    gammi.append(float(dat[2]))
    galon.append(float(dat[0]))
    galat.append(float(dat[1]))
    
garepi = []
for i in range(0, len(galat)):
    garepi.append(distance(eqlat, eqlon, galat[i], galon[i])[0])
garhyp = sqrt(array(garepi)**2 + eqdep**2)
# now plot
h7 = plt.semilogx(garepi, gammi, 'x', lw=2., color='0.25', ms=7)

####################################################################################
# now read GA MMI data at 0.05d grid
melmmifile = 'mmi_data/MEL_MMI_0.05d.csv'
lines = open(melmmifile).readlines()
melmmi = []
mellat = []
mellon = []
for line in lines:
    dat = line.strip().split(',')
    melmmi.append(float(dat[2]))
    mellon.append(float(dat[0]))
    mellat.append(float(dat[1]))
    
melrepi = []
for i in range(0, len(mellat)):
    melrepi.append(distance(eqlat, eqlon, mellat[i], mellon[i])[0])
melrhyp = sqrt(array(melrepi)**2 + eqdep**2)
# now plot
h8 = plt.semilogx(melrepi, melmmi, 'o', lw=2., markerfacecolor='None', markeredgecolor='0.25', ms=7)

####################################################################################
# now read and plot DYFI data
dyfifile = 'mmi_data/usgs_geocoded_cdi.txt'
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
h9 = plt.semilogx(dyfirepi, dyfimmi, '+', lw=2., color='0.25', ms=9)

# write DYFI MMI4SM
write_mmi_obs_raw('DYFI', dyfidict)

####################################################################################
# plot models   

# do AW07 CEUS
AW07ceus, sig = atkinson_wald_ceus_ipe(mag, rrup)

# do AW07 CA
AW07cal, sig = atkinson_wald_cal_ipe(mag, rrup)

# do Allen et al 2012 Rrup
Aea12rup, Aea12sig = allen_etal_2012_rrup_ipe(mag, rrup, eqdep)

# do Allen et al 2012 Rhypo
Aea12hypo, Aea12sig = allen_etal_2012_rhypo_ipe(mag, rrup, eqdep)

# do AWW 14 Cal
AWW14cal = atkinson_worden_wald14_cal_ipe(mag, rrup, eqdep)

# do Bea14 - W12 PGA
logB14pga = [] # tried pgv, but looked really bad!

for i in range(0, len(rjb)):
    #Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CB14imt, CY14imt \
    #= crustal_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
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
    #Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CB14imt, CY14imt \
    #= crustal_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Pea11imt, A12imt, Bea14imt, YA15imt \
    = scr_gsims(mag, eqdep, ztor, dip, rake, rrup[i], rjb[i], vs30)
    logB14pgv.append(Bea14imt['pgv'][0][0]) 

B14pgv = exp(array(logB14pgv))
# now get MMI from GMICE
B14W12mmiPGV, sig = pgm2mmi_worden12(B14pgv, 'pgv', mag, rrup)  
B14W12mmiPGVnan, sig = pgm2mmi_worden12(B14pgv, 'pgv', nan, nan)
 
# now plot models
SMbias = -0.16

plt.semilogx(rjb, AW07ceus, '-', color='c', lw=2)
h2 = plt.semilogx(rjb, AW07ceus+SMbias, '--', color='c', lw=2)
plt.semilogx(rjb, AW07cal, '-', color='limegreen', lw=2)
plt.semilogx(rjb, Aea12rup, 'b-', lw=2)
#plt.semilogx(rjb, Aea12hypo, '--', color='orange', lw=2)
plt.semilogx(rjb, B14W12mmiPGA, '-', color='orange', lw=2)
plt.semilogx(rjb, B14W12mmiPGV, '-', color='r', lw=2)
#plt.semilogx(rhypo, B14W12mmiPGVnan, '--', color='r', lw=2)

# make secondary plots to get around color issues
idx = range(4,51,5)
h1 = plt.semilogx(rjb[idx], AW07ceus[idx], '-^', color='c', lw=0.25, ms=7, mec='c')
#h2 = plt.semilogx(rjb, AW07ceus+SMbias, '--', color='c', lw=0.25)
h3 = plt.semilogx(rjb[idx], AW07cal[idx], '-v', color='limegreen', lw=0.25, ms=7, mec='limegreen')
h4 = plt.semilogx(rjb[idx], Aea12rup[idx], '-ob', lw=0.25, ms=7, mec='b')
h5 = plt.semilogx(rjb[idx], B14W12mmiPGA[idx], '-s', color='orange', lw=0.25, ms=7, mec='orange')
h6 = plt.semilogx(rjb[idx], B14W12mmiPGV[idx], '-D', color='r', lw=0.25, ms=7, mec='r')

##################################################################################
# stack data
binrhyp = hstack((dyfirhyp, garhyp, melrhyp))
binrepi = hstack((dyfirepi, garepi, melrepi))
binmmi = hstack((dyfimmi, gammi, melmmi))

# get binned stats
bins = arange(0.05, log10(maxrrup)+0.05, 0.1)
meanres, stdres, medx, outbins, nperbin = get_binned_stats_mean(bins, log10(binrepi), binmmi)

h10 = plt.errorbar(10**medx, meanres, yerr=stdres, fmt='ks', ms=7) 

plt.legend([h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0], h8[0], h9[0], h10[0]], \
           ['AW07 CEUS', 'AW07 CEUS + '+r'$\delta$', 'AW07 CA', 'Aea12 Rrup', 'Bea14-Wea12 (PGA)', 'Bea14-Wea12 (PGV)', \
            'AU MMI', 'MEL MMI', 'USGS DYFI', 'Mean MMI'], fontsize=10, loc=1, numpoints=1)

plt.grid(which='both', color='0.5')
plt.xlim([8, maxrrup])
plt.ylim([1, 7])
plt.xlabel('Epicentral Distance (km)', fontsize=14)
plt.ylabel('Macroseismic Intensity', fontsize=14)

xtic = [10, 20, 50, 100, 200]
xlab = ['10', '20', '50', '100', '200']
plt.xticks(xtic, xlab)
ylab = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
ytic = range(1, 8)
plt.yticks(ytic, ylab)

plt.savefig('moe_mmi_atten.png', format='png', dpi=300, bbox_inches='tight')

plt.show()
