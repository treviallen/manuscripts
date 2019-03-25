import shapefile
from os import path
from mapping_tools import get_field_data, distance
from shakemap_tools import parse_faultdat, make_fault_mesh
from misc_tools import get_binned_stats
import matplotlib.pyplot as plt
import datetime as dt
from numpy import arange, array, delete, ones_like, nan, where, isnan, sqrt, mean, median, std, \
                  loadtxt, log, exp, interp, unique
from mmi_tools import allen_etal_2012_rrup_ipe, allen_etal_2012_rhypo_ipe, \
                      atkinson_wald_ceus_ipe, atkinson_wald_cal_ipe, pgm2mmi_worden12, \
                      parse_usgs_dyfi_geocoded, parse_usgs_dyfi_zip, atkinson_worden_wald14_ceus_ipe, \
                      atkinson_worden_wald14_cal_ipe



maxcontfile = 'max_cont.csv'
contlookup = loadtxt(maxcontfile, delimiter=',')
mw_lu = contlookup[:,0]
rh_lu = contlookup[:,1]

####################################################################################
# now read it back in to save time
####################################################################################
# YYYYMMDDHHMN,MW,EVLO,EVLA,OBLO,OBLA,MMI,REPI,RHYP,DEP,EVENT,CLASS
# parse file here
#lines = open('mmidat_export.csv').readlines()[1:]
lines = open('mmidat_export_moe_clean.csv').readlines()[1:]

repi = []
rhyp = []
rrup = []
eqdep = []
eqdt = []
mmimw = []
mmisc = []
mmi = []
eqname = []

for line in lines:
    dat = line.strip().split(',')
    
    repi.append(float(dat[7]))
    rhyp.append(float(dat[8]))
    rrup.append(float(dat[9]))
    eqdep.append(float(dat[11]))
    eqdt.append(dt.datetime.strptime(dat[0], '%Y%m%d%H%M'))
    mmimw.append(float(dat[1]))
    mmi.append(float(dat[6]))
    eqname.append(dat[-2])
    mmisc.append(dat[-1])
    
####################################################################################
# remove data based on distance cut-off
####################################################################################
fig = plt.figure(10, figsize=(8, 8))
rcutoff = exp(interp(mmimw, mw_lu, log(rh_lu)))
#idx = where(rhyp > rcutoff)[0]

idx = []
i = 0
for rc, rh in zip(rcutoff, rhyp):
    if rh > rc:
        idx.append(i)
    i += 1

mmi = delete(mmi, idx)
rhyp = delete(rhyp, idx)
rrup = delete(rrup, idx)
mmimw = delete(mmimw, idx)
eqdep = delete(eqdep, idx)
repi = delete(repi, idx)
#rrup = rhyp

####################################################################################
# calc ipes  
####################################################################################
AW07ceus = []
AW07cal = []
Aea12rup = []
Aea12hypo = []
AWW14ceus = []
AWW14cal = []
L15 = []

c0 = 3.5
c1 = 1.05
c2 = -1.09
c3 = 1.1

for mi, mw, rr, rh, re, ed in zip(mmi, mmimw, rrup, rhyp, repi, eqdep):
    if mw > 0:
        # do AW07 CEUS
        mmipred, sig = atkinson_wald_ceus_ipe(mw, [rr])
        AW07ceus.append(mmipred[0])
        
        # do AW07 CA
        mmipred, sig = atkinson_wald_cal_ipe(mw, [rr])
        AW07cal.append(mmipred[0])
        
        # do Allen et al 2012 Rrup
        mmipred, Aea12sig = allen_etal_2012_rrup_ipe(mw, rr, ed)
        Aea12rup.append(mmipred)
        
        # do Allen et al 2012 Rhypo
        mmipred, Aea12sig = allen_etal_2012_rhypo_ipe(mw, [rh], ed)
        Aea12hypo.append(mmipred[0])
        
        # do AWW14 CEUS
        mmipred, sig = atkinson_worden_wald14_ceus_ipe(mw, [rh], [re])
        AWW14ceus.append(mmipred[0])
        
        # do AWW14 CA
        mmipred, sig = atkinson_worden_wald14_cal_ipe(mw, [rh])
        AWW14cal.append(mmipred)
        
        # Leonard 15
        L15mmi = c0 + c1 * mw + c2 * log(sqrt(rr**2 + (1+c3*exp(mw-5))**2))
        L15.append(L15mmi)

####################################################################################
# calc model residuals  
####################################################################################

modmmi = [array(AW07ceus).flatten(), \
          array(AW07cal).flatten(), \
          array(Aea12rup).flatten(), \
          array(L15).flatten(), \
          array(AWW14ceus).flatten(), \
          array(AWW14cal).flatten()]
          
modres = [mmi - x  for x in modmmi]

titles = ['AW07ceus', 'AW07cal', 'Aea12rup', 'L15', 'AWW14ceus', 'AWW14cal']

# plt residuals with distance
fig = plt.figure(1, figsize=(14, 10))
xmax = 300
# look at stats
for i, mr in enumerate(modres):
    print titles[i], mean(mr), median(mr), std(mr), len(mr)
    
    ax = plt.subplot(3, 2, i+1)
    
    plt.plot(rrup, mr, '+', c='0.5')
    plt.plot([0, xmax],[0, 0], 'k--')
    
    plt.title(titles[i])
    #plt.xlabel('Rrup')
    plt.ylabel('MMI Residual')
    
    plt.xlim([0, xmax])
    plt.ylim([-3, 3])
    
    # get stats
    bins = arange(10, 361, 10)

    medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, array(rrup), mr)
    
    # plot errors
    plt.errorbar(medx, medbin, yerr=stdbin, fmt='rs', lw=1.5)


plt.savefig('mmi_dist_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()
    
# plt residuals with mag
fig = plt.figure(2, figsize=(14, 10))
# look at stats
for i, mr in enumerate(modres):
    
    plt.subplot(3, 2, i+1)
    
    plt.plot(mmimw, mr, '+', c='0.5')
    plt.plot([2., 7],[0, 0], 'k--')
    
    plt.title(titles[i])
    #plt.xlabel('Rrup')
    plt.ylabel('MMI Residual')
    
    plt.xlim([2., 7])
    plt.ylim([-3, 3])
    
    # get stats
    bins = arange(2.5, 6.8, 0.5)

    medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, array(mmimw), mr)
    
    # plot errors
    plt.errorbar(medx, medbin, yerr=stdbin, fmt='rs', lw=1.5)

plt.savefig('mmi_mw_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()

# plt residuals with mag < 100 km
idx = where(array(rrup) <= 100.)[0]
fig = plt.figure(3, figsize=(14, 10))
# look at stats
for i, mr in enumerate(modres):
    
    plt.subplot(3, 2, i+1)
    
    plt.plot(array(mmimw)[idx], mr[idx], '+', c='0.5')
    plt.plot([2., 7],[0, 0], 'k--')
    
    plt.title(titles[i])
    #plt.xlabel('Rrup')
    plt.ylabel('MMI Residual (LT 100 km)')
    
    plt.xlim([2., 7])
    plt.ylim([-3, 3])
    
    # get stats
    bins = arange(2.5, 6.8, 0.5)

    medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, array(mmimw)[idx], mr[idx])
    
    # plot errors
    plt.errorbar(medx, medbin, yerr=stdbin, fmt='rs', lw=1.5)

plt.savefig('mmi_mw_LT100_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()

####################################################################################
# get site classes
####################################################################################

mmisc = array(mmisc)
unqclass = unique(mmisc)

# loop thru SC
n_class = []
for sc in unqclass:
    idx = where(mmisc == sc)[0]
    n_class.append(len(idx))
    
fig = plt.figure(11, figsize=(7,5))
ax = plt.subplot(111)

bar_width = 0.7
h_width = bar_width/2.
x_locs = arange(len(unqclass)) + 1
plt.bar(x_locs-h_width, n_class, bar_width, color='seagreen')
plt.xticks(x_locs)
ax.set_xticklabels(unqclass)

plt.ylabel('Number of Observations')
plt.xlabel('Modified Site Class')


plt.savefig('site_class_bar.png', fmt='png', bbox_inches='tight')

plt.show()