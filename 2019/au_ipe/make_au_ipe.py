from mapping_tools import get_field_data, distance
from misc_tools import get_binned_stats, get_binned_stats_mean
import matplotlib.pyplot as plt
import datetime as dt
from numpy import arange, array, delete, ones_like, nan, where, isnan, sqrt, mean, median, \
                  std, loadtxt, log, exp, unique, logspace, log10, hstack, savetxt, loadtxt, \
                  interp, nanmedian, nanstd, zeros_like, ones_like, isinf
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack

print 'FIX MOE SITE CLASS KLUGE'

####################################################################################
# parse max dist lookup - file made below!
####################################################################################

maxcontfile = 'max_cont.csv'
contlookup = loadtxt(maxcontfile, delimiter=',')
mw_lu = contlookup[:,0]
rh_lu = contlookup[:,1]
	
#mw_lu = 

####################################################################################
# now read it back in to save time
####################################################################################
#YYYYMMDDHHMN,MW,EVLO,EVLA,OBLO,OBLA,MMI,REPI,RHYP,DEP,EVENT,CLASS,GRA

# parse file here
#lines = open('mmidat_export.csv').readlines()[1:]
lines = open('mmi_grad.csv').readlines()[1:]
	

repi = []
rhyp = []
rrup = []
eqdep = []
eqdt = []
mw = []
mmi = []
eqname = []
siteclass = []
grad = []

for line in lines:
    dat = line.strip().split(',')
    
    eqdt.append(dt.datetime.strptime(dat[0], '%Y%m%d%H%M'))
    mw.append(float(dat[1]))
    mmi.append(float(dat[6]))
    repi.append(float(dat[7]))
    rhyp.append(float(dat[8]))
    rrup.append(float(dat[9]))
    eqdep.append(float(dat[11]))
    eqname.append(dat[12])
    siteclass.append(dat[13])
    grad.append(float(dat[-1]))
    

repi = array(repi)
mmi = array(mmi)
rhyp = array(rhyp)
rrup = array(rrup) # initial assumption
mw = array(mw)
siteclass = array(siteclass)
grad = array(grad)
eqdep = array(eqdep)

####################################################################################
# remove data based on distance cut-off
####################################################################################
fig = plt.figure(10, figsize=(8, 8))
rcutoff = exp(interp(mw, mw_lu, log(rh_lu)))
#idx = where(rhyp > rcutoff)[0]

idx = []
i = 0
for rc, rh in zip(rcutoff, rrup):
    if rh > rc:
        idx.append(i)
    i += 1

mmi = delete(mmi, idx)
rhyp = delete(rhyp, idx)
rrup = delete(rrup, idx)
mw = delete(mw, idx)
siteclass = delete(siteclass, idx)
grad = delete(grad, idx)
eqdep = delete(eqdep, idx)

#rrup = rhyp

plt.semilogx(rrup, mw, 'b+')

plt.show()

####################################################################################
# get GR for each mag bin
####################################################################################

fig = plt.figure(1, figsize=(14, 10))

ueqname = unique(array(eqname))
mbins = arange(3.5, 6.8, 0.25) 
halfbin = 0.25
logxmax = 2.

slopes = []
meanmws = []

for i, mb in enumerate(mbins):

    idx = where((mw >= (mb - halfbin)) & (mw < (mb + halfbin)))[0]
    
    # now, get stats
    bins = arange(1., logxmax, 0.1)

    medbin, stdbin, medx, binstrp, nperbin = get_binned_stats_mean(bins, log10(rrup[idx]), mmi[idx])
    
    # strip bins with LT 2 obs
    delidx = where(nperbin < 2)[0]
    medbin = delete(medbin, delidx)
    binstrp = delete(binstrp, delidx)
    stdbin = delete(stdbin, delidx)
    
    if len(binstrp) > 2:
        # now regress bins
        c = linregress(binstrp, medbin)
        
        ax = plt.subplot(4, 4, i+1)
        
        plt.plot(log10(rrup[idx]), mmi[idx], '+', c='0.5')
        plt.plot(binstrp, medbin, 'rs')
        
        xfit = array([1, logxmax])
        yfit = c[0] * xfit + c[1]
        
        plt.plot(xfit, yfit, 'k-', lw=2.0)
        
        plt.xlim([0.3, 2.5])
        plt.ylim([2, 9])
        
        plt.title(mb)
        
        slopes.append(c[0])
        meanmws.append(mean(mw[idx]))
        
    else:
        slopes.append(nan)
        meanmws.append(nan)
        
plt.suptitle('Variable (MW-dependent) slope')
meanslope = nanmedian(slopes)
print 'getting median slope'
print 'cut data GT MMI3 contour'

plt.savefig('ipe_var_slope_GR.png', fmt='png', bbox_inches='tight') 


plt.show()

####################################################################################
# get mag dependent slopes
####################################################################################

idx = ~isnan(meanmws)

plt.plot(meanmws, slopes, 'o')

ms = linregress(array(meanmws)[idx], array(slopes)[idx])

msl = array([2.5, 7.0])
mgr = ms[0] * msl + ms[1]
plt.plot(msl, mgr, 'k-', lw=2.)

plt.xlabel('Magnitude')
plt.ylabel('Greometric Spreading')
plt.show()


####################################################################################
# fix mean slope & get intercept
####################################################################################

def linear_fixed_slope(c, x):
        
    return meanslope * x + c[0]

def near_src_trunc_fixed_slope(c, x):
    from numpy import sqrt
    xx = sqrt(rref*2 + x**2)
    
    #mb = x
    #print mb
    # mb = bin mag
    #if mb < 3.75:
    #    mb = 3.75
    
    slope = ms[0] * max([mb, 3.75]) + ms[1]
    
    return slope * log10(xx) + c[0]

rref = 5.


fig = plt.figure(2, figsize=(14, 10))

intercepts = []
meanmagbin = []
meanmws = []

mbins = arange(3., 6.8, 0.25) # use smaller mag to constrain m-scaling
xmax = 350.
for i, mb in enumerate(mbins):
    
    idx = where((mw >= (mb - halfbin)) & (mw < (mb + halfbin)) & (rrup < xmax))[0]
    
    # now, get stats
    bins = arange(0.7, logxmax, 0.1)

    medbin, stdbin, medx, binstrp, nperbin = get_binned_stats_mean(bins, log10(rrup[idx]), mmi[idx])
    
    # strip bins with LT 2 obs
    delidx = where(nperbin < 2)[0]
    medbin = delete(medbin, delidx)
    binstrp = delete(binstrp, delidx)
    stdbin = delete(stdbin, delidx)
    
    # get intercept
    data = odrpack.RealData(binstrp, medbin)
    intfit = odrpack.Model(linear_fixed_slope)
    odr = odrpack.ODR(data, intfit, beta0=[5.0])
    
    #data = odrpack.RealData(10**binstrp, medbin, meta={'mag':mb})
    #intfit = odrpack.Model(near_src_trunc_fixed_slope)
    #odr = odrpack.ODR(data, intfit, beta0=[5.0])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    intercept = out.beta
    intercepts.append(intercept[0])
    meanmagbin.append(mean(mmi[idx]))
    meanmws.append(mean(mw[idx]))
    
    if len(binstrp) > 0:
        # now regress bins
        c = linregress(binstrp, medbin)
        
        ax = plt.subplot(4, 4, i+1)
        
        plt.plot(log10(rrup[idx]), mmi[idx], '+', c='0.5')
        plt.plot(binstrp, medbin, 'rs')
        
        xfit = arange(0.3, logxmax, 0.1)
        #yfit = meanslope * xfit + intercept
        
        # assume mean slope
        yfit = meanslope * log10(sqrt((10**xfit)**2 + rref**2)) + intercept
        
        # assume mag-dept slope
        slope = ms[0] * mb + ms[1]
        #yfit = slope * log10(sqrt((10**xfit)**2 + rref**2)) + intercept
        
        plt.plot(xfit, yfit, 'k-', lw=2.0)
        
        plt.xlim([0.3, 2.5])
        plt.ylim([2, 9])
        
        plt.title(mb)
        plt.suptitle('Fixed Slope')
        
plt.savefig('ipe_fixed_slope_GR.png', fmt='png', bbox_inches='tight') 
plt.show()

####################################################################################
# plt intercept with mw
####################################################################################

fig = plt.figure(2, figsize=(8, 8))

idx = isnan(meanmws) == False
plt.plot(meanmws, intercepts, 'rs')
#plt.plot(meanmws[:-2], intercepts[:-2], 'bs')

plt.xlabel('Mean Bin Mag')
plt.ylabel('MMI Intercept')

# regress mags
#m = linregress(meanmws[:-2], intercepts[:-2])
m = linregress(array(meanmws)[idx], array(intercepts)[idx])
xfit = array([2.5, 6.7])
yfit = m[0] * xfit + m[1]
plt.plot(xfit, yfit, 'k-', lw=2.0)

plt.savefig('ipe_mag_scaling.png', fmt='png', bbox_inches='tight') 
plt.show()

####################################################################################
# try calculating residual
####################################################################################

# get magscaling
mf = m[0] * mw + m[1]

# get mag-dept GR
magslope = zeros_like(mw)
for i in range(0, len(mw)):
    magslope[i] = ms[0] * max([mw[i], 3.75]) + ms[1]

# get prelim distance scaling 
print 'Using MEAN slope...'
rf = meanslope * log10(sqrt(rrup**2 + rref**2))
#rf = magslope * log10(sqrt(rrup**2 + rref**2))

# get pred mmi
pmmi = mf + rf

# get mmi residual
rmmi = mmi - pmmi

####################################################################################
# get 3.5 MMI contour for data clipping
####################################################################################
mrng = arange(2.5, 6.8, 0.1)

'''
mmi35 = 3.

cont35 = []
for mr in mrng:
    
    mcont = m[0] * mr + m[1]
    
    cont35.append(sqrt((10**((mmi35 - mcont) / meanslope))**2 - rref**2))
    
fig = plt.figure(4, figsize=(6, 6))
plt.semilogy(mrng, cont35, 'rs')
plt.ylabel('MMI IV Contour')
plt.xlabel('Mangnitude')
plt.grid(which='major')

# export look-up
cont35 = array(cont35)
cont35 = cont35.reshape(len(cont35), 1)
mrng = mrng.reshape(len(mrng), 1)
#savetxt('max_cont.csv' , hstack((mrng, cont35)), delimiter=',')#
'''

####################################################################################
# plot dist residual
####################################################################################

# plt residuals with distance
fig = plt.figure(3, figsize=(10, 4))
xmax = 350
# look at stats

plt.plot(rrup, rmmi, '+', c='0.5')
plt.plot([0, xmax],[0, 0], 'k--')

plt.title('Prelim Atten Res')
#plt.xlabel('Rrup')
plt.ylabel('MMI Residual')

plt.xlim([0, xmax])
plt.ylim([-3, 3])

# get stats
bins = arange(5, xmax+1, 10)

medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, rrup, rmmi)

# plot errors
plt.errorbar(binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

# annotate
didx = where(array(rrup) <= 300)[0] # comparable with other mod cmp
pltstr = 'Med = ' + str('%0.2f' % median(rmmi[didx])) + '\n' \
         'Std = ' + str('%0.2f' % std(rmmi[didx]))
plt.text(280, -2.35, pltstr, fontsize=14, va='center')

plt.savefig('prelim_mmi_residuals.png', fmt='png', bbox_inches='tight')  

y_at_xhinge = 0
xhinge = 105.
def fit_fixed_intercept(c, x):
    '''
    x = array fo x values
    y_at_xmax = y value at xmax
    '''
    
    xmax = 10**logxmax # hardwired in distance atten
    
    ans = c * (x - xhinge) + y_at_xhinge
    
    return ans


# fit distance atten
ridx = where((binstrp >= xhinge) & (binstrp < xmax))[0]
data = odrpack.RealData(binstrp[ridx], medbin[ridx])
print binstrp[ridx]
truncdist = odrpack.Model(fit_fixed_intercept)
odr = odrpack.ODR(data, truncdist, beta0=[0.])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, ODR = 0
out = odr.run()

c = out.beta[0] # slope of dist atten

####################################################################################
# try calculating residual with distcne correction
####################################################################################

ridx = where(rrup >= xhinge)
dc = zeros_like(rrup)
dc[ridx] = c * (rrup[ridx] - xhinge)

# get pred mmi
pmmi2 = mf + rf + dc

# get mmi residual
rmmi2 = mmi - pmmi2

####################################################################################
# re-plot dist residual
####################################################################################

# plt residuals with distance
fig = plt.figure(5, figsize=(10, 4))
#xmax = 400
# look at stats

plt.plot(rrup, rmmi2, '+', c='0.5')
plt.plot([0, xmax],[0, 0], 'k--')

plt.title('Prelim Atten Res - Long Dist Corr')
#plt.xlabel('Rrup')
plt.ylabel('MMI Residual')

plt.xlim([0, xmax])
plt.ylim([-3, 3])

# get stats
#bins = arange(10, 361, 10)

medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, rrup, rmmi2)

# plot errors
plt.errorbar(binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

# annotate
didx = where(array(rrup) <= 300)[0] # comparable with other mod cmp
pltstr = 'Med = ' + str('%0.2f' % median(rmmi2[didx])) + '\n' \
         'Std = ' + str('%0.2f' % std(rmmi2[didx]))
plt.text(280, -2.35, pltstr, fontsize=14, va='center')

plt.savefig('prelim_mmi_distcorr_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()
  
####################################################################################
# plt residual by predicted MMI BY depth
####################################################################################
    
fig = plt.figure(31, figsize=(10,5))
ax = plt.subplot(111)

didx = where(array(rrup) <= 300)[0]

plt.plot(eqdep[didx], rmmi2[didx], '+', c='0.5')
plt.plot([0, 20], [0, 0], 'k--')
plt.ylim([-3, 3])
plt.xlim([0, 20])

# get stats
bins = arange(0, 21, 1)
medbin, stdbin, medx, binstrp, nperbin = get_binned_stats_mean(bins, eqdep[didx], rmmi2[didx])

# plot errors
plt.errorbar(binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

plt.ylabel('MMI Residual')
plt.xlabel('Depth (km)')


####################################################################################
# get depth corrections
####################################################################################

from scipy.special import erf
from scipy.optimize import curve_fit

vert = 7.
def erfunc(x, mFL, b, c):
    from numpy import sqrt
    
    return mFL*erf((x-vert)/(b*sqrt(2))) + c

# guesses mFL, a, b, c = 0.5, 30, 10, 0
ef, extras = curve_fit(erfunc, medx, medbin, p0=[0., 10, 0])

xerf = arange(0, 20, 0.05)
yerf = ef[0]*erf((xerf-vert)/(ef[1]*sqrt(2))) + ef[2]
plt.plot(xerf, yerf, 'k-', lw=2.)


####################################################################################
# apply depth corrections
####################################################################################

hc = ef[0]*erf((eqdep-vert)/(ef[1]*sqrt(2))) + ef[2]
# get pred mmi
pmmi3 = mf + rf + dc + hc

# get mmi residual
rmmi3 = mmi - pmmi3

# annotate
didx = where(array(rrup) <= 300)[0] # comparable with other mod cmp
pltstr = 'Med = ' + str('%0.2f' % median(rmmi3[didx])) + '\n' \
         'Std = ' + str('%0.2f' % std(rmmi3[didx]))
plt.text(17, -2.35, pltstr, fontsize=14, va='center')

plt.savefig('au_ipe_depth_res.png', fmt='png', bbox_inches='tight')

plt.show()


####################################################################################
# export IPE coefs
####################################################################################
# c2 = ms[0] * mw + ms[1]
coefs  = 'mmi = c0 * mw + c1 + c2 * log10(sqrt(rrup**2 + rref**2)) { + c3 * (rhyp[ridx] - xh)} + h1*erf((rhyp-vert)/(h2*sqrt(2))) + h3\n'
coefs += 'c0 ' + str(m[0]) + '\n'
coefs += 'c1 ' + str(m[1]) + '\n'
coefs += 'c2 ' + str(ms[0]) + '\n'
coefs += 'c2 ' + str(ms[1]) + '\n'
coefs += 'c3 ' + str(c) + '\n'
coefs += 'rref ' + str(5.0) + '\n'
coefs += 'xh ' + str(50.0) + '\n'
coefs += 'h1 ' + str(ef[0]) + '\n'
coefs += 'h2 ' + str(ef[1]) + '\n'
coefs += 'h3 ' + str(ef[2])

f = open('au_ipe_coefs.dat', 'wb')
f.write(coefs)
f.close()

  
####################################################################################
# plot mag residuals
####################################################################################
    

# plt residuals with mag
fig = plt.figure(2, figsize=(10, 4))
# look at stats
plt.plot(mw[didx], rmmi3[didx], '+', c='0.5')
plt.plot([2., 7],[0, 0], 'k--')

plt.title('Prelim Mag Res')
#plt.xlabel('Rrup')
plt.ylabel('MMI Residual')

plt.xlim([2., 7])
plt.ylim([-3, 3])

# get stats
bins = arange(2.5, 6.8, 0.5)

medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, mw[didx], rmmi3[didx])

# plot errors
plt.errorbar(binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

plt.savefig('prelim_mmi_mw_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()

# plt residuals with mag
fig = plt.figure(2, figsize=(10, 4))
idx = where(array(rrup) <= 100.)[0]

# look at stats
plt.plot(mw[idx], rmmi3[idx], '+', c='0.5')
plt.plot([2., 7],[0, 0], 'k--')

plt.title('Prelim Mag Res')
#plt.xlabel('Rrup')
plt.ylabel('MMI Residual')

plt.xlim([2., 7])
plt.ylim([-3, 3])

# get stats
bins = arange(2.5, 6.8, 0.5)

medbin, stdbin, medx, binstrp, npb = get_binned_stats(bins, mw[idx], rmmi3[idx])

# plot errors
plt.errorbar(binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

plt.savefig('prelim_mmi_mw_LT100_residuals.png', fmt='png', bbox_inches='tight')    
plt.show()

####################################################################################
# GET RESIDUAL BY SITE CLASS
####################################################################################

unqclass = unique(siteclass)

# loop thru SC
n_class = []
med_class = []
std_class = []
x_plt = ones_like(mmi)
xmax = 300

for i, sc in enumerate(unqclass):
    
    idx = where((siteclass == sc) & (rhyp <= xmax))[0]
    n_class.append(len(idx))
    print sc, len(idx)
    med_class.append(nanmedian(rmmi3[idx]))
    std_class.append(nanstd(rmmi3[idx]))
    x_plt[idx] = x_plt[idx] * (i+1)
    
fig = plt.figure(11, figsize=(7,5))
ax = plt.subplot(111)

x_locs = arange(len(unqclass)) + 1
plt.plot(x_plt[didx], rmmi3[didx], '+', c='0.5')
plt.plot([0, 9], [0, 0], 'k--')
plt.ylim([-3, 3])
plt.xlim([0, len(unqclass)+1])

# plot errors
plt.errorbar(x_locs, med_class, yerr=std_class, fmt='rs', lw=1.5)

plt.xticks(x_locs)
ax.set_xticklabels(unqclass)

plt.ylabel('MMI Residual')
plt.xlabel('Modified Site Class')

plt.savefig('site_class_res.png', fmt='png', bbox_inches='tight')

plt.show()

####################################################################################
# plt residual by predicted MMI BY SITE CLASS
####################################################################################

fig = plt.figure(12, figsize=(14, 10))

for i, sc in enumerate(unqclass):
    ax = plt.subplot(4,2,i+1)
    
    med_class = []
    std_class = []
    
    idx = where((siteclass == sc) & (rhyp <= xmax))[0]
    
    plt.plot([2, 9], [0, 0], 'k--')
    plt.plot(pmmi3[idx], rmmi3[idx], '+', c='0.5')
    
    # get stats
    # now, get stats
    bins = arange(2, 10, 1)

    medbin, stdbin, medx, binstrp, nperbin = get_binned_stats_mean(bins, pmmi3[idx], rmmi3[idx])
    
    # plot errors
    plt.errorbar(medx, medbin, yerr=stdbin, fmt='rs', lw=1.5)

    plt.title('Site Class '+sc)
    plt.ylabel('MMI Residual')
    plt.ylim([-3, 3])

    if i == 6:
        plt.xlabel('Predicted MMI')
    
plt.savefig('site_class_with_MMI.png', fmt='png', bbox_inches='tight')
    
plt.show()

####################################################################################
# plt residual by predicted MMI BY grdient
####################################################################################


# loop thru SC
'''
n_class = []
med_class = []
std_class = []
x_plt = ones_like(mmi)

for i, sc in enumerate(unqclass):
    
    idx = where((siteclass == sc) & (rhyp <= xmax))[0]
    n_class.append(len(idx))
    print sc, len(idx)
    med_class.append(nanmedian(rmmi2[idx]))
    std_class.append(nanstd(rmmi2[idx]))
    x_plt[idx] = x_plt[idx] * (i+1)
'''
    
fig = plt.figure(11, figsize=(10,5))
ax = plt.subplot(111)

plt.semilogx(grad, rmmi3, '+', c='0.5')
plt.plot([0.0001, 1], [0, 0], 'k--')
plt.ylim([-3, 3])
plt.xlim([0.0001, 1])

# get stats
logbins = arange(-4, 1, 0.1)
didx = where(array(rrup) <= 300)[0]
medbin, stdbin, medx, binstrp, nperbin = get_binned_stats_mean(logbins, log10(grad[didx]), rmmi3[didx])

# plot errors
plt.errorbar(10**binstrp, medbin, yerr=stdbin, fmt='rs', lw=1.5)

# fit linear
idx = ~isnan(binstrp >= -3.35)
gradfit = linregress(array(binstrp)[idx], array(medbin)[idx])

gradx = array([-3.35, log10(0.25)])
grady = gradfit[0] * gradx + gradfit[1]
plt.plot(10**gradx, grady, 'k-', lw=2.)

####################################################################################
# apply gradient corrections
####################################################################################

gc = gradfit[0] * log10(grad) + gradfit[1]
# get pred mmi
pmmi4 = mf + rf + dc + hc + gc

# get mmi residual
rmmi4 = mmi - pmmi4

didx = where((array(rhyp) <= 300) & (~isinf(rmmi4)))[0] # comparable with other mod cmp
pltstr = 'Med = ' + str('%0.2f' % nanmedian(rmmi4[didx])) + '\n' \
         'Std = ' + str('%0.2f' % nanstd(rmmi4[didx]))
plt.text(0.1, -2., pltstr, fontsize=14, va='center')


plt.ylabel('MMI Residual')
plt.xlabel('Gradient (m/m)')

plt.savefig('au_ipe_grad_res.png', fmt='png', bbox_inches='tight')

plt.show()

