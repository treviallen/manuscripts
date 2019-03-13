from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import array, arange, log10, where, median, std, mean, concatenate, vstack, interp
from misc_tools import get_binned_stats_mean

plt.close('all')
# for bilinear sloping segment

def set_legend(plt, leglist, loc):
	plt.legend(leglist, loc=loc)
	leg = plt.gca().get_legend()
	ltext  = leg.get_texts()
	plt.setp(ltext, fontsize='small')


'''
plot_fault_len_wid.py

'''

# use historical interface data
usehistoric = True
plt_inter_only = False

import pickle
from shakemap_tools import *
from mapping_tools import distance, reckon
from numpy import array, sqrt, nan, isnan, arange, abs, unique, hstack, savetxt, \
                  logical_and, mean, median, std, log10, ones, logspace, exp, max, signbit
#from make_slab_fault import make_slab_matrix
from os import path, sep
from fault_tools import *
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from scipy.stats import linregress
from mag_tools import mw2m0
from misc_tools import plttext
from gmt_tools import cpt2colormap
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 
coeftxt = ''

# now load trimmed faults
ftxt = open('..//Data_Prep//trimmed_faults_dip_type.csv').readlines()[1:]
tdate = []
tlon = []
tlat = []
tdep = []
tlen = []
twid = []
tmag = []
treg = []
tdip = []
ttyp = []
tavs = []
tmxs = []
tare = []
for line in ftxt:
    dat = line.strip().split(',')
    tdate.append(dat[0])
    tmag.append(float(dat[3]))
    tlon.append((dat[4]))
    tlat.append((dat[5]))
    tdep.append(float(dat[6]))
    tlen.append(float(dat[9]))
    twid.append(float(dat[10]))
    tare.append(float(dat[11]))
    tmxs.append(float(dat[12]))
    tavs.append(float(dat[13]))
    treg.append(dat[16])
    tdip.append(float(dat[17]))
    ttyp.append(dat[-1])

inter_fwid  = []
inter_fmag = []
inter_flen  = []
inter_hdep = []
inter_ztor = []
inter_zbor = []
inter_deprat = []
inter_sdip = []
inter_sreg = []
inter_farea = []
inter_favs = []
inter_fmxs = []
inter_farea = []

intra_fwid  = []
intra_fmag = []
intra_flen  = []
intra_hdep = []
intra_ztor = []
intra_zbor = []
intra_deprat = []
intra_farea = []
intra_sdip = []
intra_sreg = []
intra_favs = []
intra_fmxs = []
intra_farea = []

outer_fwid  = []
outer_fmag = []
outer_flen  = []
outer_hdep = []
outer_ztor = []
outer_zbor = []
outer_deprat = []
outer_farea = []
outer_sdip = []
outer_sreg = []
outer_favs = []
outer_fmxs = []
outer_farea = []

trans_fwid  = []
trans_fmag = []
trans_flen  = []
trans_hdep = []
trans_ztor = []
trans_zbor = []
trans_deprat = []
trans_farea = []
trans_sdip = []
trans_sreg = []
trans_favs = []
trans_fmxs = []
trans_farea = []

hist_fwid  = []
hist_fmag = []
hist_flen  = []
hist_hdep = []
hist_ztor = []
hist_zbor = []
hist_deprat = []
hist_farea = []
hist_sdip = []
hist_sreg = []
hist_favs = []
hist_fmxs = []
hist_farea = []

tabtxt = ''
for i in range(0,len(ttyp)):
    if ttyp[i] == 'i':# and treg[i] != 'MEX':
        inter_flen.append(tlen[i])
        inter_fwid.append(twid[i])
        inter_fmag.append(tmag[i])
        inter_sdip.append(tdip[i])
        inter_sreg.append(treg[i])
        inter_farea.append(tare[i])
        inter_favs.append(tavs[i])
        inter_fmxs.append(tmxs[i]) 
    elif ttyp[i] == 's':
        intra_flen.append(tlen[i])
        intra_fwid.append(twid[i])
        intra_fmag.append(tmag[i])
        intra_sdip.append(tdip[i])
        intra_sreg.append(treg[i])
        intra_farea.append(tare[i])
        intra_favs.append(tavs[i])
        intra_fmxs.append(tmxs[i])
    elif ttyp[i] == 'o':
        outer_flen.append(tlen[i])
        outer_fwid.append(twid[i])
        outer_fmag.append(tmag[i])
        outer_sdip.append(tdip[i])
        outer_sreg.append(treg[i])
        outer_farea.append(tare[i])
        outer_favs.append(tavs[i])
        outer_fmxs.append(tmxs[i])
    elif ttyp[i] == 'h':
        hist_flen.append(tlen[i])
        hist_fwid.append(twid[i])
        hist_fmag.append(tmag[i])
        hist_sdip.append(tdip[i])
        hist_sreg.append(treg[i])
        hist_farea.append(tare[i])
        hist_favs.append(tavs[i])
        hist_fmxs.append(tmxs[i])
    elif ttyp[i] == 't':
        trans_flen.append(tlen[i])
        trans_fwid.append(twid[i])
        trans_fmag.append(tmag[i])
        trans_sdip.append(tdip[i])
        trans_sreg.append(treg[i])
        trans_farea.append(tare[i])
        trans_favs.append(tavs[i])
        trans_fmxs.append(tmxs[i])

    # get text for summary table
    tabtxt += '\t'.join((tdate[i], treg[i], tlon[i], tlat[i], str(tdep[i]), \
                       str(tmag[i]), ttyp[i], str(tlen[i]), str(twid[i]), \
                       str(tare[i]), str(tmxs[i]), str(tavs[i]), str(tdip[i]))) + '\n'

    # use historical data
    if usehistoric == True:
        reg_flen = concatenate((inter_flen, hist_flen))
        reg_fwid = concatenate((inter_fwid, hist_fwid))
        reg_fmag = concatenate((inter_fmag, hist_fmag))
        reg_farea = concatenate((inter_farea, hist_farea))
        reg_sdip = concatenate((inter_sdip, hist_sdip))
        reg_sreg = concatenate((inter_sreg, hist_sreg))
        reg_favs = concatenate((inter_favs, hist_favs))
        reg_fmxs = concatenate((inter_fmxs, hist_fmxs))
    else:
        reg_flen = inter_flen
        reg_fwid = inter_fwid
        reg_fmag = inter_fmag
        reg_farea = inter_farea
        reg_sdip = inter_sdip
        reg_sreg = inter_sreg
        reg_favs = inter_favs
        reg_fmxs = inter_fmxs
        print array(reg_sdip)

# write table
header =    '\t'.join(('ID', 'REGION', 'LON', 'LAT', 'DEP', 'MAG', \
                       'SZTYPE', 'LEN', 'WID', 'AREA', 'MXS', 'AVS', 'DIP')) + '\n'
                       
tabtxt = header + tabtxt
f = open('summary_trimmed_table.txt', 'wb')
f.write(tabtxt)
f.close()

'''
########################################################################################
'''
bins = arange(7.0, 9.7, 0.2)
meanflen, stdlen, meanmagl, newbins = get_binned_stats_mean(bins, reg_fmag, log10(reg_flen))
meanfwid, stdwid, meanmagw, newbins = get_binned_stats_mean(bins, reg_fmag, log10(reg_fwid))
meanfarea, stdarea, meanmaga, newbins = get_binned_stats_mean(bins, reg_fmag, log10(reg_farea))
meanfavs, stdavs, meanmagav, newbins = get_binned_stats_mean(bins, reg_fmag, log10(reg_favs))
meanfmxs, stdmxs, meanmagmx, newbins = get_binned_stats_mean(bins, reg_fmag, log10(reg_fmxs))
'''
########################################################################################
'''
coeftxt = 'Function,a,b,SEa,SEb,sig,Condition\n'

# make uncertainty arrays for weighting for ODR
reg_logMo = (mw2m0(array(reg_fmag)))

#sx = ones(len(reg_fmag))*0.1 # uniform uncertainty of 0.2 mu
#sx = ones(len(reg_fmag))*0.2
im = [7.5, 8.0]
isi = [0.2, 0.1]
sx = interp(reg_fmag, im, isi)
#idx = reg_fmag >= 8.0
#sx[idx] = 0.1
sy = sx

# do ODR on raw data
func  = models.polynomial(1)
#data = Data(array(inter_fmag),log10(array(inter_flen)))
data = odrpack.RealData(array(reg_fmag),log10(array(reg_flen)), sx=sx, sy=sy)

# set data
mrng = arange(7., 9.61, 0.01)

'''
########################################################################################
test linear vs ODR linear vs ODR bi-linear
'''

# regress M vs L - ODR
beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
odrlen_lin = 10**(a + b * mrng)

linodrlenres = (a + b * array(reg_fmag)) - log10(array(reg_flen))
'''
# regress M vs L - OLS
beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
olslen_lin = 10**(a + b * mrng)

linolslenres = (a + b * array(reg_fmag)) - log10(array(reg_flen))

print '\nLinear Residuals', std(linodrlenres), std(linolslenres), '\n'
'''
# export coeffs
savetxt('len1_coeffs.txt', [a, b], delimiter='\t')

coeftxt += ','.join(('log10 L1 = a + b x MW',str('%0.2f' % a),str('%0.2f' % b), \
                     str('%0.2f' % out.sd_beta[0]),str('%0.2f' % out.sd_beta[1]), \
                     str('%0.2f' % sqrt(out.res_var)))) + '\n'
                     
'''
########################################################################################
get linear M vs W
'''

# regress M vs W - ODR
data = odrpack.RealData(array(reg_fmag),log10(array(reg_fwid)), sx=sx, sy=sy)

beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
odrwid_lin = 10**(a + b * mrng)

linodrwidres = (a + b * array(reg_fmag)) - log10(array(reg_fwid))

print '\nLinear Residuals', std(linodrwidres), std(linodrwidres), '\n'

# export coeffs
savetxt('wid1_coeffs.txt', [a, b], delimiter='\t')

coeftxt += ','.join(('log10 W1 = a + b x MW',str('%0.2f' % a),str('%0.2f' % b), \
                     str('%0.2f' % out.sd_beta[0]),str('%0.2f' % out.sd_beta[1]), \
                     str('%0.2f' % sqrt(out.res_var)))) + '\n'

'''
########################################################################################
get linear M vs A
'''

# regress M vs W - ODR
data = odrpack.RealData(array(reg_fmag),log10(array(reg_farea)), sx=sx, sy=sy)

beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
odrarea_lin = 10**(a + b * mrng)

linodrareares = (a + b * array(reg_fmag)) - log10(array(reg_farea))

print '\nLinear Residuals', std(linodrareares), std(linodrareares), '\n'

# export coeffs
savetxt('area1_coeffs.txt', [a, b], delimiter='\t')

coeftxt += ','.join(('log10 S1 = a + b x MW',str('%0.2f' % a),str('%0.2f' % b), \
                     str('%0.2f' % out.sd_beta[0]),str('%0.2f' % out.sd_beta[1]), \
                     str('%0.2f' % sqrt(out.res_var)))) + '\n'


'''
########################################################################################
do max slip vs mag
''' 
if usehistoric == True:
    idx = ~isnan(reg_fmxs)
    data = odrpack.RealData(array(reg_fmag[idx]),log10(array(reg_fmxs[idx])), sx=sx, sy=sy) 
else:                                                                                      
    data = odrpack.RealData(array(reg_fmag),log10(array(reg_fmxs)), sx=sx, sy=sy)

# regress M vs Max slip
beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
odrmxs = 10**(a + b * mrng)

# export coeffs
savetxt('mxs_coeffs.txt', [a, b], delimiter='\t')

coeftxt += ','.join(('log10 Dmax = a + b x MW',str('%0.2f' % a),str('%0.2f' % b), \
                     str('%0.3f' % out.sd_beta[0]),str('%0.3f' % out.sd_beta[1]), \
                     str('%0.3f' % sqrt(out.res_var)))) + '\n'


'''
do average slip vs mag
'''
if usehistoric == True:
    idx = ~isnan(reg_favs)
    data = odrpack.RealData(array(reg_fmag[idx]),log10(array(reg_favs[idx])), sx=sx, sy=sy) 
else:                                                                                      
    data = odrpack.RealData(array(reg_fmag),log10(array(reg_favs)), sx=sx, sy=sy)

# regress M vs Max slip
beta0=[-2.7, 0.62]
odr = ODR(data,func,beta0)
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
#out.pprint()
print ' '
a = out.beta[0]
b = out.beta[1]
odravs = 10**(a + b * mrng)

# export coeffs
savetxt('avs_coeffs.txt', [a, b], delimiter='\t')

coeftxt += ','.join(('log10 Dav = a + b x MW',str('%0.3f' % a),str('%0.3f' % b), \
                     str('%0.3f' % out.sd_beta[0]),str('%0.3f' % out.sd_beta[1]), \
                     str('%0.3f' % sqrt(out.res_var)))) + '\n'

'''
########################################################################################
'''
# bi-linear function for area and width!
def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x >= hx
    xmod[idx] = 1
    return xmod
    
def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)
    
    return ans1 + ans2
    
def bilinear_reg_fix_hinge(c, x):
    from numpy import zeros_like
    hx = 9. # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)
    
    return ans1 + ans2
    
def bilinear_reg_fix_slope(c, x):
    from numpy import zeros_like
    hx = c[2] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * yarea
    
    return ans1 + ans2   


def bilinear_reg_fix_slope_level(c, x):
    from numpy import log10, zeros_like
    ylevel = ywid # log value already !!!!! THIS NEEDS TO BE FIXED !!!!!!
    ans = zeros_like(x)
    
    '''
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    '''
    ans = c[0] * x + c[1]
    idx = ans > ylevel
    ans[idx] = ylevel

    return ans 
'''    
def bilinear_reg_fix_slope_hinge(c, x):
    from numpy import zeros_like
    hx = 8.56 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * yarea
    
    return ans1 + ans2
'''
'''
########################################################################################
'''
# regress len vs mag
data2 = odrpack.RealData(array(reg_fmag),log10(array(reg_flen)), sx=sx, sy=sy)
#data2 = odrpack.RealData(array(inter_fmag),log10(array(inter_fwid)))

bilin = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data2, bilin, beta0=[1.0, -3.0, 1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
c=0
d= 0
a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3]
odr_bl_len = 10**(b + a * mrng)
ylen = b + a * hx
idx = mrng > hx
odr_bl_len[idx] = 10**(c * (mrng[idx]-hx) + ylen)
print '\nL2max 9.5', 10**(c * (9.5-hx) + ylen),'\n'

savetxt('len_coeffs.txt', [a, b, c, hx], delimiter='\t')

coeftxt += ','.join(('log10 L2 = a + b x MW',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW <= '+str('%0.3f' % hx))) + '\n'

dy = ylen - (c * hx)
coeftxt += ','.join(('log10 L2 = a + b x MW',str('%0.3f' % dy),str('%0.3f' % c), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[2]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW > '+str('%0.3f' % hx))) + '\n'
                     
odr_bl_lenpred = 10**(b + a * reg_fmag)
ylen = b + a * hx
idx = reg_fmag > hx
odr_bl_lenpred[idx] = 10**(c * (reg_fmag[idx]-hx) + ylen)

odr_bl_lenres = log10(reg_flen) - log10(odr_bl_lenpred)

#print '\nLinear Residuals', std(linodrlenres), std(linolslenres), std(odr_bl_lenres),'\n'


'''
########################################################################################
'''
# regress wid
data2 = odrpack.RealData(array(reg_fmag),log10(array(reg_fwid)), sx=sx, sy=sy)

bilin = odrpack.Model(bilinear_reg_fix_slope)
#bilin = odrpack.Model(bilinear_reg_fix_slope_hinge)
odr = odrpack.ODR(data2, bilin, beta0=[0.6, -1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
c=0
d= 0
a = out.beta[0]
b = out.beta[1]
#c = out.beta[2]
hmag = out.beta[2]
#hx = 8.56
blwid = 10**(b + a * mrng)
ywid = b + a * hmag
idx = mrng > hmag
blwid[idx] = 10**(ywid) 
savetxt('bilin_wid_coeffs.txt', [a, b, hmag], delimiter='\t') 

coeftxt += ','.join(('log10 W2 = a + b x MW',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW <= '+str('%0.3f' % hmag))) + '\n'

dy = ywid
coeftxt += ','.join(('log10 W2 = a + b x MW',str('%0.3f' % dy),'0.00', \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[2]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW > '+str('%0.3f' % hmag))) + '\n'

'''
########################################################################################
'''
# regress area vs mag
data2 = odrpack.RealData(array(reg_fmag),log10(array(reg_farea)), sx=sx, sy=sy)
#data2 = odrpack.RealData(array(inter_fmag),log10(array(inter_fwid)))

bilin = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data2, bilin, beta0=[1.0, -3.0, 1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
c=0
d= 0
a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3]
odr_bl_area = 10**(b + a * mrng)
yarea = b + a * hx
idx = mrng > hx
odr_bl_area[idx] = 10**(c * (mrng[idx]-hx) + yarea)
print '\nS2 Mh', 10**yarea
print 'S2max 9.5', 10**(c * (9.5-hx) + yarea),'\n'

savetxt('area_coeffs.txt', [a, b, c, hx], delimiter='\t')

coeftxt += ','.join(('log10 S2 = a + b x MW',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW <= '+str('%0.3f' % hx))) + '\n'

dy = yarea - (c * hx)
coeftxt += ','.join(('log10 S2 = a + b x MW',str('%0.3f' % dy),str('%0.3f' % c), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[2]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW > '+str('%0.3f' % hx))) + '\n'

'''
########################################################################################
'''
# regress width vs length
data2 = odrpack.RealData(log10(array(reg_flen)),log10(array(reg_fwid)), sx=sx, sy=sy)
#data2 = odrpack.RealData(array(inter_fmag),log10(array(inter_fwid)))

bilin = odrpack.Model(bilinear_reg_fix_slope_level)
odr = odrpack.ODR(data2, bilin, beta0=[1.0, -3.0])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
#out.pprint()
a = out.beta[0]
b = out.beta[1]

lenplt = logspace(log10(20),log10(1500),300)
odr_bl_widlen = 10**(b + a * log10(lenplt))

maxwid = ywid
idx = odr_bl_widlen > 10**maxwid
odr_bl_widlen[idx] = 10**maxwid  
savetxt('LW_coeffs.txt', [a, b, ywid], delimiter='\t')

# get max L for max W
maxlen = (ywid - b) / a
print '\nLW2max 9.5', maxlen,'\n'


coeftxt += ','.join(('log10 W = a + log10 L x b',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'L <= '+str('%0.0f' % 10**maxlen))) + '\n'

coeftxt += ','.join(('log10 W = a + log10 L x b',str('%0.3f' % maxwid),'0.00', \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'L > '+str('%0.0f' % 10**maxlen))) + '\n'

'''
########################################################################################
'''
"""
# regress avs vs mag
data2 = odrpack.RealData(array(reg_fmag),log10(array(reg_favs)), sx=sx, sy=sy)
#data2 = odrpack.RealData(array(inter_fmag),log10(array(inter_fwid)))

bilin = odrpack.Model(bilinear_reg_fix_hinge)
odr = odrpack.ODR(data2, bilin, beta0=[1.0, -3.0, 1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
#out.pprint()
c=0
d= 0
a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3]
odr_bl_avs = 10**(b + a * mrng)
yavs = b + a * hx
print odr_bl_avs
print hx
idx = mrng > hx
print idx
odr_bl_avs[idx] = 10**(c * (mrng[idx]-hx) + yavs)
savetxt('avs_coeffs.txt', [a, b, c, hx], delimiter='\t')

coeftxt += ','.join(('log10 Dav = a + b x MW',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW <= '+str('%0.3f' % hx))) + '\n'

dy = yavs - (c * hx)
coeftxt += ','.join(('log10 Dav = a + b x MW',str('%0.3f' % dy),str('%0.3f' % c), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[2]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW > '+str('%0.3f' % hx))) + '\n'

'''
########################################################################################
'''
# regress mxs vs mag
data2 = odrpack.RealData(array(reg_fmag),log10(array(reg_fmxs)), sx=sx, sy=sy)
#data2 = odrpack.RealData(array(inter_fmag),log10(array(inter_fwid)))

bilin = odrpack.Model(bilinear_reg_fix_hinge)
odr = odrpack.ODR(data2, bilin, beta0=[1.0, -3.0, 1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
out = odr.run()
#out.pprint()
c=0
d= 0
a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3]
odr_bl_mxs = 10**(b + a * mrng)
print odr_bl_mxs
print hx
ymxs = b + a * hx
idx = mrng > hx
odr_bl_mxs[idx] = 10**(c * (mrng[idx]-hx) + ymxs)
savetxt('mxs_coeffs.txt', [a, b, c, hx], delimiter='\t')

coeftxt += ','.join(('log10 Dmx = a + b x MW',str('%0.3f' % b),str('%0.3f' % a), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[0]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW <= '+str('%0.3f' % hx))) + '\n'

dy = ymxs - (c * hx)
coeftxt += ','.join(('log10 Dmx = a + b x MW',str('%0.3f' % dy),str('%0.3f' % c), \
                     str('%0.3f' % out.sd_beta[1]),str('%0.3f' % out.sd_beta[2]), \
                     str('%0.3f' % sqrt(out.res_var)), 'MW > '+str('%0.3f' % hx))) + '\n'
"""
'''
########################################################################################
'''
# export coeffs
f = open('interface_coeffs.csv', 'wb')
f.write(coeftxt)
f.close()

'''
########################################################################################
'''


#from gmt_tools import cpt2colormap

#cptfile = 'U:\\DATA\\GMT\\cpt\\GMT_no_green.cpt'
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/GMT_no_green.cpt'
cmap = cpt2colormap(cptfile, 5)[0] # colormap for data
cs = (cmap(arange(5)))

cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/gay-flag-1979.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/GMT_rainbow.cpt'
ncolours = 7
cmap2, zvals = cpt2colormap(cptfile, ncolours, rev=False)
cs2 = (cmap2(arange(ncolours)))
#cs2 = vstack(([cs2[0:4], cs2[5:]]))


'''
########################################################################################
'''
# plot length
fig = plt.figure(1, figsize=(11.5,17))
ax = plt.subplot(3, 2, 1, aspect=1)
plt.axis('equal')

if plt_inter_only == False:
   plt.semilogy(intra_fmag, intra_flen, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)
   plt.semilogy(outer_fmag, outer_flen, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)
   plt.semilogy(trans_fmag, trans_flen, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)
   #plt.errorbar(meanmag, meanlen, stdlen, 'ks')
plt.semilogy(hist_fmag, hist_flen,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.semilogy(inter_fmag, inter_flen, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)
#plt.semilogy(meanmagl, 10**meanflen, 'ks', ms=7)#, mfc='none')
   
#plot Strasser
wclen = mag2ruplen_WC94(mrng, 'rs')
strasserlen = mag2srl_St10inter(mrng)[0]
blaserlen = mag2len_Bl10rev(mrng)[0]
leonardlen = mag2len_L10(mrng, 'rs')

a0 = plt.semilogy(mrng, wclen, '-', lw=2.0, color=[cs2[0][0],cs2[0][1],cs2[0][2]])
a1 = plt.semilogy(mrng, blaserlen, '-', lw=2.0, color=[cs2[1][0],cs2[1][1],cs2[1][2]])
a6 = plt.semilogy(mrng, leonardlen, '-', lw=2.0, color=[cs2[2][0],cs2[2][1],cs2[2][2]])
a2 = plt.semilogy(mrng, strasserlen, '-', lw=2.0, color=[cs2[3][0],cs2[3][1],cs2[3][2]])
#a3 = plt.semilogy(mrng, odr_bl_len, '-', lw=2.0, color='k') # bi-linear
a3 = plt.semilogy(mrng, odrlen_lin, '-', lw=2.0, color='k') # linear

plt.ylabel('Rupture Length (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([20, 2000])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.07, ymax, '(A)', fontsize=18)

#leg.get_frame().set_alpha(0.)               
'''
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small')


'''
########################################################################################
''''''
# plot width
ax = plt.subplot(3, 2, 2, aspect=1)
plt.axis('equal')
if plt_inter_only == False:
   h1 = plt.semilogy(intra_fmag, intra_fwid, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)      
   h2 = plt.semilogy(outer_fmag, outer_fwid, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)      
   h3 = plt.semilogy(trans_fmag, trans_fwid, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)      
h4 = plt.semilogy(hist_fmag, hist_fwid,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
h0 = plt.semilogy(inter_fmag, inter_fwid, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)

wcwid = mag2rupwid_WC94(mrng, 'rs')
strasserwid = mag2wid_St10inter(mrng)[0]
blaserwid = mag2wid_Bl10rev(mrng)[0]
somerville_l = mag2wid1_So15inter(mrng)[0]
somerville_bl = mag2wid2_So15inter(mrng)[0]
leonardwid = mag2wid_L10(mrng, 'rs')

plt.semilogy(mrng, wcwid, '-', lw=2.0, color=[cs2[0][0],cs2[0][1],cs2[0][2]])
plt.semilogy(mrng, blaserwid, '-', lw=2.0, color=[cs2[1][0],cs2[1][1],cs2[1][2]])
plt.semilogy(mrng, leonardwid, '-', lw=2.0, color=[cs2[2][0],cs2[2][1],cs2[2][2]])
plt.semilogy(mrng, strasserwid, '-', lw=2.0, color=[cs2[3][0],cs2[3][1],cs2[3][2]])
plt.semilogy(mrng, somerville_l, '-', lw=2.0, color=[cs2[5][0],cs2[5][1],cs2[5][2]])
plt.semilogy(mrng, somerville_bl, '--', lw=2.0, color=[cs2[5][0],cs2[5][1],cs2[5][2]])
#plt.semilogy(mrng, odrwid, '-', color='r', lw=2.0)
plt.semilogy(mrng, blwid, '-', lw=2.0, color='k')
plt.semilogy(mrng, odrwid_lin, '--', lw=2.0, color='k') # linear
plt.ylabel('Rupture Width (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([10, 300])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.07, ymax, '(B)', fontsize=18)

# make legend 
if plt_inter_only == False:
    leg = plt.legend((h0[0], h1[0], h2[0], h3[0], h4[0]), ['Interface', \
             'Intraslab','Outer-Rise','Strike-Slip','Other Interface'],loc=4,numpoints=1,fontsize=13)
else:
    leg = plt.legend((h0[0], h4[0]), ['Interface', \
             'Other'],loc=4,numpoints=1,fontsize=10)
             

'''
########################################################################################
'''
# plot M x A
ax = plt.subplot(3, 2, 3)
wcarea = mag2area_WC94(mrng, 'rs')
strasserarea = mag2area_St10inter(mrng)[0]
murotaniarea = mag2area_Mu13inter(mrng)[0]
somervillearea = mag2area_So15inter_SS(mrng)[0]
leonardarea = mag2area_L10(mrng, 'rs')

if plt_inter_only == False:
    plt.semilogy(intra_fmag, intra_farea, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)      
    plt.semilogy(outer_fmag, outer_farea, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)      
    plt.semilogy(trans_fmag, trans_farea, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)      
plt.semilogy(hist_fmag, hist_farea,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.semilogy(inter_fmag, inter_farea, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)      

plt.semilogy(mrng, wcarea, '-', lw=2.0, color=[cs2[0][0],cs2[0][1],cs2[0][2]])
plt.semilogy(mrng, leonardarea, '-', lw=2.0, color=[cs2[2][0],cs2[2][1],cs2[2][2]])
plt.semilogy(mrng, strasserarea, '-', lw=2.0, color=[cs2[3][0],cs2[3][1],cs2[3][2]])
a4 = plt.semilogy(mrng, murotaniarea, '-', lw=2.0, color=[cs2[4][0],cs2[4][1],cs2[4][2]])
a5 = plt.semilogy(mrng, somervillearea, '-', lw=2.0, color=[cs2[5][0],cs2[5][1],cs2[5][2]])
# plot ODR area
#plt.semilogy(mrng, ordarea, 'r-', lw=2.0)
plt.semilogy(mrng, odr_bl_area, '-', lw=2.0, color='k')
plt.semilogy(mrng, odrarea_lin, '--', lw=2.0, color='k') # linear

plt.ylabel(r'Rupture Area $\mathregular{(km^{2}}$)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([500,5E5])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.07, ymax, '(C)', fontsize=18)


'''
########################################################################################
'''
# plot Max slip x M
ax = plt.subplot(3, 2, 4)
if plt_inter_only == False:
   plt.semilogy(intra_fmag, intra_fmxs, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)      
   plt.semilogy(outer_fmag, outer_fmxs, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)      
   plt.semilogy(trans_fmag, trans_fmxs, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)      
plt.semilogy(hist_fmag, hist_fmxs,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.semilogy(inter_fmag, inter_fmxs, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)      

wcmaxs = mag2maxs_WC94(mrng, 'all')[0] # RS is not well constrained
somervillemaxs = mag2maxs_So15inter_SS(mrng)[0]

plt.semilogy(mrng, wcmaxs, '-', lw=2.0, color=[cs2[0][0],cs2[0][1],cs2[0][2]])
plt.semilogy(mrng, somervillemaxs, '-', lw=2.0, color=[cs2[5][0],cs2[5][1],cs2[5][2]])
#plt.semilogy(mrng, odr_bl_mxs, '-', lw=2.0, color='k')
plt.semilogy(mrng, odrmxs, '-', lw=2.0, color='k')

plt.ylabel('Maximum Slip (m)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([.5, 100])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.07, ymax, '(D)', fontsize=18)


'''
########################################################################################
'''
# plot Average slip x M
ax = plt.subplot(3, 2, 5)
if plt_inter_only == False:
    plt.semilogy(intra_fmag, intra_favs, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)      
    plt.semilogy(outer_fmag, outer_favs, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)      
    plt.semilogy(trans_fmag, trans_favs, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)      
plt.semilogy(hist_fmag, hist_favs,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.semilogy(inter_fmag, inter_favs, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)  

murotiniavs = mag2avs_Mu13inter(mrng)[0]
wcavs = mag2avs_WC94(mrng, 'all')[0] # RS is not well constrained
somervilleavs = mag2avs_So15inter_SS(mrng)[0]
leonardavs = mag2avs_L10(mrng, 'rs')


plt.semilogy(mrng, wcavs, '-', lw=2.0, color=[cs2[0][0],cs2[0][1],cs2[0][2]])
plt.semilogy(mrng, murotiniavs, '-', lw=2.0, color=[cs2[4][0],cs2[4][1],cs2[4][2]])
plt.semilogy(mrng, somervilleavs, '-', lw=2.0, color=[cs2[5][0],cs2[5][1],cs2[5][2]])
a6 = plt.semilogy(mrng, leonardavs, '-', lw=2.0, color=[cs2[2][0],cs2[2][1],cs2[2][2]])
plt.semilogy(mrng, odravs, '-', lw=2.0, color='k')
#plt.semilogy(mrng, odr_bl_avs, '-', lw=2.0, color='k')

plt.ylabel('Average Slip (m)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([.1, 40])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.07, ymax, '(E)', fontsize=18)


'''
########################################################################################
'''
# plot L x W
ax = plt.subplot(3, 2, 6)
if plt_inter_only == False:
    plt.loglog(intra_flen, intra_fwid, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)      
    plt.loglog(outer_flen, outer_fwid, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)      
    plt.loglog(trans_flen, trans_fwid, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)      
plt.loglog(hist_flen, hist_fwid,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.loglog(inter_flen, inter_fwid, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)  

#a5 = murotiniavs = mag2avs_Mu13inter(mrng)[0]
#wcavs = mag2avs_WC94(mrng, 'all')[0] # RS is not well constrained
#somervilleavs = mag2area_So15inter(mrng)[0]

leonardlw = len2wid_L10(lenplt, 'rs')
lenplt = logspace(log10(20),log10(1500),300)
plt.loglog(lenplt, leonardlw, '-', lw=2.0, color=[cs2[2][0],cs2[2][1],cs2[2][2]])

#plt.semilogy(mrng, wcavs, 'b-', lw=2.0)
#plt.semilogy(mrng, murotiniavs, '-', color='purple', lw=2.0)
#plt.semilogy(mrng, somervilleavs, '-', color='aqua', lw=2.0)
plt.loglog(lenplt, odr_bl_widlen, '-', lw=2.0, color='k')
plt.ylim([10, 300])
plt.xlim([20, 2000])
oto = plt.loglog([10, 2000], [10, 2000], '--', color='0.45', lw=2, label='1:1')
plt.legend(loc=4, fontsize=13, numpoints=3)


plt.ylabel('Width (km)', fontsize=14)
plt.xlabel('Length (km)', fontsize=14)
#plt.legend(oto[0], '1:1', loc=4, fontsize=12)

ylims = log10(ax.get_ylim())
xlims = log10(ax.get_xlim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
xmin = 10**(0.025*(xlims[1]-xlims[0]) + xlims[0])

plt.text(xmin, ymax, '(F)', fontsize=18)


# make legend for other studies
ax = plt.subplot(3, 2, 1, aspect=1)
'''
leg = plt.legend((a0[0], a1[0], a6[0], a2[0], a4[0], a5[0], a3[0]), ['Wells & Coppersmith (1994)', \
           'Blaser et al (2010)', 'Leonard (2010)', 'Strasser et al (2010)', 'Murotani et al (2013)', \
           'Somerville et al (2015)', 'Present Study'],loc=4,fontsize=10)
'''
leg = plt.legend((a0[0], a1[0], a6[0], a2[0], a4[0], a5[0], a3[0]), ['WC94', \
           'Bea10', 'L10', 'Sea10', 'Mea13', \
           'Sea16', 'AH16'],loc=4,fontsize=13)

'''
# plot L x W
ax = plt.subplot(3, 2, 4)
plt.loglog(inter_flen, inter_fwid, 'bo')
plt.loglog(intra_flen, intra_fwid, 'ro')
plt.loglog(outer_flen, outer_fwid, 'go')
plt.semilogy(hist_flen, hist_fwid, '^', color='orange', ms=7)

plt.loglog(odrlen, widrng, 'r-', lw=2.0)
plt.ylabel('Width (km)')
plt.xlabel('Length (km)')
plt.ylim([10, 500])
plt.xlim([10, 2000])
'''
if usehistoric == True:
    #plt.savefig('misc_lw_type_hist.pdf',dpi=300,format='pdf', transparent=True)
    plt.savefig('misc_lw_type_hist.png',dpi=300,format='png',bbox_inches='tight')
    plt.savefig('misc_lw_type_hist.pdf',dpi=300,format='pdf',bbox_inches='tight')
else:
    #plt.savefig('misc_lw_type.pdf',dpi=300,format='pdf', transparent=True)
    plt.savefig('misc_lw_type.pdf',dpi=300,format='pdf',bbox_inches='tight')
plt.show()


'''
########################################################################################
# examine interface dip res
'''

# make cmap
ticks = arange(5, 31, 5)
cmap=plt.cm.bone_r
cmap=plt.cm.Greys
colstep = int(round(cmap.N / 5.))
cmaplist = [cmap(i) for i in arange(0, cmap.N, colstep)]
cmap = cmap.from_list('custom', cmaplist, cmap.N)
norm = mpl.colors.BoundaryNorm(hstack((ticks, [35])), cmap.N)

fig = plt.figure(2, figsize=(17,4))
ax = plt.subplot(1, 3, 1, aspect=1)
plt.scatter(reg_fmag, reg_flen, c=reg_sdip, s=55, cmap=cmap, norm=norm)
#plot Strasser
'''
h0 = plt.semilogy(mrng, wclen, 'b-', lw=2.0)
h1 = plt.semilogy(mrng, blaserlen, 'g-', lw=2.0)
h2 = plt.semilogy(mrng, strasserlen, '-', color='orange', lw=2.0)
'''
h3 = plt.semilogy(mrng, odr_bl_len, 'k-', lw=2.0)
cb = plt.colorbar(extend='both', ticks=ticks, boundaries=hstack((ticks, [35])))
cb.set_label('Average Dip')
plt.ylabel('Rupture Length (km)')
plt.xlabel('Magnitude')
plt.ylim([20, 2000])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plttext(ax, 7.07, ymax, '(a)', fsize='x-large')

# do legend
#plt.legend((h0[0], h1[0], h2[0], h3[0]), ['WC94','Bea10','Sea10','AH15'],loc=4)

# plot width
ax = plt.subplot(1, 3, 2)
plt.scatter(reg_fmag, reg_fwid, c=reg_sdip, s=55, cmap=cmap, norm=norm)
strasserwid = mag2wid_St10inter(mrng)[0]
blaserwid = mag2wid_Bl10rev(mrng)[0]
'''
plt.semilogy(mrng, wcwid, 'b-', lw=2.0)
plt.semilogy(mrng, blaserwid, 'g-', lw=2.0)
plt.semilogy(mrng, strasserwid, '-', color='orange', lw=2.0)
plt.semilogy(mrng, odrwid, 'r-', lw=2.0)
'''
plt.semilogy(mrng, blwid, '-', color='k', lw=2.0)
plt.ylim([10,300])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plttext(ax, 7.07, ymax, '(b)', fsize='x-large')

cb = plt.colorbar(extend='both', ticks=ticks, boundaries=hstack((ticks, [35])))
cb.set_label('Average Dip')
plt.ylabel('Rupture Width (km)')
plt.xlabel('Magnitude')
plt.show()

############################################################################
# plot area
plt.figure(3, figsize=(5.75,4.5))
ax = plt.subplot(1, 1, 1)

#plt.scatter(reg_fmag, reg_farea, c=reg_sdip, s=55, cmap=cmap, vmin=10, vmax=30, norm=norm)
plt.scatter(reg_fmag, reg_farea, c=reg_sdip, s=55, cmap=cmap, norm=norm)
plt.semilogy(mrng, odr_bl_area, 'k-', lw=2.0)
cb = plt.colorbar(extend='both', ticks=ticks, boundaries=hstack((ticks, [35])))
cb.set_label('Average Dip (Degrees)', fontsize=15)
plt.ylabel(r'Rupture Area $\mathregular{(km^2}$)', fontsize=16)
plt.xlabel('Magnitude', fontsize=16)
plt.ylim([1000,5E5])
plt.xlim([7.0, 9.6])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
#plttext(ax, 7.07, ymax, '(c)', fsize='x-large')
print '\nInterface Magnitudes:', min(reg_fmag), max(reg_fmag), '\n'

'''
# plot L x W - colour by mag
ax = plt.subplot(2, 2, 4)
lenmerge = hstack((reg_flen,intra_flen,outer_flen))
widmerge = hstack((reg_fwid,intra_fwid,outer_fwid))
magmerge = hstack((reg_fmag,intra_fmag,outer_fmag))
plt.scatter(lenmerge, widmerge, c=magmerge, s=40, cmap=plt.cm.Spectral_r)
plt.loglog(odrlen, widrng, 'r-')
plt.ylabel('Width (km)')
plt.xlabel('Length (km)')
cb = plt.colorbar()
cb.set_label('Magnitude')
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylim([10, 500])
plt.xlim([10, 2000])
'''
'''
# bin length
binwid = 0.2
binstart = 6.65
binstop = max(magmerge) + binwid/2.
pltx = []
medplt = []
meanplt = []
stdplt = []
#lenmerge = array(inter_flen) # this is a tmp fix
#magmerge = array(inter_fmag) # this is a tmp fix
for inc in arange(binstart, binstop, binwid):
    eind = logical_and(magmerge >= inc, magmerge < inc+binwid)
    pltx.append(median(magmerge[eind]))
    medplt.append(median(log10(lenmerge[eind])))
    meanplt.append(median(log10(lenmerge[eind])))
    stdplt.append(std(log10(lenmerge[eind])))

ax = plt.subplot(2, 2, 3)
plt.errorbar(pltx,medplt,yerr=stdplt,fmt='rs')
h0 = plt.plot(mrng, log10(wclen), 'b-', lw=2.0)
h1 = plt.plot(mrng, log10(blaserlen), 'g-', lw=2.0)
h2 = plt.plot(mrng, log10(strasserlen), '-', color='orange', lw=2.0)
plt.plot(mrng, log10(odrlen), 'k-', lw=2.0)
plt.ylabel('log Fault Length (km)')
plt.xlabel('Magnitude')


# bin width
pltx = []
medplt = []
meanplt = []
stdplt = []
#widmerge = array(inter_fwid) # this is a tmp fix
for inc in arange(binstart, binstop, binwid):
    eind = logical_and(magmerge >= inc, magmerge < inc+binwid)
    pltx.append(median(magmerge[eind]))
    medplt.append(median(log10(widmerge[eind])))
    meanplt.append(median(log10(widmerge[eind])))
    stdplt.append(std(log10(widmerge[eind])))

ax = plt.subplot(2, 2, 4)
plt.errorbar(pltx,medplt,yerr=stdplt,fmt='rs')
h0 = plt.plot(mrng, log10(wcwid), 'b-', lw=2.0)
h1 = plt.plot(mrng, log10(blaserwid), 'g-', lw=2.0)
h2 = plt.plot(mrng, log10(strasserwid), '-', color='orange', lw=2.0)
plt.plot(mrng, log10(odrwid), 'k-', lw=2.0)
plt.ylabel('log Fault Width (km)')
plt.xlabel('Magnitude')
plt.ylim([10, 1000])
'''
if usehistoric ==True:
    plt.savefig('lw_dip_hist.pdf',dpi=200,format='pdf', bbox_inches='tight')
    plt.savefig('lw_dip_hist.png',dpi=200,format='png', bbox_inches='tight')
else:
    plt.savefig('lw_dip.pdf',dpi=200,format='pdf', bbox_inches='tight')
plt.show()

'''
# plot L x W
ax = plt.subplot(2, 2, 3)
#plt.semilogx(inter_flen, inter_sdip, 'bo')
ax.scatter(inter_flen, inter_sdip, c=inter_fmag, s=40, cmap=plt.cm.Spectral_r)
cb = plt.colorbar()
cb.set_label('Magnitude')
ax.set_xscale('log')
plt.ylabel('Slab Dip')
plt.xlabel('Fault Length (km)')

'''

'''
########################################################################################
'''
"""
fig = plt.figure(3, figsize=(16,5.5))
ax = plt.subplot(1, 2, 1, aspect=1)
plt.scatter(reg_fmag, reg_flen, c=reg_sdip, s=70, cmap=plt.cm.Spectral_r, vmin=10, vmax=30)
#plot Strasser
h0 = plt.semilogy(mrng, wclen, 'b-', lw=2.0)
h1 = plt.semilogy(mrng, blaserlen, 'g-', lw=2.0)
h2 = plt.semilogy(mrng, strasserlen, '-', color='orange', lw=2.0)
h3 = plt.semilogy(mrng, odrlen, 'r-')
cb = plt.colorbar()
cb.set_label('Dip at Epicentre')
plt.ylabel('Fault Length (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([10, 2000])
plt.xlim([7.0, 9.6])
# do legend
leg = plt.legend((h0[0], h1[0], h2[0], h3[0]), ['Wells & Coppersmith (1994)', \
           'Blaser et al (2010)','Strasser et al (2010)','Present Study'],loc=4, fontsize=11)

'''
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small')
'''
plt.title('Interface Length', fontsize=16)

# plot width
ax = plt.subplot(1, 2, 2)
plt.scatter(reg_fmag, reg_fwid, c=reg_sdip, s=70, cmap=plt.cm.Spectral_r, vmin=10, vmax=30)
strasserwid = mag2wid_St10inter(mrng)[0]
blaserwid = mag2wid_Bl10rev(mrng)[0]
plt.semilogy(mrng, wcwid, 'b-', lw=2.0)
plt.semilogy(mrng, blaserwid, 'g-', lw=2.0)
plt.semilogy(mrng, strasserwid, '-', color='orange', lw=2.0)
h1 = plt.semilogy(mrng, odrwid, 'r-', lw=2.0)
h0 = plt.semilogy(mrng, blwid, 'r--', lw=2.0)
'''
plt.legend((h0[0], h1[0]), ['log L = a + b Mw','log W = a (Mw - 9.7)**2 + b','Present Study'],loc=4)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small')
'''
cb = plt.colorbar()
cb.set_label('Dip at Epicentre')
plt.ylabel('Fault Width (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.title('Interface Width')
plt.ylim([10, 500])
plt.xlim([7.0, 9.6])

if usehistoric == True:
    plt.savefig('interface_lw_dip_hist.pdf',dpi=300,format='pdf', bbox_inches='tight')
else:
    plt.savefig('interface_lw_dip.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()
"""
'''
########################################################################################
'''
"""
# plot by region
fig = plt.figure(4, figsize=(16, 6.5))
for i, reg in enumerate(reg_sreg):
    if reg == '':
        reg_sreg[i] = 'None'
        
regions = list(unique(reg_sreg))
ncolours = len(regions)
cmap = plt.cm.get_cmap('Spectral', ncolours)

'''
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/Paired_10.cpt'
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/cosam12.cpt'
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
ncolours = len(regions)# + 3
cmap = cpt2colormap(cptfile, ncolours)[0]
'''
cs = (cmap(arange(ncolours)))

#idx = [0, 3, 4, 5, 6, 7, 8, 9, 11]
#cs = cs[idx]

for j, reg in enumerate(regions):
    tmplen = []
    tmpwid = []
    tmpdip = []
    tmpmag = []
    for i, evreg in enumerate(reg_sreg):
        if evreg == reg:
            tmplen.append(reg_flen[i])
            tmpwid.append(reg_fwid[i])
            tmpdip.append(reg_sdip[i])
            tmpmag.append(reg_fmag[i])

    ax = plt.subplot(1, 2, 1)
    plt.semilogy(tmpmag, tmplen, 'o', ms=10, color=[cs[j][0],cs[j][1],cs[j][2]])

    ax = plt.subplot(1, 2, 2)
    plt.semilogy(tmpmag, tmpwid, 'o', ms=10, color=[cs[j][0],cs[j][1],cs[j][2]])

regions[0] = 'None' 

ax = plt.subplot(1, 2, 1)
plt.semilogy(mrng, odrlen, 'k-', lw=2.0)
plt.ylabel('Fault Length (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([30, 2000])
plt.xlim([7.0, 9.6])
plt.legend(regions,loc='lower right',numpoints=1)
#set_legend(plt, regions, 'lower right')

ax = plt.subplot(1, 2, 2)
#plt.semilogy(mrng, odrwid, 'k-', lw=2.0)
plt.semilogy(mrng, blwid, 'k-', lw=2.0)
plt.ylabel('Fault Width (km)')
plt.xlabel('Magnitude')
plt.ylim([10, 300])
plt.xlim([7.0, 9.6])
#set_legend(plt, regions, 'lower right')

if usehistoric == True:
    plt.savefig('interface_lw_region_hist.pdf',dpi=300,format='pdf', bbox_inches='tight')
else:
    plt.savefig('interface_lw_region.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()

#output = open('shakedat_merge_dip.pkl', 'wb')

# this only used to get centroid!!!
'''
# Pickle dictionary using protocol 0.
pickle.dump(shakedat, output, -1)
output.close()
'''
"""