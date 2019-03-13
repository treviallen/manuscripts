# use historical interface data
usehistoric = True
plt_inter_only = False

import pickle
from shakemap_tools import *
from mapping_tools import distance, reckon
from numpy import array, sqrt, nan, isnan, arange, abs, unique, hstack, savetxt, corrcoef, vstack, \
                  logical_and, mean, median, std, log10, ones, logspace, exp, max, signbit, log, exp, linspace
#from make_slab_fault import make_slab_matrix
from os import path, sep
from fault_tools import *
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from scipy.stats import linregress
from mag_tools import mw2m0
from misc_tools import plttext
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 
params = {'mathtext.default': 'regular' }

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
    if ttyp[i] == 'i':
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

'''
################################################################################
'''        
# read gradients from interface models
# read len model coeffs
lines = open('len1_coeffs.txt', 'rb').readlines()
lengrd1 = float(lines[1].strip())
lenc = float(lines[0].strip())

# read wid model coeffs
lines = open('wid1_coeffs.txt', 'rb').readlines()
widgrd1 = float(lines[1].strip())
widc = float(lines[0].strip())

# read area model coeffs
lines = open('area1_coeffs.txt', 'rb').readlines()
areagrd1 = float(lines[1].strip())
areac = float(lines[0].strip())

'''
# read bi-linear area coefficients
lines = open('len_coeffs.txt', 'rb').readlines()
lengrd1 = float(lines[0].strip())
lenc = float(lines[1].strip())
lengrd2 = float(lines[2].strip())
lenhmag = float(lines[3].strip())

# read bi-linear width coefficients
lines = open('bilin_wid_coeffs.txt', 'rb').readlines()
widgrd = float(lines[0].strip())
widc = float(lines[1].strip())
widhmag = float(lines[2].strip())

# read bi-linear area coefficients
lines = open('area_coeffs.txt', 'rb').readlines()
areagrd1 = float(lines[0].strip())
areac = float(lines[1].strip())
areagrd2 = float(lines[2].strip())
areahmag = float(lines[3].strip())
'''

# read mxs model coeffs
lines = open('mxs_coeffs.txt', 'rb').readlines()
mxsgrd = float(lines[1].strip())
mxsc = float(lines[0].strip())

# read mxs model coeffs
lines = open('avs_coeffs.txt', 'rb').readlines()
avsgrd = float(lines[1].strip())
avsc = float(lines[0].strip())

# read bi-linear width coefficients
lines = open('LW_coeffs.txt', 'rb').readlines()
lwgrd = float(lines[0].strip())
lwc = float(lines[1].strip())
lwhmag = float(lines[2].strip())

'''
################################################################################
'''
# def funtion
def linear_reg_fix_slope(c, x):
    from numpy import log10, zeros_like
    grd = grad # set gradient
    ans = c[0] + x * grd
    return ans
'''
################################################################################
'''
# set data
ftypes = ['Intraslab', 'Outer Rise', 'Strike-Slip']
grads = [lengrd1, widgrd1, areagrd1, mxsgrd, avsgrd, lwgrd]
fmag = [intra_fmag, outer_fmag, trans_fmag]
flen = [intra_flen, outer_flen, trans_flen]
fwid = [intra_fwid, outer_fwid, trans_fwid]
farea = [intra_farea, outer_farea, trans_farea]
fmxs = [intra_fmxs, outer_fmxs, trans_fmxs]
favs = [intra_favs, outer_favs, trans_favs]
mrng = arange(7.2, 9.61, 0.01)
omrng = arange(7.0, 8.71, 0.01)

from gmt_tools import cpt2colormap 
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/GMT_no_green.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cmap = cpt2colormap(cptfile, 6)[0]

#cmap = plt.cm.get_cmap('hsv', 5)
#cs = (cmap(arange(5)))
#cols = [cs[1], cs[2], cs[3]]

cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temp-c.cpt'
ncolours = 6
cmap2, zvals = cpt2colormap(cptfile, ncolours, rev=False)
cs2 = (cmap2(arange(ncolours)))
#cs2 = vstack(([cs2[0], cs2[2:4], cs2[5:]]))
#cols = [cs2[1], cs2[2], cs2[3]]

ncolours = 6
cmap = plt.cm.get_cmap('Spectral', ncolours)
colours = (cmap(arange(ncolours)))
cs2 = vstack(([cs2[1:], cs2[0]]))
#cs2 = vstack(([cs2[0:3], cs2[4:]]))
cols = [cs2[1], cs2[2], cs2[3]]

#cs = cs2
colours = ['b', 'cyan', 'gold', 'darkorange','r', 'indigo']
cs2 = colours
'''
################################################################################
'''
# plot length data
import matplotlib.pyplot as plt

# plot length
fig = plt.figure(1, figsize=(11.5,17))
ax = plt.subplot(3, 2, 1, aspect=1)
plt.axis('equal')
'''
s2 = plt.semilogy(intra_fmag, intra_flen, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)
s3 = plt.semilogy(outer_fmag, outer_flen, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)
s4 = plt.semilogy(trans_fmag, trans_flen, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)
s5 = plt.semilogy(hist_fmag, hist_flen,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
s1 = plt.semilogy(inter_fmag, inter_flen, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)
'''
s2 = plt.semilogy(intra_fmag, intra_flen, 'o', color=colours[1], ms=7)
s3 = plt.semilogy(outer_fmag, outer_flen, 'o', color=colours[2], ms=7)
s4 = plt.semilogy(trans_fmag, trans_flen, 'o', color=colours[3], ms=7)
s5 = plt.semilogy(hist_fmag, hist_flen,   '^', color=colours[4], ms=9)
s1 = plt.semilogy(inter_fmag, inter_flen, 'o', color=colours[0], ms=7)



# plot interface
cinter_len = 10**(lenc + lengrd1*mrng)
'''

# plot bi-linear interface
cinter_len = 10**(lenc + lengrd1*mrng)
idx = mrng > lenhmag
ylen = lenc + lengrd1 * lenhmag
cinter_len[idx] = 10**(lengrd2 * (mrng[idx]-lenhmag) + ylen)
'''
h = plt.semilogy(mrng, cinter_len, 'k-', lw=2.0)
hh = [h[0]]

# fit lengths
savec = []
tabtxt = 'Length\nType,a,SEa,b,sig,Range\n'

for i in range(0, len(fmag)):
    data = odrpack.RealData(array(fmag[i]),log10(array(flen[i])))
    grad = grads[0]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    savec.append(c)
    regmrng = arange(min(fmag[i]), max(fmag[i])+0.01, 0.01)
    clen = 10**(c + grad * regmrng)
    print i
    h = plt.semilogy(regmrng, clen, '-', color=cs2[i+1], lw=2.0)
    hh.append(h[0])
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
savetxt('other_len_coeffs.txt', savec, delimiter='\t') 

# plot other models
wc94 = mag2ruplen_WC94(omrng, 'ss')
h = plt.semilogy(omrng, wc94, '--', lw=2.0, color=cs2[4])
hh.append(h[0])

st10 = mag2srl_St10intra(omrng)[0]
h = plt.semilogy(omrng, st10, '--', lw=2.0, color=cs2[5])
hh.append(h[0])


plt.ylabel('Rupture Length (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([20, 2000])
plt.xlim([7.0, 9.])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.05, ymax, '(A)', fontsize=18)

'''
################################################################################
'''
# plot width data

ax = plt.subplot(3, 2, 2, aspect=1)
plt.axis('equal')
'''
plt.semilogy(intra_fmag, intra_fwid, 'o', color=[cs[1][0],cs[1][1],cs[1][2]], ms=7)
plt.semilogy(outer_fmag, outer_fwid, 'o', color=[cs[2][0],cs[2][1],cs[2][2]], ms=7)
plt.semilogy(trans_fmag, trans_fwid, 'o', color=[cs[3][0],cs[3][1],cs[3][2]], ms=7)
plt.semilogy(hist_fmag, hist_fwid,   '^', color=[cs[4][0],cs[4][1],cs[4][2]], ms=9)
plt.semilogy(inter_fmag, inter_fwid, 'o', color=[cs[0][0],cs[0][1],cs[0][2]], ms=7)
'''
plt.semilogy(intra_fmag, intra_fwid, 'o', color=colours[1], ms=7)
plt.semilogy(outer_fmag, outer_fwid, 'o', color=colours[2], ms=7)
plt.semilogy(trans_fmag, trans_fwid, 'o', color=colours[3], ms=7)
plt.semilogy(hist_fmag, hist_fwid,   '^', color=colours[4], ms=9)
plt.semilogy(inter_fmag, inter_fwid, 'o', color=colours[0], ms=7)

# plot interface
#mrng = arange(7.2, 8.51, 0.01)
cinter_wid = 10**(widc + widgrd1*mrng)
'''
idx = mrng > widhmag
ywid = widc + widgrd*widhmag
cinter_wid[idx] = 10**(widc + widgrd*widhmag)
'''
plt.semilogy(mrng, cinter_wid, 'k-', lw=2.0)

tabtxt += 'Width\nType,a,SEa,b,sig,Range\n'

# fit widgths
savec = []
for i in range(0, len(fmag)):
    data = odrpack.RealData(array(fmag[i]),log10(array(fwid[i])))
    grad = grads[1]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    savec.append(c)
    regmrng = arange(min(fmag[i]), max(fmag[i])+0.01, 0.01)
    cwid = 10**(c + grad * regmrng)
    
    plt.semilogy(regmrng, cwid, '-', color=cs2[i+1], lw=2.0)
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
savetxt('other_wid_coeffs.txt', savec, delimiter='\t') 

# plot other models
st10 = mag2wid_St10intra(omrng)[0]
h = plt.semilogy(omrng, st10, '--', lw=2.0, color=cs2[5])
wc94 = mag2rupwid_WC94(omrng, 'ss')
h = plt.semilogy(omrng, wc94, '--', lw=2.0, color=cs2[4])

plt.ylabel('Rupture Width (km)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([10, 300])
plt.xlim([7.0, 9.])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.05, ymax, '(B)', fontsize=18)

'''
################################################################################
'''
# plot area data

ax = plt.subplot(3, 2, 3, aspect=1)
plt.axis('equal')

plt.semilogy(intra_fmag, intra_farea, 'o', color=colours[1], ms=7)
plt.semilogy(outer_fmag, outer_farea, 'o', color=colours[2], ms=7)
plt.semilogy(trans_fmag, trans_farea, 'o', color=colours[3], ms=7)
plt.semilogy(hist_fmag, hist_farea,   '^', color=colours[4], ms=9)
plt.semilogy(inter_fmag, inter_farea, 'o', color=colours[0], ms=7)

# plot interface
cinter_area = 10**(areac + areagrd1*mrng)
'''
idx = mrng > areahmag
yarea = areac + areagrd1 * areahmag
cinter_area[idx] = 10**(areagrd2 * (mrng[idx]-areahmag) + yarea)
'''
plt.semilogy(mrng, cinter_area, 'k-', lw=2.0)

tabtxt += 'Area\nType,a,SEa,b,sig,Range\n'

# fit areagths
savec = []
for i in range(0, len(fmag)):
    data = odrpack.RealData(array(fmag[i]),log10(array(farea[i])))
    grad = grads[2]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    savec.append(c)
    regmrng = arange(min(fmag[i]), max(fmag[i])+0.01, 0.01)
    carea = 10**(c + grad * regmrng)
    
    plt.semilogy(regmrng, carea, '-', color=cs2[i+1], lw=2.0)
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
savetxt('other_area_coeffs.txt', savec, delimiter='\t') 

# plot other models
#i06 = mag2area_I06intra(omrng)[0]
#h = plt.semilogy(omrng, i06, 'm--', lw=2.0)#, color=cs2[5])
st10 = mag2area_St10intra(omrng)[0]
h = plt.semilogy(omrng, st10, '--', lw=2.0, color=cs2[5])
wc94 = mag2area_WC94(omrng, 'ss')
h = plt.semilogy(omrng, wc94, '--', lw=2.0, color=cs2[4])

#plt.ylabel(r'$\mathregular{Rupture Area (km^2)}$', fontsize=14)
plt.ylabel(r'Rupture Area $\mathregular{(km^{2}}$)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.legend((s1[0], s2[0], s3[0], s4[0], s5[0]), ['Interface', \
           'Intraslab','Outer-Rise','Strike-Slip','Other Interface'],loc=4,numpoints=1,fontsize=12)
plt.ylim([500,2E5])
plt.xlim([7.0, 9.])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.05, ymax, '(C)', fontsize=18)

'''
################################################################################
'''
# plot max slip

ax = plt.subplot(3, 2, 4, aspect=1)
plt.axis('equal')

plt.semilogy(intra_fmag, intra_fmxs, 'o', color=colours[1], ms=7)
plt.semilogy(outer_fmag, outer_fmxs, 'o', color=colours[2], ms=7)
plt.semilogy(trans_fmag, trans_fmxs, 'o', color=colours[3], ms=7)
plt.semilogy(hist_fmag, hist_fmxs,   '^', color=colours[4], ms=9)
plt.semilogy(inter_fmag, inter_fmxs, 'o', color=colours[0], ms=7)

# plot interface
mrng = arange(7.2, 9.01, 0.01)
cinter_mxs = 10**(mxsc + mxsgrd*mrng)
plt.semilogy(mrng, cinter_mxs, 'k-', lw=2.0)

tabtxt += 'Max Slip\nType,a,SEa,b,sig,Range\n'

# fit mxsgths
savec = []
for i in range(0, len(fmag)):
    data = odrpack.RealData(array(fmag[i]),log10(array(fmxs[i])))
    grad = grads[3]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    savec.append(c)
    regmrng = arange(min(fmag[i]), max(fmag[i])+0.01, 0.01)
    cmxs = 10**(c + grad * regmrng)
    
    plt.semilogy(regmrng, cmxs, '-', color=cs2[i+1], lw=2.0)
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
savetxt('other_mxs_coeffs.txt', savec, delimiter='\t') 

wc94 = mag2maxs_WC94(omrng, 'ss')[0]
h = plt.semilogy(omrng, wc94, '--', lw=2.0, color=cs2[4])

plt.ylabel('Maximum Slip (m)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([.5,100])
plt.xlim([7.0, 9.])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.05, ymax, '(D)', fontsize=18)

plt.legend(hh, ['Interface', 'Intraslab','Outer-Rise','Strike-Slip','WC94 SS','S10 Intraslab'], \
           loc=4,numpoints=1,fontsize=12)

'''
################################################################################
'''
# plot average slip

ax = plt.subplot(3, 2, 5, aspect=1)
plt.axis('equal')

plt.semilogy(intra_fmag, intra_favs, 'o', color=colours[1], ms=7)
plt.semilogy(outer_fmag, outer_favs, 'o', color=colours[2], ms=7)
plt.semilogy(trans_fmag, trans_favs, 'o', color=colours[3], ms=7)
plt.semilogy(hist_fmag, hist_favs,   '^', color=colours[4], ms=9)
plt.semilogy(inter_fmag, inter_favs, 'o', color=colours[0], ms=7)

# plot interface
cinter_avs = 10**(avsc + avsgrd*mrng)
plt.semilogy(mrng, cinter_avs, 'k-', lw=2.0)

tabtxt += 'Av Slip\nType,a,SEa,b,sig,Range\n'

# fit avsgths
savec = []
for i in range(0, len(fmag)):
    data = odrpack.RealData(array(fmag[i]),log10(array(favs[i])))
    grad = grads[4]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    savec.append(c)
    regmrng = arange(min(fmag[i]), max(fmag[i])+0.01, 0.01)
    cavs = 10**(c + grad * regmrng)
    
    plt.semilogy(regmrng, cavs, '-', color=cs2[i+1], lw=2.0)
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
savetxt('other_avs_coeffs.txt', savec, delimiter='\t') 

wc94 = mag2avs_WC94(omrng, 'ss')[0]
h = plt.semilogy(omrng, wc94, '--', lw=2.0, color=cs2[4])

plt.ylabel('Average Slip (m)', fontsize=14)
plt.xlabel('Magnitude', fontsize=14)
plt.ylim([.1,50])
plt.xlim([7.0, 9.])
ylims = log10(ax.get_ylim())
ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])
plt.text(7.05, ymax, '(E)', fontsize=18)

'''
################################################################################
'''
# plot L-W data

ax = plt.subplot(3, 2, 6, aspect=1)
plt.axis('equal')

plt.semilogy(intra_flen, intra_fwid, 'o', color=colours[1], ms=7)
plt.semilogy(outer_flen, outer_fwid, 'o', color=colours[2], ms=7)
plt.semilogy(trans_flen, trans_fwid, 'o', color=colours[3], ms=7)
plt.semilogy(hist_flen, hist_fwid,   '^', color=colours[4], ms=9)
plt.semilogy(inter_flen, inter_fwid, 'o', color=colours[0], ms=7)

# plot interface
#mrng = arange(7.2, 8.51, 0.01)
lenplt = linspace(log10(20),log10(1500),300)
cinter_lw = 10**(lwc + lwgrd*lenplt)
ywid = lwhmag
idx = cinter_lw > 10**ywid
cinter_lw[idx] = 10**(ywid)
plt.semilogy(10**lenplt, cinter_lw, 'k-', lw=2.0)

tabtxt += 'LW\nType,a,SEa,b,sig,Range\n'

# fit widgths
savec = []
for i in range(0, len(flen)):
    data = odrpack.RealData(log10(array(flen[i])),log10(array(fwid[i])))
    grad = grads[5]
    
    lin = odrpack.Model(linear_reg_fix_slope)
    odr = odrpack.ODR(data, lin, beta0=[-2.])
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    if i == 0: #in-slab
        intrac = c
    savec.append(c)
    reglenplt = linspace(log10(min(flen[i])), log10(max(flen[i]))+0.01, 200)
    clw = 10**(c + grad * reglenplt)
    
    plt.semilogy(10**reglenplt, clw, '-', color=cs2[i+1], lw=2.0)
    
    tabtxt += ','.join((ftypes[i], str('%0.2f' % c), str('%0.2f' % out.sd_beta[0]), \
                        str('%0.2f' % grad), str('%0.2f' % sqrt(out.res_var)), \
                        str('%0.1f' % min(fmag[i]))+'-'+str('%0.1f' % max(fmag[i]))))+'\n'
    
oto = plt.loglog([10, 2000], [10, 2000], '--', color='0.45', lw=2, label='1:1')
plt.legend(loc=4, fontsize=13, numpoints=3)

plt.ylabel('Width (km)', fontsize=14)
plt.xlabel('Length (km)', fontsize=14)
plt.ylim([10, 200])
plt.xlim([20, 1000])
ylims = log10(ax.get_ylim())
xlims = log10(ax.get_xlim())
ymax = 10**(1.02*(ylims[1]-ylims[0]) + ylims[0])
xmin = 10**(0.025*(xlims[1]-xlims[0]) + xlims[0])
plt.text(xmin, ymax, '(F)', fontsize=18, va='bottom')

plt.savefig('reg_other_types.pdf',dpi=300,format='pdf',bbox_inches='tight')
plt.savefig('reg_other_types.png',dpi=300,format='png',bbox_inches='tight')
plt.show()

'''
################################################################################
'''
# test LW intra vs 1:1
from scipy import stats
regres = (intrac + lwgrd*log10(intra_flen)) - log10(intra_fwid)
print '\nIntra reg res fit:', median(regres), mean(regres), std(regres), corrcoef(log10(intra_flen), regres)
gradient, intercept, r_value, p_value, std_err = stats.linregress(log10(intra_flen), regres)
print "R-squared", r_value**2


otores = log10(intra_flen) - log10(intra_fwid)
print '\nIntra reg res 1:1:', median(otores), mean(otores), std(otores), corrcoef(log10(intra_flen), otores)
gradient, intercept, r_value, p_value, std_err = stats.linregress(log10(intra_flen), otores)
print "R-squared", r_value**2


'''
################################################################################
'''
# save table txt
f = open('other_coeffs.csv', 'wb')
f.write(tabtxt)
f.close()
