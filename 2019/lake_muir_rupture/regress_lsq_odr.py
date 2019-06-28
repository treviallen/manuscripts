# -*- coding: utf-8 -*-
"""
Created on Wed May 01 14:09:02 2019

@author: u56903
"""
from numpy import array, arange, log10, argsort, copyto, zeros_like
from scipy import stats
from regressions import odr_lin_reg, lsq_lin_reg, lsq_reg_with_confidence
from misc_tools import get_log_xy_locs
from fault_tools import mag2ruplen_WC94, mag2srl_WC94, mag2srl_Cea14, mag2len_L14, mag2lsr_L14
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

# parse SRL table
eqn = []
year = []
mag = []
m0  = []
dsrl = []
vsrl = []
vdsp = []

lines = open('SRL_table.csv').readlines()

for line in lines[1:]:
    dat = line.strip().split(',')
    
    eqn.append(dat[0])
    year.append(int(dat[1]))
    mag.append(float(dat[2]))
    m0.append(float(dat[3]))
    dsrl.append(float(dat[4]))
    vsrl.append(float(dat[5]))
    vdsp.append(float(dat[6]))
    
mag = array(mag)
dsrl = array(dsrl)
vsrl = array(vsrl)

# sort
idx = argsort(mag)
mag = mag[idx]
dsrl = dsrl[idx]
vsrl = vsrl[idx]

##############################################################################
# plot VSRL data
##############################################################################

fig = plt.figure(1, figsize=(13,6))
ax = plt.subplot(121)

plt.semilogy(mag, vsrl, 's', c='seagreen', label='VSRL Data')
plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Visible SRL (km)', fontsize=14)

# plt WC94
mrng = arange(4.7, 6.81, 0.001)
wc94_lrng, sig = mag2srl_WC94(mrng, 'rs')
plt.semilogy(mrng, wc94_lrng, '-', c='dodgerblue', lw=1.5, label='Wells & Coppersmith (1994; SRL)')

# plt Cea14
c14_lrng, sig = mag2srl_Cea14(mrng, 'au')
plt.semilogy(mrng, c14_lrng, '-', c='orange', lw=1.5, label='Clark et al. (2014)')

# plt L14
l14_lrng = mag2lsr_L14(mrng, 'scrrs')
plt.semilogy(mrng, l14_lrng, '-', c='r', lw=1.5, label='Leonard (2014; SRL)')

# do regression
start = [1., -5.]
m, c = lsq_lin_reg(mag, log10(vsrl), start)
lower, upper = lsq_reg_with_confidence(mag, log10(vsrl))
plt.fill_between(mag, 10**lower, 10**upper, color='0.85', edgecolor='', zorder=0)

# do odr for completeness
om, oc = odr_lin_reg(mag, log10(vsrl), start)

print('VSRL coeffs')
print('lsq', m, c)
print('odr', om, oc)

# plot new relationship
lrng = 10**(m * mrng + c)
plt.semilogy(mrng, lrng, '-', c='k', lw=1.5, label='Present Study (LSQ)')

plt.legend(loc=4, numpoints=1, fontsize=11)
plt.text(4.55, 90, '(a)', va='top', ha='left', fontsize=18)
plt.ylim([0.1, 100])

##############################################################################
# plot DSRL data
##############################################################################

ax = plt.subplot(122)

plt.semilogy(mag, dsrl, 's', c='seagreen', label='DSRL Data')
plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Detectable SRL (km)', fontsize=14)

# plt WC94
mrng = arange(4.7, 6.81, 0.01)
wc94_lrng = mag2ruplen_WC94(mrng, 'rs')
plt.semilogy(mrng, wc94_lrng, '-', c='dodgerblue', lw=1.5, label='Wells & Coppersmith (1994; RLD)')

# plt Cea14
c14_lrng, sig = mag2srl_Cea14(mrng, 'au')
plt.semilogy(mrng, c14_lrng, '-', c='orange', lw=1.5, label='Clark et al. (2014)')

# plt L14
l14_lrng = mag2len_L14(mrng, 'scrrs')
plt.semilogy(mrng, l14_lrng, '-', c='red', lw=1.5, label='Leonard (2014; FRL)')

# do regression
start = [1., -5.]
m, c = lsq_lin_reg(mag, log10(dsrl), start)
lower, upper = lsq_reg_with_confidence(mag, log10(dsrl))
plt.fill_between(mag, 10**lower, 10**upper, color='0.85', edgecolor='', zorder=0)

# do odr for completeness
om, oc = odr_lin_reg(mag, log10(dsrl), start)
#odrsl_lrng = 10**(om * mrng + oc) # keep for later

print('\nDSRL coeffs')
print('lsq', m, c)
print('odr', om, oc)

# plot new relationship
lrng = 10**(m * mrng + c)
odrsl_lrng = 10**(m * mrng + c) # keep for later
plt.semilogy(mrng, lrng, '-', c='k', lw=1.5, label='Present Study (LSQ)')

plt.legend(loc=4, numpoints=1, fontsize=11)
plt.text(4.55, 90, '(b)', va='top', ha='left', fontsize=18)

##############################################################################
# export & show
##############################################################################

plt.savefig('2019_srl_regression.png', fmt='png', dpi=300, bbox_inches='tight')
#plt.show()

##############################################################################
# get ratio of DSRL vs VSRL
##############################################################################

#fig = plt.figure(2, figsize=(19.5,6))
fig = plt.figure(2, figsize=(13,6))

rat_vdsrl = vsrl/dsrl
idx = rat_vdsrl != 1.
'''
ax = plt.subplot(131)



plt.semilogy(mag[idx], rat_vdsrl[idx], 's', c='seagreen', label='DSRL Data')

plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('VSRL/DSRL', fontsize=14)

# do odr for completeness
start = [0.5, -2.]
om, oc = odr_lin_reg(mag[idx], log10(rat_vdsrl[idx]), start)
plt_rat = 10**(om * mrng + oc)
plt.semilogy(mrng, plt_rat, '-', c='k', lw=1.5, label='Present Study (ODR)')
plt.ylim([0.1, 1])

print('\nRatio coeffs')
print('odr1', om, oc)
'''
##############################################################################

ax = plt.subplot(121)

plt.loglog(vsrl[idx], rat_vdsrl[idx], 's', c='seagreen', label='DSRL Data')

plt.grid(which='both')
plt.xlabel('VSRL (km)', fontsize=14)
plt.ylabel('VSRL/DSRL', fontsize=14)
plt.xlim([0.1, 40])

# do odr for completeness
start = [0.5, -2.]
om, oc = odr_lin_reg(log10(vsrl[idx]), log10(rat_vdsrl[idx]), start)
plt_vsrl = arange(0.2, 37.1, 0.1)
plt_rat = 10**(om * log10(plt_vsrl) + oc)
plt.semilogy(plt_vsrl, plt_rat, '-', c='k', lw=1.5, label='Present Study (ODR)')
plt.ylim([0.1, 1])

xpos = get_log_xy_locs([0.1, 40], 0.02)
ypos = get_log_xy_locs([0.1, 1], 0.98)

# set y ticks
ticks = array([0.1, 0.2, 0.5, 1.0])
ticklabels = [str(x) for x in ticks]
ax.set_yticks(ticks)
ax.set_yticklabels(ticklabels)

plt.text(xpos, ypos, '(a)', va='top', ha='left', fontsize=18)

print('\nRatio coeffs')
print('odr2', om, oc)

##############################################################################
# correct non-drsl data

ax = plt.subplot(122)

cor_dsrl = zeros_like(dsrl)
copyto(cor_dsrl, dsrl)

fidx = rat_vdsrl == 1.
cor_factor = 10**(om * log10(vsrl) + oc)

print('\ncorr factor = 10**(',om,'* log10(VSRL) +',oc,')\n' )

# fix max to 1
cidx = cor_factor > 1.
cor_factor[cidx] = 1.
cor_dsrl[fidx] = vsrl[fidx] / cor_factor[fidx]

plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel("DSRL' (km)", fontsize=14)

plt.semilogy(mag, dsrl, 's', mfc='w', label='DSRL = VSRL')
plt.semilogy(mag[idx], cor_dsrl[idx], 's', c='darkorange', label='Observed DSRL')
plt.semilogy(mag[fidx], cor_dsrl[fidx], 's', c='dodgerblue', label='Corrected DSRL')

# plt unfixed reg
plt.semilogy(mrng, odrsl_lrng, '-', c='k', lw=1.5, label="DSRL (LSQ)")

# do lin regression
m, c = lsq_lin_reg(mag, log10(cor_dsrl), start)
fdrsl_lrng = 10**(m * mrng + c)
plt.semilogy(mrng, fdrsl_lrng, '-', c='r', lw=1.5, label="DSRL' (LSQ)")

# do odr for completeness
om, oc = odr_lin_reg(mag, log10(cor_dsrl), start)
#fdrsl_lrng = 10**(om * mrng + oc)

print('\nCorrected DSRL coeffs')
print('lsq', m, c)
print('odr', om, oc)

plt.legend(loc=4, numpoints=1, fontsize=11)
plt.ylim([1, 100])

ticks = array([1.0, 10, 100])
ticklabels = [str(x).replace('.0','') for x in ticks]
ax.set_yticks(ticks)
ax.set_yticklabels(ticklabels)

plt.text(4.55, 90, '(b)', va='top', ha='left', fontsize=18)

plt.savefig('2019_vsrl-dsrl_ratio.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
















