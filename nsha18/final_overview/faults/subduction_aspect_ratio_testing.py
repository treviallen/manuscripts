# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:53:36 2018

@author: u56903
"""
# get mag 2 len in km**2
def mag2len_AH17(mw, rtype):
    b = 0.63
    if rtype == 'inter':
        a = -2.90
        sig = 0.182
    elif rtype == 'intra':
        a = -3.03
        sig = 0.14
    elif rtype == 'outer':
        a = -2.87
        sig = 0.08
    elif rtype == 'ss':
        a = -2.81
        sig = 0.15
        
    return 10**(a + b * mw), sig
    
# get mag 2 area in km**2
def mag2area_AH17_other(mw, rtype):
    b = 0.96
    if rtype == 'intra':
        a = -3.89
        sig = 0.19
    elif rtype == 'outer':
        a = -3.89
        sig = 0.11
    elif rtype == 'ss':
        a = -4.04
        sig = 0.20
        
    return 10**(a + b * mw), sig

# get mag 2 len in km**2
def len2wid_AH17(length, rtype):
    from numpy import log10
    b = 0.74
    if rtype == 'inter':
        a = 0.39
        sig = 0.15
    elif rtype == 'intra':
        a = 0.35
        sig = 0.13
    elif rtype == 'outer':
        a = 0.04
        sig = 0.09
    elif rtype == 'ss':
        a = -0.22
        sig = 0.18
        
    return 10**(a + b * log10(length)), sig

#from fault_tools import mag2len_AH17, len2wid_AH17
import matplotlib.pyplot as plt
from numpy import arange, ones_like
import matplotlib as mpl
mpl.style.use('classic')

mrng = arange(6.5, 8.6, 0.1)
ar_ones = ones_like(mrng)

# get lens from mw
inter_lens = mag2len_AH17(mrng, 'inter')[0]
intra_lens = mag2len_AH17(mrng, 'intra')[0]

# get wids from lens
inter_wids = len2wid_AH17(inter_lens, 'inter')[0]
intra_wids = len2wid_AH17(intra_lens, 'intra')[0]

# get aspect ratio
inter_ar = inter_lens / inter_wids
intra_ar = intra_lens / intra_wids

fig = plt.figure(1, figsize=(6,6))

plt.plot(mrng, inter_ar, '-', c='darkorange', lw=2., label='Interface')
plt.plot(mrng, intra_ar, '-', c='dodgerblue', lw=2., label='Intraslab')
plt.plot(mrng, ar_ones, 'k--', label='1:1')
plt.grid()
plt.legend(loc=2)
#plt.title('Subduction Aspect Ratio')
plt.xlabel('Magnitude (MW)', fontsize=16)
plt.ylabel('Aspect Ratio', fontsize=16)
plt.savefig('Subduction_Aspect_Ratio.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()