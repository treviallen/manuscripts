# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 08:39:02 2015

@author: tallen
"""
# does not do iteration - just for plotting purposes
def calc_risk_integral(RTGM, beta, SAs, Probs):
    from scipy.stats import norm, lognorm
    from numpy import array, arange, exp, log, trapz, interp, isinf, where
    from scipy import interpolate
    from misc_tools import extrap1d
    
    FRAGILITY_AT_RTGM = 0.10 
    BETA = 0.6
    AFE4UHGM = - log( 1 - 0.02 )/ 50 # exceedance frequency for 1/2475 yrs
    TARGET_RISK = - log( 1 - 0.01 ) / 50
    
    '''
    SAs = array([ 0.1613, 0.1979, 0.2336, 0.3385, 0.4577, 0.5954, 0.7418, 0.7905, 0.9669, 1.1697])
    Probs = array([0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])
    '''
    # get uniform hazard at 1/2475
    idx = where(isinf(Probs) == False)[0]
    Probs = Probs[idx]
    SAs = SAs[idx]
    
    UHGM = exp((interp(log(AFE4UHGM), log(Probs[::-1]), log(SAs[::-1]))))
    
    # up sample hazard curve
    UPSAMPLING_FACTOR = 1.05
    SMALLEST_SA = min([min(SAs), UHGM/20])
    LARGEST_SA  = max([max(SAs), UHGM*20])
    
    upSAs = exp(arange(log(SMALLEST_SA),log(LARGEST_SA),log(UPSAMPLING_FACTOR)))
    
    f_i = interpolate.interp1d(log(SAs), log(Probs))
    f_x = extrap1d(f_i)
    upProbs = exp(f_x(log(upSAs)))
    '''
    upSAs = SAs
    upProbs = Probs
    '''
    # get fragility curve
    FragilityCurve = {}
    FragilityCurve['Median'] = RTGM / exp( norm.ppf( FRAGILITY_AT_RTGM ) * BETA )  
    FragilityCurve['PDF'] = lognorm.pdf(upSAs,BETA,scale=(FragilityCurve['Median']))
    FragilityCurve['CDF'] = lognorm.cdf(upSAs,BETA,scale=(FragilityCurve['Median']))
    FragilityCurve['SAs'] = upSAs
    FragilityCurve['Beta'] = BETA 
    
    # do risk integral
    Integrand = FragilityCurve['PDF'] * upProbs
    Risk = trapz(Integrand, upSAs)
    
    # calculate collapse probability
    CollapseProb = 1 - exp(-50 * Risk)
    
    RiskCoefficient = RTGM / UHGM
    
    return upProbs, upSAs, FragilityCurve, Integrand, CollapseProb

from numpy import array, arange, sin, cos, radians, log10, interp, log, exp
import matplotlib.pyplot as plt
import matplotlib as mpl  
from misc_tools import plttext
from os import path
import matplotlib as mpl
from misc_tools import dictlist2array
from hazard_tools import get_nsha18_city_curve
mpl.style.use('classic')


mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
cmap = plt.cm.get_cmap('bwr', 2)
cs = (cmap(arange(2)))

###############################################################################
# parse first job file to define plotting order
###############################################################################

pltT = 'PGA'

hazpath = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles'
hazcurvefile = path.join(hazpath, ''.join(('hazard_curve-mean-',pltT,'_1.csv')))

curveDict = get_nsha18_city_curve(['Darwin', 'Newcastle'], hazcurvefile)

###################################################################################
# set up figure
###################################################################################

figure = plt.figure(1,figsize=(16,12))

colours = ['r', 'b']

locs = ['Darwin, NT', 'Newcastle, NSW']
#RTGMs = [0.061755, 0.063473] # at sites above
RTGMs = [0.061248, 0.060562] # at sites above
beta = 0.6
targetProb = 1/2475.

###################################################################################


'''
Plot for hazard value
'''
ax = plt.subplot(3, 2, 1)

tstr = ''.join(str(pltT).split('.'))

    
'''
Plot hazard curve
'''
for col, siteCurve in zip(colours, curveDict):
    # interp haz at 1/2475-year
    targetSA   = exp(interp(log(targetProb), log(siteCurve['poe_probs_annual'][::-1]), log(siteCurve['imls'][::-1])))
    print targetSA
    
    upProbs, upSAs, FragilityCurve, Integrand, CollapseProb \
        = calc_risk_integral(targetSA, beta, siteCurve['imls'], siteCurve['poe_probs_annual'])
        
    plt.loglog(upSAs, upProbs, '-', lw=3., color=col)
            
# plt 2% in 50 yr
plt.loglog([0.3E-2, 0.3], [targetProb, targetProb], 'k--', lw = 1.5)

# plot vert RTGM
k = 0
for col, siteCurve in zip(colours, curveDict):
    targetSA  = exp(interp(log(targetProb), log(siteCurve['poe_probs_annual'][::-1]), log(siteCurve['imls'][::-1])))

    # get text xy
    xtxt = 10**(-0.3 * cos(radians(k*180 + 45)) + log10(targetSA) + k*0.3)
    ytxt = 10**(-0.3 * sin(radians(k*180 + 45)) + log10(targetProb) + k*0.3)
    
    if k == 0:
        xtxt = .014
        ytxt = 10**(log10(targetProb) - 0.55)
    else:
        xtxt = 0.08
        ytxt = 10**(log10(targetProb) + 0.5)
    	
    plt.loglog([targetSA, targetSA], [1E-4, 0.1], '--', lw=1.5, color=col)
    ax.annotate('MCE GM = '+ str("%0.3f" % targetSA)+' g', xy=([targetSA, targetProb]), xytext=(xtxt, ytxt), color=col, \
                arrowprops=dict(fc=col, shrink=0.05,width=2, ec=col), fontsize=12)
    k += 1

i = 0
if i == 0:
    #plt.legend(locs ,loc=0, fontsize=12)
    plt.ylabel('P[Sa > a] (per yr)', fontsize=14)
plt.ylim([1E-4, 0.02])
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(pltT)+' s) in g', fontsize=14)

if pltT == 'PGA':
    plt.title('NSHA18 PGA Hazard Curves', fontsize=14)
else:
    plt.title('NSHA18 Hazard Curves [Sa('+pltT[3:6]+' s)]', fontsize=14)

#plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
#plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
plt.grid(which='both')

###################################################################################

'''
Plot both fragility curve and risk integral
'''
k = 0
for col, siteCurve in zip(colours, curveDict):

    targetSA = exp(interp(log(targetProb), log(siteCurve['poe_probs_annual'][::-1]), log(siteCurve['imls'][::-1])))
    	
    upProbs, upSAs, FragilityCurve, Integrand, CollapseProb \
        = calc_risk_integral(targetSA, beta, siteCurve['imls'], siteCurve['poe_probs_annual'])
    
    # plot fragility curve
    ax = plt.subplot(3, 2, 3)
    plt.semilogx(upSAs, FragilityCurve['PDF'], '-', lw=3., color=col)
    plt.semilogx([targetSA, targetSA], [0, 7.], '--', lw=1.5, color=col)
    
    # plot risk integrand
    ax = plt.subplot(3, 2, 5)
    plt.semilogx(upSAs, Integrand, '-', lw=3., color=col)
    plt.semilogx([targetSA, targetSA], [0, 5E-3], '--', lw=1.5, color=col)
    
    # write prob collapse
    ptext = 'P[Collapse] = '+str("%0.2f" % (CollapseProb*100))+'% in 50 yrs'
    print ptext
    ylim = ax.get_ylim()
    plt.annotate(ptext, (0,0), fontsize=12, xytext=(0.02, .9 - k*.1), xycoords='axes fraction', color=col)
    k += 1    

###################################################################################

# make pretty                                            
ax = plt.subplot(3, 2, 1)
if i == 0:
    plt.legend(locs ,loc=1, fontsize=13)
    plt.ylabel('$P[Sa > a]$ (per yr)', fontsize=16)
plt.ylim([1E-4, 0.1])
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)

if pltT == 'PGA':
    plt.title('NSHA18 PGA Hazard Curves', fontsize=14)
else:
    plt.title('NSHA18 Hazard Curves [Sa('+pltT[3:6]+' s)]', fontsize=14)

plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

###################################################################################

ax = plt.subplot(3, 2, 3)
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
plt.title('Fragility Curves', fontsize=14)
plt.ylabel('$P_{f}(a)$', fontsize=16)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

###################################################################################

ax = plt.subplot(3, 2, 5)
plt.xlim([0.3E-2, 0.3])
plt.ylim([0, 3E-3])
plt.title('Risk Integrand', fontsize=14)
#plt.xlabel('Sa('+pltT[3:6]+' s) in g', fontsize=14)
plt.xlabel('PGA (g)', fontsize=14)
plt.ylabel('$P[Sa > a] . P_{f}(a)$', fontsize=16)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
# scale y-axis
scalefact = 10.**-5
ax.yaxis.get_major_formatter().set_powerlimits((-2, 3))

#########################################################################################
'''
Plot for RTGM 
'''
ax = plt.subplot(3, 2, 2)
    
'''
Plot hazard curve
'''
for col, siteCurve, RTGM in zip(colours, curveDict, RTGMs):
    upProbs, upSAs, FragilityCurve, Integrand, CollapseProb \
        = calc_risk_integral(RTGM, beta, siteCurve['imls'], siteCurve['poe_probs_annual'])
    plt.loglog(upSAs, upProbs, '-', lw=3., color=col)
            
# plt 2% in 50 yr
plt.loglog([0.3E-2, 0.3], [targetProb, targetProb], 'k--', lw = 1.5)

# plot vert RTGM
k = 0
for col, siteCurve, RTGM in zip(colours, curveDict, RTGMs):
    plt.loglog([RTGM, RTGM], [1E-4, 0.1], '--', lw=1.5, color=col)
    
    if k == 0:
        xtxt = .014
        ytxt = 10**(log10(targetProb) - 0.55)
    else:
        xtxt = 0.08
        ytxt = 10**(log10(targetProb) + 0.5)
        
    ax.annotate('RTGM = '+ str("%0.3f" % RTGM)+'g', xy=([RTGM, targetProb]), xytext=(xtxt, ytxt), color=col, \
                arrowprops=dict(fc=col, shrink=0.05,width=2, ec=col), fontsize=12)
    
    k += 1
    
plt.ylabel('P[Sa > a] (per yr)', fontsize=16)
plt.ylim([1E-4, 0.02])
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(pltT)+' s) in g', fontsize=14)
if pltT == 'PGA':
    plt.title('NSHA18 PGA Hazard Curves', fontsize=14)
else:
    plt.title('NSHA18 Hazard Curves [Sa('+pltT[3:6]+' s)]', fontsize=14)
    
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

###################################################################################

'''
Plot both fragility curve and risk integral
'''
k = 0
for col, siteCurve, RTGM in zip(colours, curveDict, RTGMs):
    upProbs, upSAs, FragilityCurve, Integrand, CollapseProb \
        = calc_risk_integral(RTGM, beta, siteCurve['imls'], siteCurve['poe_probs_annual'])
        
    # plot fragility curve
    ax = plt.subplot(3, 2, 4)
    plt.semilogx(upSAs, FragilityCurve['PDF'], '-', lw=3., color=col)
    plt.semilogx([RTGM, RTGM], [0, 7], '--', lw=1.5, color=col)
    
    # plot risk integrand
    ax = plt.subplot(3, 2, 6)
    plt.semilogx(upSAs, Integrand, '-', lw=3., color=col)
    plt.semilogx([RTGM, RTGM],[0, 5E-3], '--', lw=1.5, color=col)
    
    # write prob collapse
    ptext = 'P[Collapse] = '+str("%0.2f" % (CollapseProb*100))+'% in 50 yrs'
    print ptext
    ylim = ax.get_ylim()
    plt.annotate(ptext, (0,0), fontsize=12, xytext=(0.02, .9 - k*.1), xycoords='axes fraction', color=col)
    k += 1
                
###################################################################################

# make pretty                                            
ax = plt.subplot(3, 2, 2)
if i == 0:
    #plt.legend(locs ,loc=0, fontsize=12)
    plt.ylabel('$P[Sa > a]$ (per yr)', fontsize=16)
plt.ylim([1E-4, 0.1])
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)

if pltT == 'PGA':
    plt.title('NSHA18 PGA Hazard Curves', fontsize=14)
else:
    plt.title('NSHA18 Hazard Curves [Sa('+pltT[3:6]+' s)]', fontsize=14)
    
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)

###################################################################################

ax = plt.subplot(3, 2, 4)
plt.xlim([0.3E-2, 0.3])
#plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
plt.ylabel('$P_{f}(a)$', fontsize=16)
plt.title('Fragility Curves', fontsize=14)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
# annotate

###################################################################################

ax = plt.subplot(3, 2, 6)
plt.xlim([0.3E-2, 0.3])
plt.ylim([0, 3E-3])
plt.title('Risk Integrand', fontsize=14)
#plt.xlabel('Sa('+pltT[3:6]+' s) in g', fontsize=14)
plt.xlabel('PGA (g)', fontsize=14)
plt.ylabel('$P[Sa > a] . P_{f}(a)$', fontsize=16)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
# scale y-axis
scalefact = 10.**-5
ax.yaxis.get_major_formatter().set_powerlimits((-2, 3))


plt.savefig('darwin_sydney_rtgm.png', bbox_inches='tight', format='png', dpi=300)
plt.show()
