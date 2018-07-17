# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 08:39:02 2015

@author: tallen
"""
# does not do iteration - just for plotting purposes
def calc_risk_integral(RTGM, beta, SAs, Probs):
    from scipy.stats import norm, lognorm
    from numpy import array, arange, exp, log, trapz, interp
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
    UHGM = exp((interp(log(AFE4UHGM), log(Probs[::-1]), log(SAs[::-1]))))
    
    # up sample hazard curve
    UPSAMPLING_FACTOR = 1.05
    SMALLEST_SA = min([min(SAs), UHGM/10])
    LARGEST_SA = max([max(SAs), UHGM*10])
    
    
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

from numpy import array, arange, sin, cos, radians, log10
import matplotlib.pyplot as plt
import matplotlib as mpl  
from misc_tools import plttext

pltT = [0.2]

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
cmap = plt.cm.get_cmap('bwr', 2)
cs = (cmap(arange(2)))

figure = plt.figure(1,figsize=(16,12))

locs = ['Masset, BC', 'Saguenay, QC']
prefs = [40, 436]
RTGMs = [0.760928978, 0.72434229] # at sites above
beta = 0.6

probs = array([0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])

for i, t in enumerate(pltT):
    '''
    Plot for hazard value
    '''
    ax = plt.subplot(3, 2, 1)
    
    tstr = ''.join(str(pltT[0]).split('.'))
    hazfile = ''.join(('NBCC2015Loc_mean_hazcurves_',tstr,'.csv'))
    print hazfile

    # read risk coef
    lat = []
    lon = []
    ref = []
    haz = []
    place = []
    lines = open(hazfile).readlines()[4:]
    for line in lines:
        dat = line.strip().split(',')
        lon.append(float(dat[0]))
        lat.append(float(dat[1]))
        ref.append(float(dat[13]))
        place.append(dat[14])
        haz.append(array([float(x) for x in dat[3:13]]))
        
    '''
    Plot hazard curve
    '''
    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                print ref[j], place[j]
                col = [cs[k][0],cs[k][1],cs[k][2]]
                upProbs, upSAs, FragilityCurve, Integrand, CollapseProb = calc_risk_integral(haz[j][7], beta, haz[j], probs)
                plt.loglog(upSAs, upProbs, '-', lw=3., color=col)
                
    # plt 2% in 50 yr
    maphaz = 1/2475.
    plt.loglog([2*10**-2, 5], [maphaz, maphaz], 'k--', lw = 2)
    
    
    # plot vert RTGM
    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                col = [cs[k][0],cs[k][1],cs[k][2]]
                # get text xy
                xtxt = 10**(-0.3 * cos(radians(k*180 + 45)) + log10(haz[j][7]) + k*0.3)
                ytxt = 10**(-0.3 * sin(radians(k*180 + 45)) + log10(maphaz) + k*0.3)
                
                if k == 0:
                    xtxt = 1.2
                    ytxt = 10**(log10(maphaz) + 0.5)
                else:
                    xtxt = .14
                    ytxt = 10**(log10(maphaz) - 0.55)
                	
                plt.loglog([haz[j][7], haz[j][7]], [probs[-1], probs[0]], '--', lw=1.5, color=col)
                ax.annotate('MCE GM = '+ str("%0.2f" % haz[j][7])+' g', xy=([haz[j][7], maphaz]), xytext=(xtxt, ytxt), color=col, \
                            arrowprops=dict(fc=col, shrink=0.05,width=2, ec=col), fontsize=12)
                
    if i == 0:
        #plt.legend(locs ,loc=0, fontsize=12)
        plt.ylabel('P[Sa > a] (per yr)', fontsize=14)
    plt.ylim([probs[-1], probs[0]])
    plt.xlim([3*10**-2, 2])
    #plt.xlabel('Sa('+str(pltT[0])+' s) in g', fontsize=14)
    plt.title('2015 NBCC Hazard Curves [Sa('+str(pltT[0])+' s)]', fontsize=15)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    
    '''
    Plot both fragility curve and risk integral
    '''

    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                upProbs, upSAs, FragilityCurve, Integrand, CollapseProb = calc_risk_integral(haz[j][7], beta, haz[j], probs)
                
                col = [cs[k][0],cs[k][1],cs[k][2]]
                
                # plot fragility curve
                ax = plt.subplot(3, 2, 3)
                plt.semilogx(upSAs, FragilityCurve['PDF'], '-', lw=3., color=col)
                plt.semilogx([haz[j][7], haz[j][7]], [0, 0.5], '--', lw=1.5, color=col)
                
                # plot risk integrand
                ax = plt.subplot(3, 2, 5)
                plt.semilogx(upSAs, Integrand, '-', lw=3., color=col)
                plt.semilogx([haz[j][7], haz[j][7]], [0., 0.00035], '--', lw=1.5, color=col)
                
                # write prob collapse
                ptext = 'P[Collapse] = '+str("%0.2f" % (CollapseProb*100))+'% in 50 yrs'
                print ptext
                ylim = ax.get_ylim()
                plt.annotate(ptext, (0,0), fontsize=12, xytext=(0.02, .9 - k*.1), xycoords='axes fraction', color=col)
                    
    # make pretty                                            
    ax = plt.subplot(3, 2, 1)
    if i == 0:
        plt.legend(locs ,loc=0, fontsize=12)
        plt.ylabel('$P[Sa > a]$ (per yr)', fontsize=14)
    plt.ylim([probs[-1], probs[0]])
    plt.xlim([3*10**-2, 5])
    #plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
    plt.title('2015 NBCC Hazard Curves [Sa('+str(pltT[0])+' s)]', fontsize=13)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    
    
    ax = plt.subplot(3, 2, 3)
    plt.xlim([3*10**-2, 5])
    #plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
    plt.title('Fragility Curves', fontsize=13)
    plt.ylabel('$f_{capacity}(a)$', fontsize=14)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    
    ax = plt.subplot(3, 2, 5)
    plt.xlim([3*10**-2, 5])
    plt.ylim([0, 3.5*10**-4])
    plt.title('Risk Integrand', fontsize=13)
    plt.xlabel('Sa('+str(pltT[0])+' s) in g', fontsize=12)
    plt.ylabel('$P[Sa > a] . f_{capacity}(a)$', fontsize=14)
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
    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                print ref[j], place[j]
                col = [cs[k][0],cs[k][1],cs[k][2]]
                upProbs, upSAs, FragilityCurve, Integrand, CollapseProb = calc_risk_integral(RTGMs[k], beta, haz[j], probs)
                plt.loglog(upSAs, upProbs, '-', lw=3., color=col)
                
    # plt 2% in 50 yr
    maphaz = 1/2475.
    plt.loglog([2*10**-2, 5], [maphaz, maphaz], 'k--', lw = 2)
    
    # plot vert RTGM
    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                col = [cs[k][0],cs[k][1],cs[k][2]]
                plt.loglog([RTGMs[k], RTGMs[k]], [probs[-1], probs[0]], '--', lw=1.5, color=col)
                
                if k == 0:
                    xtxt = 1.2
                    ytxt = 10**(log10(maphaz) + 0.5)
                else:
                    xtxt = .14
                    ytxt = 10**(log10(maphaz) - 0.55)
                    
                ax.annotate('RTGM = '+ str("%0.2f" % RTGMs[k])+'g', xy=([RTGMs[k], maphaz]), xytext=(xtxt, ytxt), color=col, \
                            arrowprops=dict(fc=col, shrink=0.05,width=2, ec=col), fontsize=12)
                
    plt.ylabel('P[Sa > a] (per yr)', fontsize=14)
    plt.ylim([probs[-1], probs[0]])
    plt.xlim([3*10**-2, 2])
    #plt.xlabel('Sa('+str(pltT[0])+' s) in g', fontsize=14)
    plt.title('2015 NBCC Hazard Curves [Sa('+str(pltT[0])+' s)]', fontsize=15)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    
    '''
    Plot both fragility curve and risk integral
    '''

    for k, pr in enumerate(prefs):
        for j, r in enumerate(ref):        
            if r == pr:
                upProbs, upSAs, FragilityCurve, Integrand, CollapseProb = calc_risk_integral(RTGMs[k], beta, haz[j], probs)
                print place[j], 'CollapseProb:', CollapseProb
                
                col = [cs[k][0],cs[k][1],cs[k][2]]
                
                # plot fragility curve
                ax = plt.subplot(3, 2, 4)
                plt.semilogx(upSAs, FragilityCurve['PDF'], '-', lw=3., color=col)
                plt.semilogx([RTGMs[k], RTGMs[k]], [0, 0.6], '--', lw=1.5, color=col)
                
                # plot risk integrand
                ax = plt.subplot(3, 2, 6)
                plt.semilogx(upSAs, Integrand, '-', lw=3., color=col)
                plt.semilogx([RTGMs[k], RTGMs[k]], [0., 0.0004], '--', lw=1.5, color=col)
                
                # write prob collapse
                ptext = 'P[Collapse] = '+str("%0.2f" % (CollapseProb*100))+'% in 50 yrs'
                print ptext
                ylim = ax.get_ylim()
                plt.annotate(ptext, (0,0), fontsize=12, xytext=(0.02, .9 - k*.1), xycoords='axes fraction', color=col)
                    
    # make pretty                                            
    ax = plt.subplot(3, 2, 2)
    if i == 0:
        #plt.legend(locs ,loc=0, fontsize=12)
        plt.ylabel('$P[Sa > a]$ (per yr)', fontsize=14)
    plt.ylim([probs[-1], probs[0]])
    plt.xlim([3*10**-2, 5])
    #plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
    plt.title('2015 NBCC Hazard Curves [Sa('+str(pltT[0])+' s)]', fontsize=13)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    
    ax = plt.subplot(3, 2, 4)
    plt.xlim([3*10**-2, 5])
    #plt.xlabel('Sa('+str(t)+' s) in g', fontsize=11)
    plt.ylabel('$f_{capacity}(a)$', fontsize=14)
    plt.title('Fragility Curves', fontsize=13)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    # annotate
    
    
    ax = plt.subplot(3, 2, 6)
    plt.xlim([3*10**-2, 5])
    plt.ylim([0, 3.5*10**-4])
    plt.title('Risk Integrand', fontsize=13)
    plt.xlabel('Sa('+str(pltT[0])+' s) in g', fontsize=12)
    plt.ylabel('$P[Sa > a] . f_{capacity}(a)$', fontsize=14)
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    # scale y-axis
    scalefact = 10.**-5
    ax.yaxis.get_major_formatter().set_powerlimits((-2, 3))


plt.savefig('masset_saguenay.png', format='png', dpi=150)
plt.show()
