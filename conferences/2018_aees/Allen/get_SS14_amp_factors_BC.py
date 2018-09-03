# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:50:20 2015

@author: tallen

Code generates a python dictionary of amplification factors
relative to B/C site conditions

To access amp factors for a given PGAref, site class, type:

	> C_ampfactors[PGAref][Vs30]

	e.g. 
	
	IN > C_ampfactors['0.1g']['vs1100']
	
	OUT > array([ 0.84782554,  0.75796444,  0.69797956,  0.64357284,  0.62994067,
                0.62824066,  0.6411572 ,  0.69124828,  0.85792907,  0.67049578]) 
	
Outputs CSV file of amp factors relative to B/C site conditions

"""

#from boore_atkinson_site_2008 import boore_atkinson_siteamp
from atkinson_boore_site_2006 import atkinson_boore_siteamp
from seyhan_stewart_2014 import seyhan_stewart_siteamp
from numpy import array, arange, where, log, exp, interp, hstack
from gmt_tools import cpt2colormap

# periods for NBCC2015 (PGA = 0; PGV = -1)     
T = [0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, \
     1.5, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 4., 4.5, 5., 0., -1.] # in s
#T = [0.1, 0.2, 0.3, 0.5, 1., 2., 5., 10., 0., -1.] # in s - NBCC


# Set Vs30 array
vs30 = [115, 180, 270, 412, 560, 760, 1100] # in m/s 

# frist get B/C pga4nl
target_pga4nl_BC = array([0.1, 0.2, 0.3, 0.4, 0.5]) # in g

# get B/C linear ground-motions that provide target_pga4nl_C values - AW rounding here led to minor differences in final factors
refVs30 = 760.  # reset to B/C
refT    = 0.0 # i.e. PGA
refPGA  = 0.1 # in g
pga4nl_BC = target_pga4nl_BC / seyhan_stewart_siteamp(refVs30, refT, refPGA) # B/C to C conversion factor for PGA of 1.208 based on BA08

#pga4nl_BC = target_pga4nl_C / 1.  # for testing purposes - needed to match Gail's BA08 C->B/C factors

##########################################################################
'''
Get BA08 amp factors for site class BC - this loop is not required for 
NBCC2015 site class tables, but kept for testing purposes
'''

# set empty dictionary
BC_ampfactors = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    # loop thru Vs30
    for v in vs30:
        # loop thru periods
        amp = []
        for t in T:
            if v < 1000:
                # get amp factor for each period and append to amp
                amp.append(seyhan_stewart_siteamp(v, t, pga))
                
            else:
                amp.append(atkinson_boore_siteamp(v, t, pga))
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
    
    # add velocity dictionary to pga4nl dictionary
    BC_ampfactors[str(target_pga4nl_BC[i])+'g'] = vdict
    
BC_ampfactors['periods'] = array(T)

##########################################################################
'''
Re-cast amp factors in terms of site class C 
This gets appropriate amp factors for classes C-E
'''
# set empty dictionary
C_ampfactors  = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    for v in vs30:
        # loop thru periods
        amp = []
        for t in T:
            if v < 1000:
                # get amp factor for each period relative to C and append to amp
                amp.append(seyhan_stewart_siteamp(v, t, pga) / 
                           seyhan_stewart_siteamp(760., t, pga))
                           
            else:
                amp.append(atkinson_boore_siteamp(v, t, pga))
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
    
    # add velocity dictionary to amp factor dictionary
    C_ampfactors[str(target_pga4nl_BC[i])+'g'] = vdict

# add periods to amp factor dictionary 
C_ampfactors['periods'] = array(T)

##########################################################################

##########################################################################
'''
Get CSV text for output to match NBCC tables
'''

# set params
ascii_no = range(65, 70) # ascii chars A-E
site_class = ['B', 'B/C', 'C', 'C/D', 'D', 'D/E', 'E'] # convert ascii number to char
rev_vs30 = vs30[::-1] # reverse vs30 order for tables
PGAref_header = ','+','.join(['PGAref='+str(x)+'g' for x in target_pga4nl_BC])+'\n'
NBCCtxt = ''

# loop thru period
for i, t in enumerate(T):
    # get period header
    if t == 0.0:
        T_header = ','.join(('Site Class', 'Values of F(T) for T = PGA\n'))
    elif t == -1.0:
        T_header = ','.join(('Site Class', 'Values of F(T) for T = PGV\n'))
    else:
        T_header = ','.join(('Site Class', 'Values of F(T) for T = '+str(t)+' s\n'))
    
    # add headers to NBCC text
    NBCCtxt += T_header + PGAref_header
    
    # add amp factors for each site class
    for j, sc in enumerate(site_class):
        # loop thru PGAref and add amp factor
        site_txt = ''
        for pga in target_pga4nl_BC:
            # round to 2 decimal points for output
            site_txt += ',' + str('%0.2f' % C_ampfactors[str(pga)+'g']['vs'+str(rev_vs30[j])][i])
        
        # add site class txt to NBCC text
        NBCCtxt += sc + site_txt + '\n'
    NBCCtxt += '\n\n'

# write NBCC factors to file
f = open('SS14_ampfactors_BC.csv', 'wb')
f.write(NBCCtxt)
f.close()

##########################################################################
'''
plt amp factors
'''

import matplotlib.pyplot as plt
import matplotlib as mpl

pltpga = [0.1, 0.4]

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14 

numcols = 7
cmap, zvals = cpt2colormap('/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt', numcols)
cs = (cmap(arange(numcols)))

# remove yellow
#idx=array([0, 1, 3, 4, 5])
#cs = cs[idx]

legtxt = ['Site Class B', 'Site Class B/C', 'Site Class C', 'Site Class C/D', 'Site Class D', \
          'Site Class D/E', 'Site Class E'] #, 'F(T)', 'Fa & Fv']
figure = plt.figure(1,figsize=(19,8))

# Boarcherdt 1994 for 0.1 and 0.4 g
Fa01 = [0.6, 0.7, 1.0, 1.5, 1.5]
Fa04 = [1.1, 1.0, 1.0, 0.9, 0.9]
Fv01 = [0.4, 0.6, 1.0, 2.0, 2.0]
Fv04 = [0.6, 0.7, 1.0, 1.6, 1.6]
FaT  = [0.1, 0.5]
FvT  = [0.4, 2.0]


# Finn & Whiteman for 0.1 and 0.4 g - Relative to site Sa(0.2 s)
Fa01 = array([0.7, 0.8, 1.0, 1.3, 2.1])  
Fa04 = array([1.1, 1.0, 1.0, 0.9, 0.9])  
Fv01 = array([0.4, 0.6, 1.0, 2.0, 2.0])  
Fv04 = array([0.6, 0.7, 1.0, 1.6, 1.6])  
FaT = [0.1, 0.5]
FvT = [0.4, 2.0] 


for i, pga in enumerate(pltpga):
    ax = plt.subplot(1, 2, i+1)
    plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    # loop thru vs
    for j, v in enumerate(rev_vs30):
        pltamps = C_ampfactors[str(pga)+'g']['vs'+str(v)][:-2]
        plt.semilogx(T[:-2], pltamps, '-', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])        
        	
    plt.ylim([0.5, 7.])
    plt.title('PGAref = ' + str(pga) + ' g', fontsize=18)
    plt.xlabel('Period (s)', fontsize=16)
    plt.ylabel('Amplification Factor (Relative to B/C)', fontsize=16)
    
    '''
    if i == 0:
        plt.legend(legtxt ,loc=2, fontsize=14)
    else:
        h1 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=3.)
        h2 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=3.)
        plt.legend((h1[0], h2[0]),('F(T)', 'Fa & Fv') ,loc=2, fontsize=14)
    '''
    # at John's request
    plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=4.)
    plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=4.)
    plt.legend(legtxt ,loc=2, fontsize=13)
    
    '''    
    for j, v in enumerate(rev_vs30):	
        # plt Fa, Fv
        if i == 0:
            plt.semilogx(FaT, [Fa01[j], Fa01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv01[j], Fv01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
        else:                                              
            plt.semilogx(FaT, [Fa04[j], Fa04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv04[j], Fv04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
    '''
            
    # set new xticls
    ax.set_xticks([0.1, 0.3, 1., 3., 10])
    xlabels = ['0.1', '0.3', '1.0', '3.0', '10']
    ax.set_xticklabels(xlabels)
    plt.xlim([0.05, 5.])
        
plt.savefig('SS14_F(T)_factors_BC.png', format='png', bbox_inches='tight', dpi=150)

plt.show()