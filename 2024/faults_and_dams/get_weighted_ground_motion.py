# start of plotting routine
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp, unique, vstack, nan
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from calc_oq_gmpes import nsha23_gsims
from misc_tools import get_log_xy_locs, get_mpl2_colourlist

# set site details
vs30 = 760.
dep = 10.
dip = 45.
rake = 90.
ztor = 0
###################################################################
# parse lines here
csvfile = 'Data_for_Trevor_16102024.csv'
lines = open(csvfile).readlines()

rrup = []
mag = []
dom = []
for line in lines[1:]:
    dat = line.strip().split(',')
    rrup.append(float(dat[3])/1000)
    mag.append(float(dat[4]))
    dom.append(dat[7])

###################################################################

rjb = rrup

Tplot = [0]
# wts for [A12, AB06, D15, S09nc, ngae] 
nsha23nc_wts = [0.29, 0.15, 0.17, 0.29, 0.10]

# wts for [A12, AB06, D15, ESHM20, NGAE, S09YC, S09NC] 
nsha23c_wts = [0.24, 0.13, 0.14, 0.10, 0.09, 0.16, 0.14]

txtnc = ','.join(lines[0].split(',')[0:8]) + ',A12_med,AB06_med,D15_med,S09nc_med,NGAE_med,A12_std,AB06_std,D15_std,S09nc_std,NGAE_std,A12_wt,AB06_wt,D15_wt,S09nc_wt,NGAE_wt\n'
txtc  = ','.join(lines[0].split(',')[0:8]) + ',A12_med,AB06_med,D15_med,ESHM_med,NGAE_med,S09yc_med,S09nc_med,A12_std,AB06_std,D15_std,ESHM_std,NGAE_std,S09yc_std,S09nc_std,A12_wt,AB06_wt,D15_wt,ESHM_wt,NGAE_wt,S09yc_wt,S09nc_wt\n'

# loop thru periods
for j, t in enumerate(Tplot):
    for i,r in enumerate(rrup):

        # get ground motion estimates from GMPEs
        vs30 = 760.
        AB06imt, Sea09imt, Sea09YCimt, A12imt, D15imt, NGAEimt, ESHM20Cimt \
                 = nsha23_gsims(mag[i], dep, ztor, dip, rake, rrup[i], rjb[i], vs30)
        
        if t == 0.0:
            if dom[i] == '1':
                datArray = [A12imt['pga'][0][0], AB06imt['pga'][0][0], D15imt['pga'][0][0], ESHM20Cimt['pga'][0][0], NGAEimt['sa'][0], Sea09YCimt['pga'][0][0], Sea09imt['pga'][0][0], \
                            A12imt['sig'][0], AB06imt['sig'][0], D15imt['sig'][0], ESHM20Cimt['sig'][0], NGAEimt['sig'][0], Sea09YCimt['sig'][0], Sea09imt['sig'][0]]
                newline = ','.join(lines[i+1].split(',')[0:8]) + ',' + ','.join(str(x) for x in datArray) + ',' + ','.join(str(x) for x in nsha23c_wts) + '\n'
                txtc += newline
            else:
                datArray = [A12imt['pga'][0], AB06imt['pga'][0], D15imt['pga'][0], Sea09imt['pga'][0], NGAEimt['sa'][0], \
                            A12imt['sig'][0], AB06imt['sig'][0], D15imt['sig'][0], Sea09imt['sig'][0], NGAEimt['sig'][0]]
                newline = ','.join(lines[i+1].split(',')[0:8]) + ',' + ','.join(str(x) for x in datArray) + ',' + ','.join(str(x) for x in nsha23nc_wts) + '\n'
                txtnc += newline
            
        else:
            
            A12r = interp(t, A12imt['per'], A12imt['sa'])
            A12s = interp(t, A12imt['per'], A12imt['sig'])
            AB06r = interp(t, AB06imt['per'], AB06imt['sa'])
            AB06s = interp(t, AB06imt['per'], AB06imt['sig'])
            Sea09NCr = interp(t, Sea09imt['per'], Sea09imt['sa'])
            Sea09NCs = interp(t, Sea09imt['per'], Sea09imt['sig'])
            Sea09YCr = interp(t, Sea09YCimt['per'], Sea09YCimt['sa'])
            Sea09YCs = interp(t, Sea09YCimt['per'], Sea09YCimt['sig'])
            D15r = interp(t, D15imt['per'], D15imt['sa'])
            D15s = interp(t, D15imt['per'], D15imt['sig'])
            NGAEr = interp(t, NGAEimt['per'], NGAEimt['sa'])
            NGAEs = interp(t, NGAEimt['per'], NGAEimt['sig'])
            ESHM20Cr = interp(t, ESHM20Cimt['per'], ESHM20Cimt['sa'])
            ESHM20Cs = interp(t, ESHM20Cimt['per'], ESHM20Cimt['sig'])
            
            if dom[i] == '1':
                datArray = [A12r, AB06r, D15r, ESHM20Cr, NGAEr, Sea09YCr, Sea09r, \
                            A12s, AB06s, D15s, ESHM20Cs, NGAEs, Sea09YCs, Sea09s]
                newline = ','.join(lines[i+1].split(',')[0:8]) + ',' + ','.join(str(x) for x in datArray) + ',' + ','.join(str(x) for x in nsha23c_wts) + '\n'
                txtc += newline
            else:
                datArray = [A12r, AB06r, D15r, Sea09r, NGAEr, \
                            A12s, AB06s, D15s, Sea09s, NGAEs]
                newline = ','.join(lines[i+1].split(',')[0:8]) + ',' + ','.join(str(x) for x in datArray) + ',' + ','.join(str(x) for x in nsha23nc_wts) + '\n'
                txtnc += newline
            
    
    '''
    # get mean NSHA model
    wtmods = nsha_wts * exp(vstack((Aea16r, A12r, AB03r, AB06r, AB11r, Bea14r, Sea09r)))
    wtmods = wtmods.sum(axis=0)
    '''
    
'''
write files
'''
if t == 0.0:
	imt = 'PGA'
else:
	imt = '1.0'
	
f = open('nsha23_gmms_noncratonic_dams_'+imt+'.csv', 'w')
f.write(txtnc)
f.close()

f = open('nsha23_gmms_cratonic_dams_'+imt+'.csv', 'w')
f.write(txtc)
f.close()
    