from os import path
from numpy import array, exp, log, nan, mean, std
from mmi_tools import mmi2pgm_worden12, mmi2pgm_atkinson07, mmi2pgm_dangkua11_ena, mmi2pgm_caprio15

'''
Replace localities file with national grid and make csvs for each MMI level

Use average MMI values to make maps
'''

w12_periods = ['pgv', 'pga', 'sa03', 'sa10'] 
a07_periods = ['pgv', 'pga', 'sa03', 'sa10', 'sa20']
c15_periods = ['pgv', 'pga']

target_mmi = array([5.0])

# non-published probabilities
probs1 = [0.0200, 0.01375, 0.0100, 0.00445, 0.0021, 0.0010, 0.0005, 0.000404, \
         0.0002, 0.0001]

# published probabilities         
probs2 = [0.0200, 0.01375, 0.0100, 0.00445, 0.0021, 0.0010, 0.0005, 0.000404]
         
years = 1
###############################################################################################
# first get PGAref 
'''
a)	0.8 PGA where the ratio Sa(0.2)/PGA < 2.0, and 
b)	1 PGA otherwise.
'''         
###############################################################################################
print 'Calculating PGA ref...'
# get 1/2475 0.2 s value
hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_SA0.2all_000thperc.sol')
lines = open(hcfile).readlines()[2:]
ref02 = []
grdref = []
lon = []
lat = []
for line in lines:
    dat = line.strip().split()
    ref02.append(float(dat[9]))
    grdref.append(dat[-1])
    lon.append(float(dat[0]))
    lat.append(float(dat[1]))

# get 1/2475 PGA value
hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_PGAall_000thperc.sol')
lines = open(hcfile).readlines()[2:]
ref00 = []
for line in lines:
        dat = line.strip().split()
        ref00.append(float(dat[9]))
        
        
pga_ref = array(ref02) / array(ref00)

###############################################################################################
# calc MMI return T
###############################################################################################
def get_mmi_return_T(line, probs1):
    from scipy.interpolate import interp1d
    from misc_tools import extrap1d
    
    dat = line.strip().split()
    lat.append(float(dat[1]))
    lon.append(float(dat[0]))
    haz_curve = array([float(x) for x in dat[2:12]])
    
    # do interpolation for each period
    f_i = interp1d(log(haz_curve), log(probs1))
    f_x = extrap1d(f_i)        
    mmi_prob = 1. / exp(f_x(log(target_gm)))[0] # MMI return period
    
    return mmi_prob
    
def get_mmi_return_T_pub(line, probs2):
    from scipy.interpolate import interp1d
    from misc_tools import extrap1d
    
    dat = line.strip().split()
    lat.append(float(dat[1]))
    lon.append(float(dat[0]))
    haz_curve = array([float(x) for x in dat[2:10]])
    
    # do interpolation for each period
    f_i = interp1d(log(haz_curve), log(probs2))
    f_x = extrap1d(f_i)        
    mmi_prob = 1. / exp(f_x(log(target_gm)))[0] # MMI return period
    
    return mmi_prob

###############################################################################################
# get W12 probs
###############################################################################################

print 'Getting W12 probabilities...'
w12_dict = {'pgv':[], 'pga':[], 'sa03':[], 'sa10':[]}
for period in w12_periods:
    # set target GM
    target_gm = mmi2pgm_worden12(target_mmi, period, nan, nan)[0]
    print period, target_gm[0]
    
    if period == 'pgv':
        target_gm /= 100.
    else:
        target_gm /= 981.
    
    prob_50 = []
    per = period.replace('sa', '')
    try:        
        # for PGA/PGV
        try:
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per.upper()+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        except:
            per = 'SA'+per[0]+'.'+per[1]
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        
        
        for line in lines:
            dat = line.strip().split()
            
            mmi_prob = get_mmi_return_T(line, probs1) # MMI return period
    
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
        
            # fill dict
            w12_dict[period].append(prob_50)
            
    # read alt file
    except:
        # for PGA/PGV
        print per
        per = period.replace('SA', '')
        per = period.replace('sa', '')
        per = 'Sa'+per[0]+'.'+per[1]
        hcfile = path.join('..','2015_grid','solfiles','GSC2015_seismichazardgrid_'+per+'.txt')
        print hcfile, per
        
        for line in lines:
            
            mmi_prob = get_mmi_return_T_pub(line, probs2) # MMI return period
            
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
    
            # fill dict
            w12_dict[period].append(prob_50)
#    except:
#        print per, 'file not available - use file with alt probabilities' 
        
###############################################################################################
# get AK07 probs
###############################################################################################

print 'Getting AK07 probabilities...'       
a07_dict = {'pgv':[], 'pga':[], 'sa03':[], 'sa10':[], 'sa20':[]}
for period in a07_periods:
    # set target GM
    target_gm = mmi2pgm_atkinson07(target_mmi, period)[0]
    print period, target_gm 
    
    if period == 'pgv':
        target_gm /= 100.
    else:
        target_gm /= 981.
    
    # parse national grid file
    prob_50 = []
    per = period.replace('sa', '')
    try:        
        # for PGA/PGV
        try:
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per.upper()+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        except:
            per = 'SA'+per[0]+'.'+per[1]
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        
        for line in lines:
            dat = line.strip().split()
            
            mmi_prob = get_mmi_return_T(line, probs1) # MMI return period
    
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
        
            # fill dict
            a07_dict[period].append(prob_50)
            
    # read alt file
    except:
        # for PGA/PGV
        per = period.replace('SA', '')
        per = period.replace('sa', '')
        per = 'Sa'+per[0]+'.'+per[1]
        hcfile = path.join('..','2015_grid','solfiles','GSC2015_seismichazardgrid_'+per+'.txt')
        
        lines = open(hcfile).readlines()[2:]
        
        for line in lines:
            
            mmi_prob = get_mmi_return_T_pub(line, probs2) # MMI return period
            
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
    
            # fill dict
            a07_dict[period].append(prob_50)
            
###############################################################################################
# do Dan-Cra11 ENA probs
###############################################################################################

print 'Getting DC11 probabilities...'       
d11_dict = {'pgv':[], 'pga':[], 'sa03':[], 'sa10':[], 'sa20':[]}
for period in a07_periods:
    # set target GM
    target_gm = mmi2pgm_dangkua11_ena(target_mmi, period, nan, nan)[0]
    print period, target_gm 
    
    if period == 'pgv':
        target_gm /= 100.
    else:
        target_gm /= 981.
    
    # parse national grid file
    prob_50 = []
    per = period.replace('sa', '')
    try:        
        # for PGA/PGV
        try:
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per.upper()+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        except:
            per = 'SA'+per[0]+'.'+per[1]
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
                
        for line in lines:
            dat = line.strip().split()
            
            mmi_prob = get_mmi_return_T(line, probs1) # MMI return period
    
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
        
            # fill dict
            d11_dict[period].append(prob_50)
            
    # read alt file
    except:
        # for PGA/PGV
        per = period.replace('SA', '')
        per = period.replace('sa', '')
        per = 'Sa'+per[0]+'.'+per[1]
        hcfile = path.join('..','2015_grid','solfiles','GSC2015_seismichazardgrid_'+per+'.txt')
        
        lines = open(hcfile).readlines()[2:]
        
        for line in lines:
            
            mmi_prob = get_mmi_return_T_pub(line, probs2) # MMI return period
            
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
    
            # fill dict
            d11_dict[period].append(prob_50)

###############################################################################################
# do Caprio et al 2015 probs
###############################################################################################
print 'Getting C15 probabilities...'       
c15_dict = {'pgv':[], 'pga':[]}
for period in c15_periods:
    # set target GM
    target_gm = mmi2pgm_caprio15(target_mmi, period)[0]
    print period, target_gm[0]	
    
    if period == 'pgv':
        target_gm /= 100.
    else:
        target_gm /= 981.
    
    # parse national grid file
    prob_50 = []
    per = period.replace('sa', '')
    try:        
        # for PGA/PGV
        try:
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per.upper()+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
        except:
            per = 'SA'+per[0]+'.'+per[1]
            hcfile = path.join('..','2015_grid','solfiles','NBCC2015_Canada_'+per+'all_000thperc.sol')
            lines = open(hcfile).readlines()[2:]
                
        for line in lines:
            dat = line.strip().split()
            
            mmi_prob = get_mmi_return_T(line, probs1) # MMI return period
    
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
        
            # fill dict
            c15_dict[period].append(prob_50)
            
    # read alt file
    except:
        # for PGA/PGV
        per = period.replace('SA', '')
        per = period.replace('sa', '')
        per = 'Sa'+per[0]+'.'+per[1]
        hcfile = path.join('..','2015_grid','solfiles','GSC2015_seismichazardgrid_'+per+'.txt')
        
        lines = open(hcfile).readlines()[2:]
        
        for line in lines:
            
            mmi_prob = get_mmi_return_T_pub(line, probs2) # MMI return period
            
            prob_50 = (1. - exp(-years/array(mmi_prob))) * 100. # in %
    
            # fill dict
            c15_dict[period].append(prob_50)
            
###############################################################################################
# get weights for mmi eqns
###############################################################################################

w12sigs = array([0.73, 0.65, 0.84, 0.80]) #pga, pgv, sa0.3, sa1.0
w12wts  = w12sigs / sum(w12sigs)
iw12wt  = 1/w12wts
iw12wt  = iw12wt/sum(iw12wt)

c15sigs = array([0.70, 0.90]) #pga, pgv
c15wts  = c15sigs / sum(c15sigs)
ic15wt  = 1/c15wts
ic15wt  = ic15wt/sum(ic15wt)

# wt W12 and C15 equally
iwt = hstack((iw12wt*0.5, ic15wt*0.5))

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    from numpy import average, sqrt
    
    wtaverage = average(values, weights=weights)
    variance = average((values-wtaverage)**2, weights=weights)  # Fast and numerically precise
    return (wtaverage, sqrt(variance))

###############################################################################################
# write to text
###############################################################################################

# write to file
txt = 'GRD_REF,LON,LAT,W12_PGV,W12_PGA,W12_SA0.3,W12_SA1.0,A07_PGV,A07_PGA,A07_SA0.3,A07_SA1.0,A07_SA2.0,DC11_ENA_PGV,DC11_ENA_PGA,DC11_ENA_SA0.3,DC11_ENA_SA1.0,DC11_ENA_SA2.0,Cea15_PGV,Cea15_PGA,MEAN_WEST,STD_WEST,MEAN_EAST,STD_EAST,SA2.0/PGA\n'
for gr in grdref:
    i = int(gr) - 1
    # calc weighted average & std western - just use W12 & C15
    meanwest = sum(array([w12_dict['pga'][i], w12_dict['pgv'][i], \
                      w12_dict['sa03'][i], w12_dict['sa10'][i], \
                      c15_dict['pga'][i], c15_dict['pgv'][i]]) * iwt)
                     
    values = array([w12_dict['pga'][i], w12_dict['pgv'][i], \
                    w12_dict['sa03'][i], w12_dict['sa10'][i], \
                    c15_dict['pga'][i], c15_dict['pgv'][i]])
        
    stdwest  = weighted_avg_and_std(values, iwt)[1]

    # calc average & std eastern - use AK07 & DC11 (PGA & PGV only)
    txt += ','.join((gr, str(lon[i]), str(lat[i]), \
                    str(w12_dict['pgv'][i]), str(w12_dict['pga'][i]), \
                    str(w12_dict['sa03'][i]), str(w12_dict['sa10'][i]), \
                    str(a07_dict['pgv'][i]), str(a07_dict['pga'][i]), \
                    str(a07_dict['sa03'][i]), str(a07_dict['sa10'][i]), str(a07_dict['sa20'][i]), \
                    str(d11_dict['pgv'][i]), str(d11_dict['pga'][i]), \
                    str(d11_dict['sa03'][i]), str(d11_dict['sa10'][i]), str(d11_dict['sa20'][i]), \
                    str(c15_dict['pgv'][i]), str(c15_dict['pga'][i]), \
                    str(meanwest), str(stdwest), str(meaneast), str(stdeast), \
                    str(pga_ref[i]))) + '\n'

f = open('MMI'+str(target_mmi[0])+'_'+str(years)+'-yr.grid.csv', 'wb')
f.write(txt)
f.close()

