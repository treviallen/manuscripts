def calc_nac_gmm_spectra(mag, rhyp, dep, region, pgmTrue):
    from numpy import loadtxt, log10, log
    
    if pgmTrue == True:
        if region == 'BS':
            coeffile = 'ncc_pgm_gmm_coeffs.BS.csv'
        elif region == 'NGH':
            coeffile = 'ncc_pgm_gmm_coeffs.NGH.csv'
        elif region == 'OBE':
            coeffile = 'ncc_pgm_gmm_coeffs.OBE.csv'
    else:
        if region == 'BS':
            coeffile = 'ncc_gmm_coeffs.BS.csv'
        elif region == 'NGH':
            coeffile = 'ncc_gmm_coeffs.NGH.csv'
        elif region == 'OBE':
            coeffile = 'ncc_gmm_coeffs.OBE.csv'
    
    coeffs = loadtxt(coeffile, delimiter=',', skiprows=2)  
    
    T  = coeffs[:,0] #[2:-2]
    c0 = coeffs[:,1] #[2:-2]
    c1 = coeffs[:,2] #[2:-2]
    c2 = coeffs[:,3] #[2:-2]
    c3 = coeffs[:,4] #[2:-2]
    c4 = coeffs[:,5] #[2:-2]
    d0 = coeffs[:,6] #[2:-2]
    d1 = coeffs[:,7] #[2:-2]
    d2 = coeffs[:,8] #[2:-2]
    d3 = coeffs[:,9] #[2:-2]
    n0 = coeffs[:,10]
    hx = coeffs[:,11] #[2:-2]
    
    logdep = log10(dep)
#    lnsa = c0 + c1*(mag-6)**2 + c2*(mag-6) - c3*log10(rhyp) - c4*rhyp # \
#               + (d0 + d1*logdep**3 + d2*logdep**2 + d3*logdep)
    '''
    Rhyp <= hx: ln Y = c0 + c1*(M-6)**2 + c2*(M-6) + \
    (c3*hx +  c4*(log10(Rhyp)-hx)) + (d0 + d1*log10(h)**3 + d2*log10(h)**2 + d3*log10(h))
    '''
    mag_term = c0 + c1*(mag-6)**2 + c2*(mag-6)
    
    if log10(rhyp) < hx[0]:
        atten_term = c3 * log10(rhyp)
    else:
        hy = c3 * hx[0]
        atten_term = c4 * (log10(rhyp)-hx[0]) + hy
        
    near_field_term = n0 * (log10(rhyp)-hx[0])
    
    dep_term = d0 + d1*logdep**3 + d2*logdep**2 + d3*logdep
    
    lnsa = mag_term + atten_term + dep_term
           
    A19imt = {'per':T, 'sa':lnsa}

    return A19imt


################################################################################
# parse shapefile and filter stdict, and regress
################################################################################
'''
BS = Banda Sea
NGH = New Guinea Highlands
OBW = Oceanic Basin-Extended Margin West
OBE = Oceanic Basin-Extended Margin East
'''
# load shape
import shapefile
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
from calc_oq_gmpes import get_station_vs30
from misc_tools import get_binned_stats, get_binned_stats_meanx
from numpy import array, arange, exp, log, interp, vstack, nan, isnan, log10, polyfit, isfinite, where
from misc_tools import savitzky_golay
import scipy.odr.odrpack as odrpack
from scipy.stats import linregress
import pickle
from sys import argv

pgmTrue = argv[1] # True if calculating pga & pgv coeffs
if pgmTrue == 'True':
    pgmTrue = True
else:
    pgmTrue = False

print('Loading pkl file...')
stdict = pickle.load(open("stdict_ampfact.pkl", "rb" ))

shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')
zone_group = get_field_data(sf, 'ZONE_GROUP', 'str')

# loop thru zones
i = 0
vs30 = []
res_stack = []

print('Getting residuals...')
for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
    for i, sd in enumerate(stdict):
        if sd['rhyp'] < 2000.:
            pt = Point(sd['eqlo'], sd['eqla'])
            
            if zgroup == 'BS':
                mmin = 5.25
            else:
                mmin = 5.75
                
            if zgroup == 'OBW':
                zgroup = 'OBE'
            
            if zgroup == 'BS': # or zgroup == 'NGH':
                if pt.within(poly) and sd['mag'] >= mmin and sd['rhyp'] > 500 and sd['rhyp'] < 1600:
                    A19imt = calc_nac_gmm_spectra(sd['mag'], sd['rhyp'], sd['dep'], zgroup, pgmTrue)
                    
                    if pgmTrue == True:
                        lnAmp = array([log(sd['pgv']), log(sd['pga'])])
                    else:
                        lnAmp = interp(log(A19imt['per']), log(sd['per']), log(sd['geom']))
                    
                    stdict[i]['lnRes'] = lnAmp - A19imt['sa']
                    stdict[i]['lnSA'] = A19imt['sa']
                    # returns: vs30, isproxy, usgsvs, asscmvs, kvs, stla, stlo
                    stdict[i]['vs30'] = get_station_vs30(sd['sta'])[2]
                    
                    if len(res_stack) == 0:
                        res_stack = array([stdict[i]['lnRes']])
                    else:
                        res_stack = vstack((res_stack, [stdict[i]['lnRes']]))
                    
                    # build vs30 array
                    vs30.append(stdict[i]['vs30'])

vs30 = array(vs30)            
###############################################################################
# build residual data
###############################################################################
print('Plotting...')

# set periods
Tplt = A19imt['per']
fig = plt.figure(1, figsize=(18,10))

d0_array = []
d1_array = []
d2_array = []
d3_array = []
coefs1 = []
t_medx = []
t_logmedamp = []

for i, T in enumerate(Tplt):
    
    # get res data
    Yres = res_stack[:,i]
    
    if i < 25:
        ax = plt.subplot(5,5,i+1)	
        plt.semilogx([100, 1000],[0,0], 'k--', lw=0.5)
        
        plt.semilogx(vs30, Yres, '+', c='0.7', ms=5)
        plt.ylabel(str(T))
        
        if i >= 15:
           plt.xlabel('Vs30')
           
        plt.ylim([-4, 4])
        plt.xlim([200, 1000])
        xticks = [200, 500, 1000]
        ax.set_xticks(xticks)
        ax.set_xticklabels(list([str(x) for x in xticks]))
    
    # bin stats
    bins = arange(2.1, 3, 0.075)
    #logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(vs30), Yres)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats_meanx(bins, log10(vs30), Yres) # medx = meanx
    
    #t_medx.append(medx)
    t_medx.append(medx)
    t_logmedamp.append(logmedamp)
    
    # fit cubic
    #d1, d2, d3, d0 = polyfit(array(medx), array(logmedamp), 3) # binned data
    #d0_array.append(d0)
    #d1_array.append(d1)
    #d2_array.append(d2)
    #d3_array.append(d3)
    
    plt.semilogx(10**medx, logmedamp, 'rs', ms=6)
    
    ###############################################################################
    # fit data
    ###############################################################################
    def fit_site_amp(c, x):
        return c[0] + c[1] / (x - log10(150))
    
    #data = odrpack.RealData(medx, logmedamp)
    idx = where(isfinite(Yres))[0]
    data = odrpack.RealData(log10(vs30[idx]), Yres[idx])
    
    sitefit = odrpack.Model(fit_site_amp)
    odr = odrpack.ODR(data, sitefit, beta0=[0.1, 0.])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    c = out.beta
    print(c)
    coefs1.append(c)
    
    # plot fit
    dx = arange(2., 3, 0.01)
    dy = c[0] + c[1] / (dx - log10(150))
    #dy = d0 + d1*dx**3 + d2*dx**2 + d3*dx
    #print(dx, dy)
    plt.plot(10**dx, dy, 'g-', lw=1.5)

###############################################################################
# regress coeffs
###############################################################################
sg_window = 11
sg_poly = 2

#linreg = linregress(log10(Tplt), log10(array(coefs1)[:,1]))
if pgmTrue == True:
    smooth_c1 = array(coefs1)[:,1] # slope - do not smooth
else:
    smooth_c1 = savitzky_golay(array(coefs1)[:,1], sg_window, sg_poly) # slope

# refit corner
refit_c = []
coefs2 = []
for i, T in enumerate(Tplt):
    #c1fix = 10**(linreg[1] + linreg[0]*log10(T))
    c1fix = smooth_c1[i]
        
    def fit_site_amp(c, x):
        #print(yoff)
        return  c[0] + c1fix / (x - log10(150))
        
    data = odrpack.RealData(t_medx[i], t_logmedamp[i])
    
    sitefit = odrpack.Model(fit_site_amp)
    odr = odrpack.ODR(data, sitefit, beta0=[0.1])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    c = out.beta
    print(c)
    coefs2.append(c)
    
    # plot fitted fit
    dy = c[0] + c1fix / (dx - log10(150))
    
    if i < 25:
        ax = plt.subplot(5,5,i+1)
        # plot fit
        
        plt.plot(10**dx, dy, 'r-', lw=1.5)
    
plt.savefig('nac_attenuation_site_amp.png', fmt='png', bbox_inches='tight')       
plt.show()

###############################################################################
# plot coefs
###############################################################################

fig = plt.figure(2, figsize=(12, 7))
ax = plt.subplot(1,2,1)

plt.semilogx(Tplt, array(coefs1)[:,0], 'ro')
plt.semilogx(Tplt, array(coefs2)[:,0], 'bo') # re-fitted coefs
	
ax = plt.subplot(1,2,2)

plt.loglog(Tplt, array(coefs1)[:,1], 'ro')
plt.semilogx(Tplt, smooth_c1, 'go')
	
#c1fix = 10**(linreg[1] + linreg[0]*log10(Tplt))
#plt.loglog(Tplt, c1fix, 'k-')

plt.savefig('site_amp_coeffs.png', fmt='png', bbox_inches='tight')
	
plt.show()     
###############################################################################
# write coefs
###############################################################################

txt = 'NAC Site Class Model: c0 + c1 / (log10(VS30) - log10(150))\n'

for i, t in enumerate(Tplt):
    txt += ','.join((str(t), str('%0.5f' % coefs2[i]), str('%0.5f' % smooth_c1[i]))) + '\n'

if pgmTrue == True:
    f = open('nac_pgm_site_amp_coeffs.csv', 'w')
else:
    f = open('nac_site_amp_coeffs.csv', 'w')
f.write(txt)
f.close()    
