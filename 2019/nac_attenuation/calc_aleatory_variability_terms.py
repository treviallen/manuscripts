def calc_nac_gmm_spectra(mag, rhyp, dep, vs30, region):
    from numpy import loadtxt, log10, log
    
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
    
    # get site coefs
    sitefile = 'nac_site_amp_coeffs.csv'
    coeffs = loadtxt(sitefile, delimiter=',', skiprows=1)  
    
    T  = coeffs[:,0]
    s0 = coeffs[:,1]
    s1 = coeffs[:,2]
    	
    site_term = s0 + s1 / (log10(vs30) - log10(150))
    
    lnsa = mag_term + atten_term + dep_term + near_field_term + site_term
           
    A19imt = {'per':T, 'sa':lnsa}

    return A19imt

def return_dep_terms(T):
    from numpy import loadtxt, log10, log, logspace, interp
    
    coeffile = 'ncc_gmm_coeffs.BS.csv'
    
    coeffs = loadtxt(coeffile, delimiter=',', skiprows=2)  
    
    Ts  = coeffs[:,0]
    c0 = coeffs[:,1]
    c1 = coeffs[:,2]
    c2 = coeffs[:,3]
    c3 = coeffs[:,4]
    c4 = coeffs[:,5]
    d0 = coeffs[:,6]
    d1 = coeffs[:,7]
    d2 = coeffs[:,8]
    d3 = coeffs[:,9]
    hx = coeffs[:,10]
    
    dep_terms = []
    deps = logspace(log10(5), log10(750), 100)
    for dep in deps:
        logdep = log10(dep)
        
        dep_term = d0 + d1*logdep**3 + d2*logdep**2 + d3*logdep
        
        dep_terms.append(interp(T, Ts, dep_term))
    
    return deps, dep_terms
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
from misc_tools import get_binned_stats, get_log_xy_locs
from numpy import array, arange, exp, log, interp, vstack, nan, isnan, log10, polyfit, \
                  percentile, logspace, unique, mean, where, std, sqrt
from misc_tools import dictlist2array 
import scipy.odr.odrpack as odrpack
from mpl_toolkits.basemap import Basemap
from scipy.stats import linregress
import pickle

print('Loading pkl file...')
stdict = pickle.load(open("stdict.pkl", "rb" ))

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
ev_stack = []
mag_stack = []
rhyp_stack = []
lat_stack = []
lon_stack = []
dep_stack = []

print('Getting residuals...')
idx = []
for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
    for i, sd in enumerate(stdict):
        pt = Point(sd['eqlo'], sd['eqla'])
        
        if zgroup == 'BS':
            mmin = 5.25
        if zgroup == 'BS':
            if pt.within(poly) and sd['mag'] >= mmin and sd['rhyp'] > 500 and sd['rhyp'] < 1500 and sd['network'] != 'OA':
                idx.append(i)
                
                # get site vs30
                vs30 = get_station_vs30(sd['sta'])[2] # use USGS
                if isnan(vs30):
                    vs30 = 450
            
                A19imt = calc_nac_gmm_spectra(sd['mag'], sd['rhyp'], sd['dep'], vs30, zgroup)
                
                lnAmp = interp(log(A19imt['per']), log(sd['per']), log(sd['geom']))
                
                stdict[i]['lnRes'] = lnAmp - A19imt['sa']
                stdict[i]['lnSA'] = A19imt['sa']
                stdict[i]['vs30'] = vs30
                
                if len(res_stack) == 0:
                    res_stack = array([stdict[i]['lnRes']])
                else:
                    res_stack = vstack((res_stack, [stdict[i]['lnRes']]))
                ev_stack.append(sd['ev'])
                mag_stack.append(sd['mag'])
                rhyp_stack.append(sd['rhyp'])
                lat_stack.append(sd['eqla'])
                lon_stack.append(sd['eqlo'])
                dep_stack.append(sd['dep'])
                
                # build vs30 array
                #vs30.append(stdict[i]['vs30'])

gmmT = A19imt['per']

stdict = stdict[array(idx)]
ev_stack = array(ev_stack)
mag_stack = array(mag_stack)
lat_stack = array(lat_stack)
lon_stack = array(lon_stack)
dep_stack = array(dep_stack)

###############################################################################
# get inter-event
###############################################################################
print('Getting inter-event term...')
unique_ev = unique(ev_stack)

inter_mag = []
inter_lon = []
inter_lat = []
inter_dep = []
inter_ev = []

for uev in unique_ev:
    idx = where(ev_stack == uev)[0]
    
    if len(inter_ev) == 0:
        inter_ev = array(mean(res_stack[idx], axis=0))
    else:
        inter_ev = vstack((inter_ev, mean(res_stack[idx], axis=0)))
    
    inter_mag.append(mag_stack[idx][0])
    inter_lon.append(lon_stack[idx][0])
    inter_lat.append(lat_stack[idx][0])
    inter_dep.append(dep_stack[idx][0])

###############################################################################
# map tau
###############################################################################
"""
urcrnrlat = -1
llcrnrlat = -12.
urcrnrlon = 135.
llcrnrlon = 120
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(3, figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=1000.)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.fillcontinents(color='w',lake_color='0.9')
m.drawmapboundary(fill_color='0.9')
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

x, y = m(inter_lon, inter_lat)

# get average event residual
res = []
norm = mpl.colors.Normalize(vmin=-1.5, vmax=1.5)

# get event term for 0.5 s
for ie, im in zip(inter_ev, inter_mag):
    res.append(interp(log(0.5), log(gmmT), ie))
m.scatter(x, y, c=res, marker='o', s=50, cmap='seismic', norm=norm)

#norm = mpl.colors.Normalize(vmin=0, vmax=250)
#m.scatter(x, y, c=inter_dep, marker='o', s=50, cmap='Spectral', norm=norm)

plt.savefig('map_inter-event_residuals.png', fmt='png', bbox_inches='tight')       

plt.show()
"""
###############################################################################
# plot tau vs dep (for 0.2 s)
###############################################################################
# get event term for 0.5 s
res = []
for ie, im in zip(inter_ev, inter_mag):
    res.append(interp(log(0.5), log(gmmT), ie))

fig = plt.figure(4, figsize=(8, 5))
plt.plot(inter_dep, res, 'rs')
plt.plot([0, 700], [0, 0], 'k--')
plt.ylim([-2, 2])
plt.show()

###############################################################################
# get tau
###############################################################################
print('Getting tau...')
do scatter with depth - show deepest are badest...
# set periods
Trng = array([0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 9., 10.]) #, 10.]) # secs; PGA = 0.01; PGV = -99
Tplt = array([0.2, 0.5, 2.0, 5.0])
fig = plt.figure(1, figsize=(14,7))
props = dict(boxstyle='round', facecolor='w', alpha=1)
tau = []

i = 0
for j, T in enumerate(Trng):
    res = []
    mag = []
    
    for ie, im in zip(inter_ev, inter_mag):
        res.append(interp(log(T), log(gmmT), ie))
        mag.append(im)
    
    tau.append(std(res))
    print('T', T, 'tau', std(res))
    
    # check if in plt periods
    idx = where(Tplt == T)[0]
    if len(idx) > 0:
        ax = plt.subplot(2,2,i+1)
    
        plt.plot(mag, res, 'o', c='dodgerblue', mec='blue', ms=6)
        plt.plot([5, 1000],[0,0], 'k--', lw=1.)
            
        if i >= 2:
           plt.xlabel('Moment Magnitude', fontsize=16)
        if i == 0 or i == 2:
           plt.ylabel('Inter-Event Terms\n(ln Residual)', fontsize=16)
        
        plt.xlim([5, 7.7])
        plt.ylim([-3, 3])
        pertxt = ' '.join(('T =', str(T), 's'))
        xpos = get_log_xy_locs([5, 7.7], 0.96)
        ypos = (6*0.92) - 3
        plt.text(xpos, ypos, pertxt, ha='right', va='top', fontsize=16, bbox=props)
        
        i += 1

plt.savefig('inter-event_residuals.png', fmt='png', bbox_inches='tight')       
plt.show()

###############################################################################
# get phi
###############################################################################

intra_ev = []
for uev, ie in zip(unique_ev, inter_ev):
    idx = where(ev_stack == uev)[0]
    
    if len(intra_ev) == 0:
        intra_ev = array(res_stack[idx] - ie)
    else:
        intra_ev = vstack((intra_ev, (res_stack[idx] - ie)))

###############################################################################
# plot phi
###############################################################################

print('Getting phi...')

# set periods
Tplt = [0.2, 0.5, 2.0, 5.0]
fig = plt.figure(2, figsize=(14,7))
props = dict(boxstyle='round', facecolor='w', alpha=1)

phi = []
i = 0
for j, T in enumerate(Trng):
    res = []
    rhyp = []
    
    for ie, rs in zip(intra_ev, rhyp_stack):
        res.append(interp(log(T), log(gmmT), ie))
        rhyp.append(rs)
    
    phi.append(std(res))
    print('T', T, 'phi', std(res))
    
    # check if in plt periods
    idx = where(Tplt == T)[0]
    if len(idx) > 0:
        ax = plt.subplot(2,2,i+1)
    
        plt.semilogx(rhyp, res, 'o', c='dodgerblue', mec='blue', ms=6)
        plt.semilogx([1, 1500],[0,0], 'k--', lw=1.)
        
        if i >= 2:
           plt.xlabel('Hypocentral Distance', fontsize=16)
        if i == 0 or i == 2:
           plt.ylabel('Intra-Event\n(ln Residual)', fontsize=16)
        
        plt.xlim([500, 1600])
        plt.ylim([-3, 3])
        pertxt = ' '.join(('T =', str(T), 's'))
        xpos = get_log_xy_locs([500, 1600], 0.04)
        ypos = (6*0.92) - 3
        plt.text(xpos, ypos, pertxt, ha='left', va='top', fontsize=16, bbox=props)
        
        i += 1

plt.savefig('intra-event_residuals.png', fmt='png', bbox_inches='tight')       
plt.show()
    
###############################################################################
# plot phi & tau with T
###############################################################################

fig = plt.figure(5, figsize=(8, 5))

sigma = sqrt(array(tau)**2 + array(phi)**2)
plt.semilogx(Trng, sigma, 'r-', lw=2, label=r'$\sigma$')
plt.semilogx(Trng, tau, 'r-', lw=1, label=r'$\tau$')
plt.semilogx(Trng, phi, 'r--', lw=1, label=r'$\phi$')

plt.legend(loc=1, fontsize=16)
plt.xlabel('Period (s)', fontsize=16)
plt.ylabel('Standard Deviations (ln Units)', fontsize=16)
plt.xlim([0.07, 10])
plt.ylim([0.3, 1.])
plt.savefig('standard_deviations.png', fmt='png', bbox_inches='tight')       
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    