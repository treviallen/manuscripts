def calc_nac_gmm_spectra(mag, rhyp, dep, region):
    from numpy import loadtxt, log10, log
    
    if region == 'BS':
        coeffile = 'ncc_gmm_coeffs.BS.csv'
    elif region == 'NGH':
        coeffile = 'ncc_gmm_coeffs.NGH.csv'
    elif region == 'OBE':
        coeffile = 'ncc_gmm_coeffs.OBE.csv'
    
    coeffs = loadtxt(coeffile, delimiter=',', skiprows=2)  
    
    T  = coeffs[:,0]
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
from numpy import array, exp, log, interp, vstack
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

print('Starting inversion...')
for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
    for i, sd in enumerate(stdict):
        pt = Point(sd['eqlo'], sd['eqla'])
        
        if zgroup == 'BS':
            mmin = 5.25
        else:
            mmin = 5.75
            
        if zgroup == 'OBW':
            zgroup = 'OBE'

        if pt.within(poly) and sd['mag'] >= mmin:
            A19imt = calc_nac_gmm_spectra(sd['mag'], sd['rhyp'], sd['dep'], zgroup)
            
            lnAmp = interp(log(A19imt['per']), log(sd['per']), log(sd['geom']))
            
            stdict[i]['lnRes'] = lnAmp - A19imt['sa']
            stdict[i]['lnSA'] = A19imt['sa']
            stdict[i]['vs30'] = get_station_vs30(sd['sta'])[0]
            
            if i == 0:
                res_stack = array([stdict[i]['lnRes']])
            else:
                res_stack = vstack((res_stack, [stdict[i]['lnRes']]))
            
            # build vs30 array
            vs30.append(stdict[i]['vs30'])
            
###############################################################################
# build residual data
###############################################################################
# set periods
Tplt = sd['per']
fig = plt.figure(1, figsize=(18,10))

for i, T in enumerate(Tplt):
    ax = plt.subplot(4,5,i+1)
    
    # get res data
    Yres = res_stack[:,i]
    
    plt.semilogx(vs30, Yres, '+', c='0.7', ms=5)
    plt.ylabel = str(T)
    
    if i > 15:
       plt.xlabel('Vs30') 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    