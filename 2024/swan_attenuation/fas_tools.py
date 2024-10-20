from numpy import array, arange, nan, exp, log, log10, pi, logspace, interp


# get allen_etal_2006 for horizontal-component GMs in SWWA
def get_allen_etal_2006(mw, rhyp):
    '''
    # coeffs are: f, c1, c2, c3, c4, sigma, Q
    # Values describe Fourier acceleration spectra in (mm/sec)
    '''
    a06_coeffs = array([[0.78, 1.115, 1.489, 0.0390, 0.00000, 0.27, nan],
                        [0.88, 1.244, 1.477, 0.0206, 0.00000, 0.24, nan],
                        [1.07, 1.357, 1.461, 0.0044, 0.00089, 0.23, 469],
                        [1.17, 1.466, 1.437, -0.0106, 0.00094, 0.23, 485],
                        [1.47, 1.571, 1.410, -0.0227, 0.00108, 0.20, 526],
                        [1.66, 1.650, 1.375, -0.0310, 0.00117, 0.19, 551],
                        [1.95, 1.720, 1.337, -0.0343, 0.00130, 0.18, 585],
                        [2.34, 1.772, 1.316, -0.0324, 0.00146, 0.18, 626],
                        [2.83, 1.817, 1.304, -0.0259, 0.00164, 0.21, 672],
                        [3.32, 1.844, 1.298, -0.0164, 0.00182, 0.24, 712],
                        [3.91, 1.854, 1.295, -0.0056, 0.00201, 0.28, 757],
                        [4.59, 1.852, 1.282, 0.0039, 0.00223, 0.30, 803],
                        [5.47, 1.840, 1.242, 0.0107, 0.00249, 0.31, 857],
                        [6.45, 1.822, 1.195, 0.0146, 0.00276, 0.30, 911],
                        [7.62, 1.806, 1.151, 0.0157, 0.00307, 0.28, 969],
                        [9.08, 1.810, 1.114, 0.0152, 0.00342, 0.30, 1034],
                        [10.74, 1.810, 1.074, 0.0141, 0.00381, 0.27, 1100],
                        [12.70, 1.814, 1.035, 0.0149, 0.00423, 0.35, 1170],
                        [15.04, 1.801, 1.001, 0.0176, 0.00471, 0.33, 1246],
                        [17.87, 1.767, 0.984, 0.0263, 0.00525, 0.32, 1328],
                        [21.09, 1.718, 0.935, 0.0354, 0.00582, 0.30, 1412],
                        [25.00, 1.624, 1.007, 0.0561, 0.00648, 0.28, 1504]])
                        
                        
    freqs = a06_coeffs[:,0]
    c1 = a06_coeffs[:,1]
    c2 = a06_coeffs[:,2]
    c3 = a06_coeffs[:,3]
    c4 = a06_coeffs[:,4]
    
    if rhyp <= 80.:
        log_fas = c1 + c2*(mw - 4.0) + c3*(mw - 4.0)**2 \
                - log10(rhyp) - c4 * rhyp
    else:
        log_fas = c1 + c2*(mw - 4.0) + c3*(mw - 4.0)**2 \
                - log10(80.0) - 0.5 * log10(rhyp/80.0) - c4 * rhyp
                
    fas = 10**log_fas # in mm/s
    
    # apply rough H2V correction
    fas /= 10**0.18 
    
    # convert to displacement in mm-s
    fds = fas / (2 * pi * freqs)**2
    
    return freqs, fas, fds # approx vertical comp
    
# get a range of distances for Allen et al 2006 for SWWA in mm-s
def get_dist_atten_allen_2006(mw, rhyps, f_interp):
    '''
    f_interp = interpolation frequency
    rhyps = array of distances
    '''
    
    # loop through distances
    interp_fds = []
    for rhyp in rhyps:
        freqs, fas, fds = get_allen_etal_2006(mw, rhyp)
        
        interp_fds.append(exp(interp(log(f_interp), log(freqs), log(fds))))
    
    interp_fds = array(interp_fds)
    
    return interp_fds
    
'''
Get Allen 2007 for SEA
'''

# get allen_etal_2007 for vertical-component GMs in SEA
def get_allen_etal_2007(mw, rhyp):
    '''
    # coeffs are: f, c1, c2, c3, c4, sigma, Q
    # Values describe Fourier acceleration spectra in (mm/sec)
    '''
    a07_coeffs = array([[0.78, 1.113, 1.438, 0.1643, 0.00006, 0.35, 5211],
                        [0.98, 1.295, 1.457, 0.1482, 0.00009, 0.32, 4082],
                        [1.37, 1.505, 1.470, 0.0902, 0.00018, 0.27, 2880],
                        [1.56, 1.585, 1.453, 0.0521, 0.00023, 0.26, 2527],
                        [1.95, 1.731, 1.390, -0.0105, 0.00036, 0.25, 2061],
                        [2.54, 1.907, 1.287, -0.0529, 0.00058, 0.25, 1665],
                        [3.13, 2.030, 1.203, -0.0743, 0.00082, 0.26, 1440],
                        [4.10, 2.172, 1.079, -0.1144, 0.00126, 0.27, 1237],
                        [5.08, 2.255, 0.967, -0.1499, 0.00170, 0.28, 1135],
                        [6.45, 2.341, 0.901, -0.1487, 0.00228, 0.29, 1073],
                        [8.01, 2.408, 0.777, -0.2059, 0.00286, 0.30, 1062],
                        [9.96, 2.458, 0.739, -0.2131, 0.00345, 0.31, 1095],
                        [12.70, 2.443, 0.726, -0.1783, 0.00402, 0.33, 1197],
                        [15.82, 2.358, 0.690, -0.1604, 0.00438, 0.35, 1368],
                        [19.92, 2.196, 0.667, -0.1316, 0.00453, 0.38, 1667]])
                        
    freqs = a07_coeffs[:,0]
    c1 = a07_coeffs[:,1]
    c2 = a07_coeffs[:,2]
    c3 = a07_coeffs[:,3]
    c4 = a07_coeffs[:,4]
    
    if rhyp <= 90.:
        log_fas = c1 + c2*(mw - 4.0) + c3*(mw - 4.0)**2 \
                - 1.3 * log10(rhyp) - c4 * rhyp
    elif rhyp > 90 and rhyp <= 160:
        log_fas = c1 + c2*(mw - 4.0) + c3*(mw - 4.0)**2 \
                - 1.3 * log10(90.0) + 0.1 * log10(rhyp/90.0) - c4 * rhyp
                
    else:
        log_fas = c1 + c2*(mw - 4.0) + c3*(mw - 4.0)**2 \
                - 1.3 * log10(90.0) + 0.1 * log10(160.0/90.0) \
                - 1.6 * log10(rhyp/160.0) - c4 * rhyp
                
    fas = 10**log_fas # in mm/s
    
    # convert to displacement in mm-s
    fds = fas / (2 * pi * freqs)**2
    
    return freqs, fas, fds
    
# get a range of distances for Allen et al 2006 for SEA in mm-s
def get_dist_atten_allen_2007(mw, rhyps, f_interp):
    '''
    f_interp = interpolation frequency
    rhyps = array of distances
    '''
    
    # loop through distances
    interp_fds = []
    for rhyp in rhyps:
        freqs, fas, fds = get_allen_etal_2007(mw, rhyp)
        
        interp_fds.append(exp(interp(log(f_interp), log(freqs), log(fds))))
    
    interp_fds = array(interp_fds)
    
    return interp_fds

'''
Get Atkinson 2004
'''

def get_atkinson_2004(mw, rhyp):
    '''
    # coeffs are: f, c1, c2, c3, c4, sigma, Q
    # Values describe Fourier acceleration spectra in (cm/sec)
    # coeffs from errrata
    '''
    a04_coeffs = array([[0.20, -0.118, 0.986, 0.2294, 0.00000, 0.46, 241, nan],
                        [0.25, -0.122, 1.066, 0.2658, 0.00000, 0.39, 253, nan],
                        [0.32, -0.073, 1.166, 0.2493, -0.00003, 0.38, 291, 3388],
                        [0.40, 0.025, 1.605, 0.4199, -0.00006, 0.32, 358, 2424],
                        [0.50, 0.144, 1.773, 0.3918, -0.00009, 0.31, 443, 1955],
                        [0.63, 0.257, 1.746, 0.2559, -0.00017, 0.26, 555, 1378],
                        [0.79, 0.427, 1.782, 0.2148, -0.00026, 0.24, 702, 1140],
                        [1.00, 0.580, 1.619, 0.0978, -0.00035, 0.22, 1313, 1055],
                        [1.26, 0.775, 1.560, 0.0336, -0.00045, 0.22, 1532, 1004],
                        [1.59, 0.926, 1.525, 0.0204, -0.00054, 0.24, 1620, 1058],
                        [2.00, 1.031, 1.459, 0.0075, -0.00063, 0.25, 1693, 1161],
                        [2.51, 1.183, 1.411, -0.0090, -0.00080, 0.26, 1679, 1155],
                        [3.16, 1.267, 1.326, -0.0145, -0.00096, 0.27, 1651, 1219],
                        [3.98, 1.349, 1.245, -0.0199, -0.00118, 0.28, 1577, 1250],
                        [5.01, 1.442, 1.193, -0.0185, -0.00144, 0.29, 1523, 1284],
                        [6.31, 1.476, 1.111, -0.0185, -0.00164, 0.31, 1450, 1425],
                        [7.94, 1.504, 1.033, -0.0085, -0.00185, 0.32, 1383, 1588],
                        [10.00, 1.524, 0.941, -0.0164, -0.00204, 0.32, 1277, 1811],
                        [12.59, 1.497, 0.918, 0.0248, -0.00214, 0.33, 880, 2176],
                        [15.85, 1.487, 0.837, 0.0163, -0.00244, 0.35, 862, 2401],
                        [19.95, 1.302, 0.671, -0.0350, -0.00271, 0.49, 821, 2722]])
                       
    freqs = a04_coeffs[:,0]
    c1 = a04_coeffs[:,1]
    c2 = a04_coeffs[:,2]
    c3 = a04_coeffs[:,3]
    c4 = a04_coeffs[:,4]
    	
    # set m1 
    m1 = 0.36 + 0.91 * mw # Eqn 12
    
    if rhyp <= 70.:
        log_fas = c1 + c2*(m1 - 4.0) + c3*(m1 - 4.0)**2 \
                - 1.3 * log10(rhyp) + c4 * rhyp
    elif rhyp > 70 and rhyp <= 140:
        log_fas = c1 + c2*(m1 - 4.0) + c3*(m1 - 4.0)**2 \
                - 1.3 * log10(70.0) + 0.2 * log10(rhyp/70.0) + c4 * rhyp
                
    else:
        log_fas = c1 + c2*(m1 - 4.0) + c3*(m1 - 4.0)**2 \
                - 1.3 * log10(70.0) + 0.2 * log10(140.0/70.0) \
                - 0.5 * log10(rhyp/140.0) + c4 * rhyp
                
    fas = 10 * (10**log_fas) # in mm/s
    
    # convert to displacement in mm-s
    fds = fas / (2 * pi * freqs)**2
    
    return freqs, fas, fds
    
# get a range of distances for Allen et al 2006 for SEA in mm-s
def get_dist_atten_atkinson_2004(mw, rhyps, f_interp):
    '''
    f_interp = interpolation frequency
    rhyps = array of distances
    '''
    
    # loop through distances
    interp_fds = []
    for rhyp in rhyps:
        freqs, fas, fds = get_atkinson_2004(mw, rhyp)
        
        interp_fds.append(exp(interp(log(f_interp), log(freqs), log(fds))))
    
    interp_fds = array(interp_fds)
    
    return interp_fds

'''
Get Atkinson & Boore (2014)
'''
def get_atkinson_boore(mw, rhyp):
    from misc_tools import mw2m0
    
    freqs = logspace(-1,log10(25), 30) 
    vs = 3.7 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    #rho = 2.8
    C = 4. * pi * rho * vsm**3 * 1000 / (0.55 * 2.0 * 0.71)
    
    m0 = mw2m0(mw) # in N-m
    sd = 6*10E6 # MPa
    r0 = (7. * m0 / (16. * sd))**(1/3)
    f0 = 2.34 * vsm / (2 * pi * r0)

    omega0 = m0 / C
    src_spec = omega0 / (1 + (freqs / f0)**2)
        
    k0 = 0.005
    k_term = exp(-1 * pi * freqs * k0)
    
    # q term
    qf = 525 * freqs**0.45
    q_term = -(pi * freqs) / (2.3 * qf * vs)
    #print(q_term)
    
    #fds = src_spec
    
    if rhyp <= 50.:
        log_fds = log10(src_spec) - 1.3 * log10(rhyp) + q_term*rhyp
    else:
        log_fds = log10(src_spec) - 1.3 * log10(50.0) - 0.5 * log10(rhyp/50.0) + q_term*rhyp
        
    fds = 10**log_fds * 10 # convert from cm-s to mm-s
    
    return freqs, fds
    
# get a range of distances for Atkinson & Boore (2014) in mm-s
def get_dist_atten_atkinson_boore_2014(mw, rhyps, f_interp):
    '''
    f_interp = interpolation frequency
    rhyps = array of distances
    '''
    
    # loop through distances
    interp_fds = []
    for rhyp in rhyps:
        freqs, fds = get_atkinson_boore(mw, rhyp)
        
        interp_fds.append(exp(interp(log(f_interp), log(freqs), log(fds))))
    
    interp_fds = array(interp_fds) 
    
    return interp_fds