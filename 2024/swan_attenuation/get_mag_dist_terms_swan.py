def get_distance_term(rhyp, c):
    from numpy import sqrt, log10
    '''
    c = frequency-specific coeffs
    '''
    
    if rhyp <= c['r1']:
        D1 = sqrt(rhyp**2 + c['nref']**2)
        distterm = c['nc0s'] * log10(D1) #+ c['nc1s']
    
    # set mid-field
    elif rhyp > c['r1']:
        D1 = sqrt(c['r1']**2 + c['nref']**2)
        distterm = c['nc0s'] * log10(D1) \
                   + c['mc0f'] * log10(rhyp / c['r1']) + c['mc1f'] * (rhyp - c['r1'])
    
    '''
    # set far-field
    elif rhyp > c['r2']:
        D1 = sqrt(c['r1']**2 + c['nref']**2)
        distterm = c['nc0s'] * log10(D1) \
                   + c['mc0f'] * log10(c['r2'] / c['r1']) + c['mc1fs'] * (c['r2'] - c['r1']) \
                   + c['fc0s'] * log10(rhyp / c['r2']) + c['fc1s'] * (rhyp - c['r2']) \
                   + c['fc2s'] * (log10(rhyp) - log10(c['r2']))
    '''                  
    return distterm
    
def get_magnitude_term(mw, c):
    
    return c['magc0s'] * mw + c['magc1s']
    
def parse_kappa_data():
    from numpy import array, loadtxt
    
    kapdat = []
    # read parameter file
    lines = open('site_kappa.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        kap = {'sta':dat[0], 'kappa0': float(dat[1]), 'cnt': float(dat[2])}
    
        kapdat.append(kap)
    
    return kapdat    
    
def get_kappa_term(sta, freqs):
    from numpy import exp, log10, isnan, pi
    
    # get kappas               
    kapdat = parse_kappa_data()#	get distance independent kappa
    
    kappa = kapdat[-1]['kappa0'] # default kappa
    
    # get site kappa
    for kap in kapdat:
        if kap['sta'] == sta:
            if not isnan(kap['kappa0']):
                kappa = kap['kappa0']
    
    return log10(exp(-1 * pi * freqs * kappa))