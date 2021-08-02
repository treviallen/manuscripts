
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import scipy.odr.odrpack as odrpack

import get_average_Q_list, extract_Q_freq

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x >= hx
    xmod[idx] = 1
    return xmod


def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod


def bilinear_reg_free(c, x):
    """ sets up bilinear equation with a free hinge position"""
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yhinge = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yhinge)

    return ans1 + ans2


def bilinear_reg_fix(c, x):
    """ sets up bilinear equation with a fixed hinge position"""
    from numpy import zeros_like
    hxfix = 4.25 #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hxfix)
    modx_hi = highside(x, hxfix)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hxfix + c[1]
    ans2 = modx_hi * (c[2] * (x-hxfix) + yarea)

    return ans1 + ans2


def get_freq_Q_log(station, elat, elon, evdp):
    """ run functions to get Q and freq as lists"""

    freqs = extract_Q_freq.get_freq_list()
    freqs_log = get_average_Q_list.convert_list_to_log(freqs)
    Q_list, stdevQ_list = get_average_Q_list.get_all_Qs_centre(elat, elon, station[0], station[1], evdp)
    Q_list_log = get_average_Q_list.convert_list_to_log(Q_list)
    
    return freqs_log, Q_list_log, stdevQ_list


def plot_Qs(axes, stations, elat, elon, evdp):
    """
    function plots frequecy vs Q for two average Q methods
    """
    for i,ax in enumerate(axes.flatten()):
        
        # x = freqs, returned in log10 space fro -1 to 1 (0.1 to 10 in linear)
        # y = Q, returned in log space for log domain bilinear regression
        x, y, stdev = get_freq_Q_log(stations[i], elat, elon, evdp)
        #stdev_log = get_average_Q_list.convert_list_to_log(stdev)
        y_range, x_range = bilinear_regression(x, y)
        yerr = stdev
        ax.errorbar(x[12:], y[12:], yerr[12:], c='gray', 
                            capsize=2, elinewidth=1, markeredgewidth=1)
        ax.plot(x[12:], y[12:], 'o', markersize='3', c='gray')
        ax.plot(x_range, y_range, '-', c='red', lw=2,label='Automatic Bilinear')

        ax.set_xlabel('log Frequency (Hz)')
        ax.set_ylabel('log Q')

    plt.savefig('bilinear_q.png', fmt='png', bbox_inches='tight', dpi=300)
    plt.show()


def bilinear_regression(x, y):

    data = odrpack.RealData(x[12:], y[12:])

    x_range = np.arange(-0.5, 1.0, step=0.01) # x range 

    bilin_reg = odrpack.Model(bilinear_reg_free)
    odr = odrpack.ODR(data, bilin_reg, beta0=[0.4, 3.0, 0.3, -0.5])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
    out = odr.run()
    print('\nbilinear auto\n')
    out.pprint()

    a = out.beta[0]
    b = out.beta[1]
    c = out.beta[2]
    hx = out.beta[3] # x hinge point

    y_range = b + a * x_range # get y values from bilinear
    yhinge = b + a * hx
    idx = x_range > hx
    y_range[idx] = c * (x_range[idx]-hx) + yhinge

    return y_range, x_range




#TESTING PARAMETERS
elat, elon = -25.579, 129.832
evdp = 0.0

# random locations - can go through stations and plot them 
# later, but not important at the moment.  
stations = np.array([[-22.643, 114.234],
                     [-20.557, 139.605],
                     [-13.957, 143.174],
                     [-25.037,128.296]])


fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
plot_Qs(axes, stations, elat, elon, evdp)













