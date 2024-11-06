from numpy import interp, log, exp, loadtxt, logspace, pi, sqrt, vstack, nanmean, nan, where
from misc_tools import listdir_extension, savitzky_golay, get_mpl2_colourlist
from mapping_tools import distance
from os import path
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
warnings.filterwarnings('ignore')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42


interp_freqs = logspace(-1.,2,100)[:-10]
# use best eq loc from small-mag event (https://earthquakes.ga.gov.au/event/ga2024svdpoz)
eqla = -32.34
eqlo = 150.87
eqdp = 3.0
	
def get_smoothed_fft_spectra(freqs, disp_amps):
    
    # smooth spectra
    if len(freqs) > 10000:
        sw = 201
    elif len(freqs) > 5000:
        sw = 101
    elif len(freqs) > 1000:
        sw = 51
    elif len(freqs) > 500:
        sw = 21
    else:
        sw = 5
    
    # smooth spectra
    if sw > 11:
        smoothed_disp = exp(savitzky_golay(log(disp_amps), sw, 3))
    else:
        smoothed_disp = exp(savitzky_golay(log(disp_amps), sw, 2))
        
    # log-log interpolate to "interp_freqs"
    smoothed_interp_disp = exp(interp(log(interp_freqs), log(freqs), log(smoothed_disp), \
                               left=nan, right=nan))
    
    return smoothed_disp, smoothed_interp_disp
    
fig = plt.figure(1,figsize = (13,6))
cols = get_mpl2_colourlist()

# look for fas files
fasfiles = listdir_extension('fds','fds')

# loop through fas files
zamps = []
hamps = []
zi = 0
hi = 0
for ff in fasfiles:
    fasfile = path.join('fds',ff)
    
    lines = open(fasfile).readlines()
    stlo = float(lines[4].split('\t')[-1])
    stla = float(lines[3].split('\t')[-1])
    date = lines[1].strip().split('\t')[-1][:-6]
    sta = lines[2].strip().split('\t')[-1]
    comp = lines[0].strip().split('.')[-2]
    
    rngkm = distance(eqla, eqlo, stla, stlo)[0]
    rhyp = sqrt(rngkm**2 + eqdp**2)
    
    data = loadtxt(fasfile, delimiter='\t', skiprows=22)
    freqs = data[:,0]
    vel_amps = data[:,1] * 2 * pi * freqs
    #vel_amps = data[:,1]
    	
    smoothed_vel, smoothed_interp_vel = get_smoothed_fft_spectra(freqs, vel_amps)
    if comp == 'EHZ':
        idx = where(interp_freqs < 0.5)[0]
        smoothed_interp_vel[idx] = nan
        
    smoothed_interp_vel *= 1000 # convert to mm
    
    label = ' '.join((date,sta,comp))
    
    if ff.endswith('Z.fds'):
        plt.subplot(122)
        plt.loglog(interp_freqs, smoothed_interp_vel, '-', c=cols[zi], label=label)
        zi += 1
        
        if len(zamps) == 0:
            zamps = log(smoothed_interp_vel)
        else:
            zamps = vstack((zamps, log(smoothed_interp_vel)))
    else:
        plt.subplot(121)
        plt.loglog(interp_freqs, smoothed_interp_vel, '-', c=cols[hi], label=label)
        hi += 1
        
        if len(hamps) == 0:
            hamps = log(smoothed_interp_vel)
        else:
            hamps = vstack((hamps, log(smoothed_interp_vel)))
            
mean_hamps = exp(nanmean(hamps, axis=0))
mean_zamps = exp(nanmean(zamps, axis=0))
idx = where(interp_freqs >= 0.5)[0]

# finish plots
plt.subplot(121)
plt.loglog(interp_freqs, mean_hamps, 'k-', lw=2.5, label='Mean Horizontal')
plt.legend(loc=3, fontsize=10)
plt.grid(which='both')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Horizontal Fourier Velocity Amplitude (mm)')
plt.ylim([1E-3, 20])
plt.xlim([0.1, 60])

plt.subplot(122)
plt.loglog(interp_freqs[idx], mean_zamps[idx], 'k-', lw=2.5, label='Mean Vertical')
plt.legend(loc=3, fontsize=10)
plt.grid(which='both')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Vertical Fourier Velocity Amplitude (mm)')
plt.ylim([1E-3, 20])
plt.xlim([0.1, 60])

plt.savefig('stacked_fourier_velocity.png', dpi=300, fmt='png', bbox_inches='tight')

plt.show()