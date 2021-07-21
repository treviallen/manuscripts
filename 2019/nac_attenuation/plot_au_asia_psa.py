# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder):
    print(stn)
    from fnmatch import filter
    from os import path, walk, system
    from numpy import sqrt

    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            
            # check strong-motion
            if filename.find(stn) >= 0 and filename.find('NE.psa') >= 0:
                efile = path.join(root, filename)
                tfile = efile.split('NE.psa')
                nfile = ''.join((tfile[0],'NN.psa',tfile[1]))
            # check velocity sensors
            elif filename.find(stn) >= 0 and filename.find('HE.psa') >= 0:
                efile = path.join(root, filename)
                tfile = efile.split('HE.psa')
                nfile = ''.join((tfile[0],'HN.psa',tfile[1]))
                
            # check velocity sensors
            elif filename.find(stn) >= 0 and filename.find('H1.psa') >= 0:
                efile = path.join(root, filename)
                tfile = efile.split('H1.psa')
                nfile = ''.join((tfile[0],'H2.psa',tfile[1]))
                
            '''
            # get Z file
            elif filename.find(stn) >= 0 and filename.find('NZ.psa') >= 0:
                zfile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('HZ.psa') >= 0:
                zfile = path.join(root, filename)
            '''
    try:
        try:
            # read data 
            T, SAe = read_psa(efile)
            T, SAn = read_psa(nfile)
            #print(efile, nfile)
            
            # read psa deatails
            esta, esps, erhyp, epga, epgv, mag, dep, stlo, stla = read_psa_details(efile)
            nsta, nsps, nrhyp, npga, npgv, mag, dep, stlo, stla = read_psa_details(nfile)
            
            pga = max([epga, npga]) / (1000. * 9.81)
            rhyp = erhyp
            
            # get geometric mean and convert to g
            geomean = exp((log(SAe) + log(SAn)) / 2.) / 9.81 # convert from m/s**2 to g  
            #geomean = exp(sqrt(log(SAe) * log(SAn))) / 9.81 # convert from m/s**2 to g
        # just E-comp
        except:
            # read data
            print(efile)
            T, geomean = read_psa(efile)
            geomean = geomean / 9.81
            
            # read psa deatails
            sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(efile)
            pga = pga / (1000. * 9.81)
        
    except:
        # read data
        print(zfile)
        T, geomean = read_psa(zfile)
        geomean = geomean / 9.81
        
        # read psa deatails
        sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(zfile)
        pga = pga / (1000. * 9.81)
        
    return T, geomean, pga, rhyp

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array

    lines = open(psafile).readlines()[24:]  # ignore first 23 lines

    SA = []
    T = []
    for line in lines:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA)
    
# get other details for plotting
def read_psa_details(psafile):
    lines = open(psafile).readlines()[0:24]

    for line in lines:
        dat = line.strip().split('\t')
        if line.find('STA') >= 0:
            sta = dat[1]
        if line.find('SR') >= 0:
            sps = float(dat[1].split()[0])
        if line.find('RHYP') >= 0:
            rhyp = float(dat[1].split()[0])
        if line.find('PGA') >= 0:
            pga = float(dat[1].split()[0]) / 9.81
        if line.find('PGV') >= 0:
            pgv = float(dat[1].split()[0])
        if line.find('EQMAG') >= 0:
            mag = float(dat[1].split()[0])
        if line.find('EQDEP') >= 0:
            dep = float(dat[1].split()[0])
        if line.find('STLO') >= 0:
            stlo = dat[1].split()[0]
        if line.find('STLA') >= 0:
            stla = dat[1].split()[0]
            
    return sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort, isnan
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')
from misc_tools import get_mpl2_colourlist, get_log_xy_locs
#from calc_oq_gmpes import scr_gsims, get_station_vs30, nga_east_mean, adjust_gmm_with_nga_east, inslab_gsims
from fnmatch import filter
from os import path, walk, system
from numpy import sqrt

folder = 'psa_interface'
prefix = 'hsd_spectra'
colTrue = 'True'
#colTrue = argv[3]

# set event details

ii = 1
fig = plt.figure(ii, figsize=(9, 9))
ax = plt.subplot(111)
cs = get_mpl2_colourlist()
syms = ['o', '^', 's', 'd', 'v', '<', 'h', '>', 'p']

usites = ['AU.DRS', 'AU.KDU', 'MS.KAPK', 'MY.IPM', '2009-09-30T10.16.MY.KUM', '2005-03-28T16.09.MY.KUM']
udists = []
lolatxt = ''
for i, stn in enumerate(usites):
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn) >= 0:
                if filename.find('NE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('H1') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NH') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HHH') >= 0:
                    psafile = path.join(root, filename)
                
    # get record details
    sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(psafile)
    print (sta, mag, rhyp)
    udists.append(rhyp)
    
    # get spectra
    T, geomean, pga, rhyp = get_site_geomean(stn, folder)
    
    # make label text
    part1 = psafile.split('/')[-1].split('T')[0]
    part2 = '.'.join(psafile.split('/')[-1].split('.')[2:4])
    part3 = '$\mathregular{M_W}$ '+str(mag)
    part4 = '$\mathregular{R_{hyp}}$ '+str('%0.0f' % rhyp)+' km'
    
    label = '; '.join((part1+' '+part2, part3, part4))
    
    # now plot!
    plt.loglog(T, geomean, syms[i], ls='-', c=cs[i], lw=2, \
               ms=8, mec=cs[i], mfc='w', mew=2, markevery=15, label=label)

plt.grid(which='both')
plt.legend(loc=3, fontsize=15, numpoints=1) 
plt.xlabel('Period (s)', fontsize=20)
plt.ylabel('Spectral Acceleration (g)', fontsize=20) 
ax.tick_params(labelsize=18)  

plt.savefig('figures/au_asia_spectra.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

