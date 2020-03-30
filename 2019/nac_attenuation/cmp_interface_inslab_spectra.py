def calc_nac_gmm_spectra(mag, rhyp, dep, vs30, region):
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
    
    # get site coefs
    sitefile = 'nac_site_amp_coeffs.csv'
    coeffs = loadtxt(sitefile, delimiter=',', skiprows=1)  
    
    T  = coeffs[:,0]
    s0 = coeffs[:,1]
    s1 = coeffs[:,2]
    	
    site_term = s0 + s1 / (log10(vs30) - log10(150))
    
    lnsa = mag_term + atten_term + dep_term + site_term
           
    A19imt = {'per':T, 'sa':lnsa}

    return A19imt


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

# calculate GMPEs and make subplots
from misc_tools import get_mpl2_colourlist
colours = get_mpl2_colourlist()
def makesubplt(colidx, fig, plt, sta, sps, mag, dep, ztor, dip, rake, rhyp, vs30, datestr, netsta):
    import matplotlib 
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8)
    
    col = colours[colidx]
    
    rrup = rhyp
    rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
    
    # get recorded process_waves.py psa data
    T, geomean, pga, rhyp = get_site_geomean(sta, folder)
    #label = '; '.join((datestr + ' ' + netsta, r'$\mathregular{M_{W}}$ '+str(mag), r'$\mathregular{R_{hyp}}$ '+str(rhyp)+' km'))
    label = '; '.join((datestr + ' ' + netsta, '$M_W$ '+str(mag), '$R_{hyp}$ '+str(rhyp)+' km'))
    plt.loglog(T, geomean, lw=2.0, color=col, label=label)

    if i == 1:
        plt.xlabel('Period (s)', fontsize=14)
        plt.ylabel('Spectral Acceleration (g)', fontsize=14)
            
        plt.xlim([0.02, 10])
        plt.grid(which='both', color='0.75')

    #plt.legend(['Yea97', 'AB06','A12','A19 (BS)', 'A19 (NGH)', 'A19 (OB)', 'Data'],loc=3, fontsize=7.)
    return plt

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort, hstack
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

from calc_oq_gmpes import inslab_gsims, scr_gsims, get_station_vs30
from gmt_tools import cpt2colormap, remove_last_cmap_colour

folder = 'psa_interface'

# set event details
'''
mag  = 5.0
dep = 11.
'''
ztor = 10. # guess
rake = 90. # USGS CMT
dip  = 30.

# set site details
#vs30 = 760

ii = 1
fig = plt.figure(ii, figsize=(7, 7))
#cmap = plt.cm.get_cmap('Spectral', 7)
ncols = 7
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
#cptfile = 'U:\\DATA\\GMT\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

from fnmatch import filter
from os import path, walk, system

# build sites
sites = []
ufiles = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.find('psa') >= 0:
            tmpsite = filename.split('.')[3]
            tmpcomp = filename.split('.')[4]
            #sites.append('.'.join((tmpsite,tmpcomp)))
            sites.append(tmpsite)
            ufiles.append(filename[0:-8])

ufiles = unique(array(ufiles))

# if two H components, rename channel
#usites = unique(sites)
'''
usites = []
for site1 in sites:
    cnt = 0
    for site2 in sites:
        if site1[0:-1] == site2[0:-1]:
            cnt += 1
            
    if cnt == 2:
        usites.append(site1[0:-1]+'H')
    else:
        usites.append(site1)
'''        
#usites = unique(sites)
usites = array(sites)
udists = []
datelist = []
lolatxt = ''
for stn in ufiles:
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
    udists.append(rhyp)
    lolatxt += '\t'.join((stlo, stla, sta))+'\n'
    datelist.append(psafile.split('/')[1][0:10])

f = open('staloc.txt', 'w')
f.write(lolatxt)
f.close()
    
# now sort by distance
udists = array(udists)
idx=argsort(array(datelist))
udists = udists[idx]
usites = usites[idx]

ufiles = ufiles[idx]
ufiles = hstack((ufiles[0], ufiles[-1], ufiles[2:-1], ufiles[1]))

# loop thru sites ordered by distance
i = 0
for stn in ufiles:
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
    print(stn, psafile)
    datestr = psafile.split('/')[1][0:10]
    netsta = '.'.join((psafile.split('/')[1].split('.')[2:4]))
    sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(psafile)
    
    # temp fix
    #mag = 4.57 # 2012-07-20
    #mag = 5.12 # 2012-06-19
    #mag = 4.9 # 2012-06-19 from Hadi
    #print('mag', mag)
    
    # now plot
    if stn != 'CDNM.HNH':
        print('rhyp', rhyp)
        vs30 = get_station_vs30(stn)[0]
        plt = makesubplt(i, fig, plt, stn, sps, mag, dep, ztor, dip, rake, rhyp, vs30, datestr, netsta)
        plt.ylim([1e-5, 0.05])
                
        i += 1
        
plt.legend(loc=3, fontsize=10)
        
plt.savefig('cmp_inter_intra_spectra.png', format='png', dpi=150, bbox_inches='tight')
plt.show()

