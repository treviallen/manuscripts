def calc_nac_gmm_spectra(mag, rhyp, dep):
    from numpy import loadtxt, log10, log
    
    #coeffile = '//Users//trev//Documents//Earthquake_Data//2017_NAC_GMM//ncc_gmm_coeffs.BS.csv'
    coeffile = 'ncc_gmm_coeffs.BS.csv'
    
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
                
            # get Z file
            elif filename.find(stn) >= 0 and filename.find('NZ.psa') >= 0:
                zfile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('HZ.psa') >= 0:
                zfile = path.join(root, filename)
            

    try:
        try:
            # read data 
            T, SAe = read_psa(efile)
            T, SAn = read_psa(nfile)
            print(efile, nfile)
            
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
def makesubplt(i, fig, plt, sta, sps, mag, dep, ztor, dip, rake, rhyp, vs30):
    import matplotlib 
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    
    rrup = rhyp
    rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
    print('rhyp2',rhyp)
    # get ground motion estimates from GMPEs
    Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, Aea15imt, Zea16imt \
             = inslab_gsims(mag, dep, ztor, dip, rake, rhyp, rjb, vs30)
             
    Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt \
             = scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
             
    A19imt = calc_nac_gmm_spectra(mag, rhyp, dep) # use rrup
    #print(exp(A19imt['sa']))

    ax = plt.subplot(3, 3, i)
    if colTrue == 'True':
        '''
        plt.loglog(Yea97imt['per'], exp(Yea97imt['sa']), '.-' , lw=1., color=cs[0])
        plt.loglog(AB03imt['per'], exp(AB03imt['sa']),     '--', lw=1.,  color=cs[1])
        plt.loglog(Gea05imt['per'], exp(Gea05imt['sa']),   '-.', lw=1.,  color=cs[2])
        plt.loglog(Zea06imt['per'], exp(Zea06imt['sa']), '-' , lw=1.,  color=cs[3])
        plt.loglog(MP10imt['per'], exp(MP10imt['sa']), '-' , lw=2.5, color=cs[4])
        plt.loglog(Aea15imt['per'], exp(Aea15imt['sa']),     '-.', lw=2.5, color=cs[5])
        plt.loglog(Bea14imt['per'], exp(Bea14imt['sa']), '--', lw=2.5, color=cs[6])
        plt.loglog(YA15imt['per'], exp(YA15imt['sa']),   ':' , lw=2.5, color=cs[7])
        '''
        plt.loglog(Yea97imt['per'][:-1], exp(Yea97imt['sa'][:-1]), '-' , lw=1.5, color=cs[0])
        plt.loglog(AB03imt['per'][:-1], exp(AB03imt['sa'][:-1]),     '-' , lw=1.5, color=cs[1])
        plt.loglog(A12imt['per'], exp(A12imt['sa']),   '-' , lw=1.5, color=cs[2])
        plt.loglog(Zea06imt['per'], exp(Zea06imt['sa']), '-' , lw=1.5, color=cs[3])
        plt.loglog(AB06imt['per'], exp(AB06imt['sa']), '-' , lw=1.5, color=cs[4])
        plt.loglog(Aea15imt['per'], exp(Aea15imt['sa']),     '-' , lw=1.5, color=cs[5])
        plt.loglog(A19imt['per'], exp(A19imt['sa']),     '-' , lw=1.5, color=cs[6])
    else:
        plt.loglog(Yea97imt['per'], exp(Yea97imt['sa']), '-', lw=1.5,  color='0.35')
        plt.loglog(AB03imt['per'], exp(AB03imt['sa']),     '--', lw=1.5, color='0.15')
        plt.loglog(Gea05imt['per'], exp(Gea05imt['sa']),   '-.', lw=1.5, color='0.15')
        plt.loglog(CY08imt['per'], exp(CY08imt['sa']),   '-', lw=1.5,  color='0.65')
        plt.loglog(Zea06imt['per'], exp(Zea06imt['sa']), '--', lw=1.5,  color='0.65')
        plt.loglog(MP10imt['per'], exp(MP10imt['sa']), '-.', lw=1.5,  color='0.65')
        plt.loglog(Aea15imt['per'], exp(Aea15imt['sa']),     '--', lw=1.5,  color='0.35')
        
    # get recorded process_waves.py psa data
    T, geomean, pga, rhyp = get_site_geomean(sta, folder)
    plt.loglog(T, geomean, lw=1.5, color='k')

    if i >= 7:
        plt.xlabel('Period (s)', fontsize=9)
    if i == 1 or i == 4 or i == 7:
        plt.ylabel('Spectral Acceleration (g)', fontsize=9)
        
    plt.xlim([0.02, 10])
    #plt.title(' '.join((sta+'; MW =',str("%0.1f" % mag)+'; Rrup =',str(rrup),'km')), fontsize=9)
    plt.title(' '.join((sta+'; Rhyp =',str(rhyp),'km')), fontsize=9)
    plt.grid(which='both', color='0.75')

    if i == 1:
        plt.legend(['Yea97', 'AB03','A12imt','Zea06','AB06','Aea15', 'A19', 'Data'],loc=3, fontsize=7.)
        '''
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
        '''

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

from calc_oq_gmpes import inslab_gsims, scr_gsims
from gmt_tools import cpt2colormap, remove_last_cmap_colour

folder = argv[1]
prefix = argv[2]
colTrue = 'True'
#colTrue = argv[3]

# set event details
'''
mag  = 5.0
dep = 11.
'''
ztor = 10. # guess
rake = 90. # USGS CMT
dip  = 30.

# set site details
vs30 = 760

ii = 1
fig = plt.figure(ii, figsize=(10, 10))
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
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.find('psa') >= 0:
            tmpsite = filename.split('.')[3]
            tmpcomp = filename.split('.')[4]
            #sites.append('.'.join((tmpsite,tmpcomp)))
            sites.append(tmpsite)

# if two H components, rename channel
#usites = unique(sites)
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
        
usites = unique(usites)

udists = []
lolatxt = ''
for stn in usites:
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn) >= 0:
                if filename.find('NE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HE') >= 0:
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

f = open('staloc.txt', 'w')
f.write(lolatxt)
f.close()
    
    
# now sort by distance
udists = array(udists)
idx=argsort(array(udists))
udists = udists[idx]
usites = usites[idx]

# loop thru sites ordered by distance
i = 0
for stn in usites:
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn) >= 0:
                if filename.find('NE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HE') >= 0:
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
    sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(psafile)
    
    # temp fix
    #mag = 4.57 # 2012-07-20
    #mag = 5.12 # 2012-06-19
    #mag = 4.9 # 2012-06-19 from Hadi
    print('mag', mag)
    
    # now plot
    # make sub plot
    if stn != 'CDNM.HNH':
        i += 1
        print('rhyp', rhyp)
        makesubplt(i, fig, plt, stn, sps, mag, dep, ztor, dip, rake, rhyp, vs30)
        if ii == 1:
            if i <= 4:
                plt.ylim([1e-5, 1])
            else:
                plt.ylim([1e-6, 0.01])
                
        elif  rhyp > 1000.:
            plt.ylim([1e-6, 0.01])
        
        if i == 9:
            i = 0
            plt.savefig(prefix+'_'+str(ii)+'_spectra.png', format='png', dpi=150, bbox_inches='tight')
            ii += 1
            fig = plt.figure(ii, figsize=(10, 10))

plt.savefig(prefix+'_'+str(ii)+'_spectra.png', format='png', dpi=150, bbox_inches='tight')
plt.show()

