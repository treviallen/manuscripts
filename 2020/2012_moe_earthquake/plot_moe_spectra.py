
# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder):
    print(stn)
    from fnmatch import filter
    from os import path, walk, system
    from numpy import sqrt
    
    zcomp = False

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
        
        zcomp = True
        
    return T, geomean, pga, rhyp, zcomp

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array
    
    
    lines = open(psafile).readlines() 
    if lines[22].startswith('---'):
        stidx = 23
    else:
        stidx = 24

    SA = []
    T = []
    for line in lines[stidx:]:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA)
    
# get other details for plotting
def read_psa_details(psafile):
    lines = open(psafile).readlines()
    
    if lines[22].startswith('---'):
        enidx = 23
    else:
        enidx = 24

    for line in lines[:enidx]:
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
    
    rrup = rhyp # subtract for misc distance to rup
    rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
    
    # get ground motion estimates from GMPEs
    Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt , SP16imt \
             = scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
    
    Tea19 = tang2019_cam_gsim(mag, dep, rrup, vs30)
    
    # adjust some GMMs using Seyhan & Stewart 2014
    A12imt = adjust_gmm_with_SS14(A12imt, 820., vs30)
    Sea09imt = adjust_gmm_with_SS14(Sea09imt, 865., vs30)
    
    # get mean USGS NGA-E
    nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup, vs30)
    #plt.loglog(nga_e_imt['per'], exp(nga_e_imt['sa']),'--' , lw=1.5, color=cs[5]) # for testing amp factors
    nga_e_imt_hold = nga_e_imt
    
    # adjust NGA-E from 3000 -> target
    nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
    
    #print(nga_e_imt['per'])
    #print(exp(nga_e_imt['sa']) / exp(nga_e_imt_hold['sa']))
    
    ax = plt.subplot(3, 3, i)
    if colTrue == 'True':
        plt.loglog(AB06imt['per'], exp(AB06imt['sa']), '-' , lw=1.5, color=cs[0])
        plt.loglog(Sea09imt['per'], exp(Sea09imt['sa']), '-' , lw=1.5, color=cs[1])
        plt.loglog(A12imt['per'], exp(A12imt['sa']),'-' , lw=1.5, color=cs[2])
        plt.loglog(Bea14imt['per'], exp(Bea14imt['sa']),'-' , lw=1.5, color=cs[3])
        #plt.loglog(YA15imt['per'], exp(YA15imt['sa']),'-' , lw=1.5, color=cs[4])
        plt.loglog(nga_e_imt['per'], exp(nga_e_imt['sa']),'-' , lw=1.5, color=cs[5])
        plt.loglog(Tea19['per'], exp(Tea19['sa']),'-' , lw=1.5, color=cs[6])
        
        
    # get recorded process_waves.py psa data
    T, geomean, pga, rhyp, zcomp = get_site_geomean(sta, folder)
    plt.loglog(T, geomean, lw=1.5, color='k')

    if i >= 7:
        plt.xlabel('Period (s)', fontsize=9)
    if i == 1 or i == 4 or i == 7:
        plt.ylabel('Spectral Acceleration (g)', fontsize=9)
    
    plt.xlim([0.01, 10])
    #plt.title(' '.join((sta+'; MW =',str("%0.1f" % mag)+'; Rrup =',str(rrup),'km')), fontsize=9)
    if zcomp == True:
        sta += '*'
    plt.title(' '.join((sta+'; Rhyp =',str(rhyp),'km; VS30 =', str(int(round(vs30))))), fontsize=8)
    plt.grid(which='both', color='0.75')
    
    if i == 1:
        #plt.legend(['Yea97', 'AB06','A12imt','Aea16', 'A19 (BS)', 'A19 (NGH)', 'A19 (OB)','Data'],loc=3, fontsize=7.)
        plt.legend(['AB06','Sea09', 'A12','Bea14', 'NGA-E', 'Tea19', 'Data'],loc=3, fontsize=9.)
        
    return ax

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort, isnan
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

from calc_oq_gmpes import inslab_gsims, scr_gsims, tang2019_cam_gsim, \
                          nga_east_mean, get_station_vs30, adjust_gmm_with_SS14, \
                          adjust_gmm_with_nga_east
from data_fmt_tools import return_sta_data
from mapping_tools import distance
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from misc_tools import get_mpl2_colourlist

folder = argv[1]
prefix = argv[2]
colTrue = 'True'
#colTrue = argv[3]

# set event details
if prefix.startswith('201206'):
    mag  = 5.16
    eqdep = 18.0
    eqlat = -38.259
    eqlon = 146.290

elif prefix.startswith('201207'):
    mag  = 4.40
    eqdep = 12.41
    eqlat = -38.231
    eqlon = 146.215

ztor = 8. # guess
rake = 90. # USGS CMT
dip  = 30.

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

ii = 1
fig = plt.figure(ii, figsize=(10, 10))
#cmap = plt.cm.get_cmap('Spectral', 7)

ncols = 8
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
#cptfile = 'U:\\DATA\\GMT\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

#cs = get_mpl2_colourlist()

from fnmatch import filter
from os import path, walk, system

# build sites
sites = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.find('psa') >= 0:
            tmpsite = filename.split('.')[1]
            tmpcomp = filename.split('.')[2]
            #sites.append('.'.join((tmpsite,tmpcomp)))
            sites.append(tmpsite)

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
usites = unique(sites)

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
    sta, sps, rhyp, pga, pgv, mag, dep, stlo, stla = read_psa_details(psafile)
    
    # temp fix
    # set event details
    if prefix.startswith('201206'):
        mag  = 5.0
        #dep = 11.
    elif prefix.startswith('201207'):
        mag  = 4.4
        #dep = 11.

    
    # now plot
    if stn != 'CDNMX.HNH':
        i += 1
        print('rhyp', rhyp)
        vs30, isproxy, usgsvs, asscmvs, kvs, stla, stlo = get_station_vs30(stn)
        #sta_dat = return_sta_data(stn)
        
        # now calculate distance on the fly
        repi = distance(eqlat, eqlon, stla, stlo)[0]
        rhyp = sqrt(repi**2 + eqdep**2)
        
        if isnan(vs30):
            vs30 = 760.
        
        ax = makesubplt(i, fig, plt, stn, sps, mag, eqdep, ztor, dip, rake, rhyp, vs30)
        if ii == 1:
            if prefix.startswith('201206'):
                if rhyp <= 100:
                    plt.ylim([1e-4, .2])
                else:
                    plt.ylim([1e-5, 0.02])
                    
            elif prefix.startswith('201207'):
                if rhyp <= 20:
                    plt.ylim([1e-4, .3])
                elif rhyp <= 50:
                    plt.ylim([4e-5, .1])
                else:
                    plt.ylim([1e-5, 0.03])
                
        elif rhyp > 1000.:
            plt.ylim([2e-5, 0.01])
        
        # add letter
        ylims = ax.get_ylim()
        plt.text(0.008, ylims[1]*1.1, letters[i-1], va='bottom', ha ='right', fontsize=13)
        
        if i == 9:
            i = 0
            plt.savefig(prefix+'_'+str(ii)+'_spectra.png', format='png', dpi=150, bbox_inches='tight')
            ii += 1
            fig = plt.figure(ii, figsize=(10, 10))

plt.savefig(prefix+'_'+str(ii)+'_spectra.png', format='png', dpi=150, bbox_inches='tight')
plt.show()

