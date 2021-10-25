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
    
    # get near source correction
    if log10(rhyp) < hx[0]:
        near_field_term = n0 * (log10(rhyp)-hx[0])
    else:
        near_field_term = 0.0
    
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
def makesubplt(i, fig, plt, sta, sps, mag, dep, ztor, dip, rake, rhyp, vs30, datestr):
    import matplotlib 
    matplotlib.rc('xtick', labelsize=12) 
    matplotlib.rc('ytick', labelsize=12)
    props = dict(boxstyle='square', facecolor='w', alpha=1)
    
    #dates = ['1995-12-25', '1995-12-25', '2019-06-24', '2019-06-24']
    
    rrup = rhyp
    rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
    
    A19imt_BS = calc_nac_gmm_spectra(mag, rhyp, dep, vs30, 'BS') # use rrup
    
    nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup, vs30)
    
    Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06imt, MP10imt, Aea16imt, Zea16imt, Kea20imt \
            = inslab_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
        
    # adjust NGA-E from 3000 -> target
    nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
    
    ax = plt.subplot(2, 2, i)
    if colTrue == 'True':
        plt.loglog(nga_e_imt['per'], exp(nga_e_imt['sa']),'--' , lw=1.5, color=cs[0])
        plt.loglog(Kea20imt['per'], exp(Kea20imt['sa']),'-.' , lw=2., color=cs[1])
        plt.loglog(A19imt_BS['per'], exp(A19imt_BS['sa']),'-' , lw=1.5, color=cs[2])
        
        
    # get recorded process_waves.py psa data
    T, geomean, pga, rhyp = get_site_geomean(sta, folder)
    plt.loglog(T, geomean, lw=1.5, color='k')

    if i >= 3:
        plt.xlabel('Period (s)', fontsize=14)
    if i == 1 or i == 3 or i == 7:
        plt.ylabel('Spectral Acceleration (g)', fontsize=14)
        
    plt.xlim([0.02, 10])
    plt.ylim([1e-4, .1])
    plt.grid(which='both', color='0.75')

    if i == 1:
        #plt.legend(['Yea97', 'AB06','A12imt','Aea16', 'A19 (BS)', 'A19 (NGH)', 'A19 (OB)','Data'],loc=3, fontsize=7.)
        plt.legend(['Goulet et al (2017)', 'Kuehn et al (2020)', 'Present Study', 'Geometric Mean'],loc=4, fontsize=10)
        
    txtbox = datestr+'\nStation: AU.'+sta+'\n'+'$\mathregular{M_W}$: '+str(mag) \
             +'\n$\mathregular{R_{hyp}}$: '+str(int(round(rhyp)))+' km'\
             +'\n$\mathregular{h_z}$: '+str(int(round(dep)))+' km'\
             +'\n'+'$\mathregular{V_{S30}}$: '+ str(int(round(vs30))) +' m/s'
    print(txtbox)
    
    xtxt = get_log_xy_locs(ax.get_xlim(), 0.03)
    ytxt = get_log_xy_locs(ax.get_ylim(), 0.03)
    
    lspace = (2, 2, 1, 1)
    plt.text(xtxt, ytxt, txtbox, size=11, ha='left', va='bottom', weight='normal') #, linespacing=lspace) #, bbox=props)
    
    # add letter
    ylims = ax.get_ylim()
    plt.text(0.017, ylims[1]*1.25, letters[i-1], va='bottom', ha ='right', fontsize=17)

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
from calc_oq_gmpes import scr_gsims, get_station_vs30, nga_east_mean, adjust_gmm_with_nga_east, inslab_gsims
from gmt_tools import cpt2colormap, remove_last_cmap_colour

folder = 'hsd_plot_spectra'
prefix = 'hsd_spectra'
colTrue = 'True'
#colTrue = argv[3]

# set event details

ii = 1
fig = plt.figure(ii, figsize=(10, 10))
#cmap = plt.cm.get_cmap('Spectral', 7)
ncols = 7
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'
#cptfile = 'U:\\DATA\\GMT\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))
cs = get_mpl2_colourlist()

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

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
    ztor = 100 # dummy
    dip = 30 # dummy
    rake = -90
    
    # now plot
    if stn != 'CDNM.HNH':
        i += 1
        print('rhyp', rhyp)
        vs30 = get_station_vs30(stn)[2] # use USGS
        if isnan(vs30):
            vs30 = 760
        
        datestr = psafile.split('/')[-1].split('T')[0]    
        makesubplt(i, fig, plt, stn, sps, mag, dep, ztor, dip, rake, rhyp, vs30, datestr)
        if ii == 1:
            if rhyp <= 900:
                plt.ylim([1e-4, .1])
            else:
                plt.ylim([2e-5, 0.01])
                
        elif  rhyp > 1000.:
            plt.ylim([2e-5, 0.01])
            
        if i == 9:
            i = 0
            plt.savefig('event_psa/'+prefix+'_'+str(ii)+'_spectra.png', format='png', dpi=300, bbox_inches='tight')
            ii += 1
            fig = plt.figure(ii, figsize=(10, 10))

plt.savefig('figures/'+prefix+'_spectra.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig_12.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

