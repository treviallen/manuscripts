# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder):

    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn+'_HNE') >= 0:
                efile = path.join(root, filename)
            elif filename.find(stn+'_SLE') >= 0:
                efile = path.join(root, filename)

            if filename.find(stn+'_HNN') >= 0:
                nfile = path.join(root, filename)
            elif filename.find(stn+'_SLN') >= 0:
                nfile = path.join(root, filename)

    # read data
    T, SAe = read_psa(efile)
    T, SAn = read_psa(nfile)

    # get geometric mean and convert to g
    #geomean = exp((log(SAe) + log(SAn)) / 2.) / 100. # convert from %g to g

    # get max of H component
    maxsa = []
    for i in range(0,len(SAe)):
        maxsa.append(max([SAe[i], SAn[i]]) / 100.)

    return T, maxsa

def get_pkl_geomean(stndat):

    # read data
    try:
        # get geometric mean and convert to g
        geomean = exp((log(stndat['epsa']) + log(stndat['npsa'])) / 2.) # in g

    except:
        geomean = stndat['zpsa']  # if necessary

    return geomean

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):

    lines = open(psafile).readlines()[0:4]  # only use first 4 lines

    SA = []
    T = []
    for line in lines:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA)

'''
start main
'''
# start of plotting routine
from numpy import arange, array, ceil, sqrt, exp, log, logspace, interp, log10
import matplotlib.pyplot as plt
import pickle
from os import path, walk
from fnmatch import filter
from calc_oq_gmpes import crustal_gsims, interface_gsims #, get_T_index
from seyhan_stewart_2014 import seyhan_stewart_siteamp
import warnings

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
warnings.filterwarnings("ignore")

cmap = plt.cm.get_cmap('Spectral', 7)
cs = (cmap(arange(7)))
colTrue = True
ztor=5
ztor=2 # for MW 7.8 event
vs30 = 1500.
rjb = logspace(0,2.8,30)
ymax = -99999
VS30=arange(100,1500,100)  #use to correct Sa

# open fault distance file
fstn = []
fjb = []
fdist = open('cnsn_hg_rrup.txt', 'rb').readlines()[1:]
for line in fdist:
    dat = line.strip().split('\t')
    fstn.append(dat[0])
    fjb.append(float(dat[5]))

for root, dirnames, filenames in walk('HG_miniseed/2012/302/0302'):
   for filename in filter(filenames,'*_2.pkl'): 
        print filename
        pkl_file = open(path.join(root, filename),'rb')
        psa_dat = pickle.load(pkl_file)
        pkl_file.close()       
        
        plt.figure(figsize=(10, 10))
        titles = ['PGA','Sa(0.2)','Sa(0.5)','Sa(1.0)']
        
        # get event details
        mag = psa_dat[0]['Mw']
        dep = psa_dat[0]['dep']
        mag = psa_dat[0]['Mw']
        rake = psa_dat[0]['rake']
        dip = psa_dat[0]['dip']
        eqdate = psa_dat[0]['eqtime']
        repi = psa_dat[0]['dist']  # is this repi or rhyp?
        rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
          
        Tplot = [0.0, 0.2, 0.5, 1.0]
                    
        # loop thru periods
        for j, t in enumerate(Tplot):
            ax = plt.subplot(2, 2, j+1)
            Zea06c = []
            Aea14c = []
            AA13c = []
            Bea14c=[]
            Zea06i = []
            AA13i = []
            AA13cJBdist = []
            
            for ii, r in enumerate(rrup):
                jb = rjb[ii]
        
                # get ground motion estimates from GMPEs
                Zea06, CB08, CY08, BA11, Aea14, AA13, Bea14 \
                    = crustal_gsims(mag, dep, ztor, dip, rake, r, jb, vs30)
        
                if t == 0.0:
                    Zea06c.append(Zea06['pga'][0][0])                    
                    Aea14c.append(Aea14['pga'][0][0])
                    Bea14c.append(Bea14['pga'][0][0])
                    AA13c.append(AA13['pga'])
                else:
                    # interpolate log values to correct period
                    Zea06c.append(interp(t, Zea06['per'], Zea06['sa']))
                    Aea14c.append(interp(t, Aea14['per'], Aea14['sa']))
                    Bea14c.append(interp(t, Bea14['per'], Bea14['sa']))
                    AA13c.append(interp(t, AA13['per'], AA13['sa']))
                AA13cJBdist.append(AA13['eqivJB'])
                
                # get ground motion estimates from GMPEs
                Yea97imt, AB03imt, Zea06imt, AA13imt\
                    = interface_gsims(mag, dep, ztor, dip, rake, r, jb, vs30)                   
        
                if t == 0.0:                    
                    Zea06i.append(Zea06imt['pga'][0][0])
                    AA13i.append(AA13imt['pga'])

                else:
                    # interpolate log values to correct period
                    Zea06i.append(interp(t, Zea06imt['per'], Zea06imt['sa']))
                    AA13i.append(interp(t, AA13imt['per'], AA13imt['sa']))
                    
                # get max val for setting y lims
                if ii == 0:
                    ymax = exp(max([Zea06c, AA13c, Aea14c, Zea06i, AA13i]))
                   
            h0 = plt.loglog(rjb, exp(Zea06c),'-', lw=2., color=[cs[0][0],cs[0][1],cs[0][2]])
            h1 = plt.loglog(AA13cJBdist, exp(AA13c), lw=2., color=[cs[1][0],cs[1][1],cs[1][2]])           
            h2 = plt.loglog(rjb, exp(Aea14c), lw=2., color=[cs[2][0],cs[2][1],cs[2][2]])
            h3 = plt.loglog(rjb, exp(Zea06i),'-', lw=2., color=[cs[5][0],cs[5][1],cs[5][2]])
            h4 = plt.loglog(rjb, exp(AA13i),'-', lw=2., color=[cs[6][0],cs[6][1],cs[6][2]])
            h5 = plt.loglog(rjb, exp(Bea14c), lw=2., color=[cs[4][0],cs[4][1],cs[4][2]])
            # now plot recorded data
            for stndat in psa_dat: 
                Sa_corr=[]

                #uncomment if event info is in pickle file already  
                mag=stndat['Mw']
                dep=stndat['dep']
                rake=stndat['rake']
                dip=stndat['dip']
                eqdate=stndat['eqtime']
                repi=stndat['dist']
                    
                # get station Rjb for plots
                for s in range(0,len(fstn)):
                    if fstn[s] == stndat['stn']:
                        srjb = fjb[s]
                
                # get pga_r in g - assume Akkar etal 2014? model
                pga_r = exp(Aea14['pga'][0][0])
                for val in VS30: 
                    if t == 0.0:
                         try:
                             siteval = max([stndat['epga'], stndat['npga']])
                             
                         except:
                             siteval = stndat['zpga']
                             
                    else:
                         geomean = get_pkl_geomean(stndat)
                     
                         # interpolate to period of interest
                         tmpa = []
                         for g in geomean:
                             tmpa.append(g[0])
        
                         geomean = array(tmpa)
                         siteval = interp(t, stndat['periods'], geomean)
                         
                         # check max value for plotting
                         if siteval > ymax:
                             ymax = siteval
                
                # change marker and correct - soil sites ste to C/D
                    if stndat['zchan'] == 'HNZ':   
                        # correct from 360 - 760 m/s
                        siteval /= seyhan_stewart_siteamp(360., t, pga_r)
    
                        # correct from 760 - array
                        siteval *= seyhan_stewart_siteamp(val, t, pga_r)

                        h7 = plt.loglog(srjb, siteval, 'r+', markersize=10, mew=1.5)
                    else:
                        if stndat['stn'] != 'LIB':
                            # correct from 1500 - 760 m/s
                            siteval /= seyhan_stewart_siteamp(1500., t, pga_r)
        
                            # correct from 760 -array
                            siteval *= seyhan_stewart_siteamp(val, t, pga_r)   

                            h6 = plt.loglog(srjb, siteval, 'k+', markersize=10, mew=1.5)
                    
                    Sa_corr.append(siteval)
                    if t==0.2:
                        stndat['Sa_corr_0.2']=Sa_corr
                    if t==0.5:
                        stndat['Sa_corr_0.5']=Sa_corr
                    if t==1.0:
                        stndat['Sa_corr_1.0']=Sa_corr
                    
                    stndat['Vs30_corr']=VS30
            
            if j == 1:
                plt.legend((h0[0], h1[0], h2[0],h5[0], h3[0], h4[0], h6[0], h7[0]),\
                           ['Zea06crust', 'AA13wc','Aea13','Bea14crust','Zea06inter', 'AA13inter', 'Rock', 'Soil'], \
                           loc=0, fontsize=9, numpoints=1)
            
            plt.xlabel('RJB (km)', fontsize=14)
            plt.ylabel('Spectral Acceleration (g)', fontsize=14)
            plt.xlim([10, 500])
            # set ylim
            yexp = ceil(log10(ymax))
            if t == 0.2:
                yexp = 0
            plt.ylim([10**(yexp-4), 10**yexp])
            
            plt.suptitle(''.join(['Mw',str(mag),' ',str(eqdate)]),fontsize=14)
            plt.title(titles[j], fontsize=14)
            plt.grid(which='both', color='0.5')
            
            pklfile = open(path.join(root, filename),'wb')
            pickle.dump(psa_dat, pklfile, -1)
            pklfile.close()  
            
            plt.savefig(filename[:-6]+'att_vs30.pdf', format='pdf', dpi=150)
            plt.show()