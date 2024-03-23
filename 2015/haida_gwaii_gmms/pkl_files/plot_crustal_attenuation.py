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
from numpy import arange, array, sqrt, exp, log, logspace, interp
import matplotlib.pyplot as plt
import pickle
from os import path, walk
from fnmatch import filter
from mapping_tools import distance
import datetime as dt
from calc_oq_gmpes import crustal_gsims #, get_T_index
import warnings

plt.rcParams['pdf.fonttype'] = 42
warnings.filterwarnings("ignore")

# build station list
stns = []
stnlon = []
stnlat = []
stnfile = 'canstn20140116.txt'
stnlines = open(stnfile, 'rb').readlines()
for line in stnlines:
    dat = line.strip().split('|')
    if len(dat) > 4:
        stns.append(dat[1].strip())
        stnlon.append(float(dat[4].strip()))
        stnlat.append(float(dat[3].strip()))

cmap = plt.cm.get_cmap('Spectral', 9)
cs = (cmap(arange(9)))
colTrue = True
'''
# set event details
pklfile2 = 'HG_CMT.pkl'

PGCCMT = pickle.load(open(pklfile2,'rb'))
'''
# set site details
vs30 = [760.]
ztor=5

for root, dirnames, filenames in walk('HG_miniseed/1998'):
   for filename in filter(filenames,'*_2.pkl'): 
        print filename
        psa_dat2=[]
        pkl_file = open(path.join(root, filename),'rb')
        psa_dat = pickle.load(pkl_file)
        pkl_file.close()       
        
        plt.figure(figsize=(10, 10))
        titles = ['PGA','SA(0.2)','SA(0.5)','SA(2.0)']
   
        for stndat in psa_dat:
            #comment once info is in pkl
            wftime=stndat['starttime']
            wfy=wftime.year
            wfm=wftime.month
            wfd=wftime.day
            wfh=wftime.hour
            wfmi=wftime.minute 
            wfs=wftime.second
            # make datetime object
            wft = dt.datetime(wfy, wfm, wfd, wfh, wfmi, wfs) # use 0.0 if second doesnt exist
            '''
            # get max/min times to return (+/- 5 mins)
            mindate = wft - dt.timedelta(minutes=5)
            maxdate = wft + dt.timedelta(minutes=5)
            
            for event in PGCCMT:
                if event['Mw']!='':
                   # set PCG datetime
                    cmttime = dt.datetime(int(event['Year']),int(event['Month']),\
                    int(event['Day']),int(event['Hour']),int(event['Min']),\
                    int(event['Sec']))
                    if cmttime >= mindate and cmttime <= maxdate:
                        mag=event['Mw']
                        dep=event['Depth']
                        rake=event['Rake1']
                        dip=event['Dip1']
                        eqlat=event['Latitude']
                        eqlon=event['Longitude']
                        mech=event['Mech']
                        eqdate=dt.datetime(int(event['Year']),int(event['Month']),\
                        int(event['Day']),int(event['Hour']),int(event['Min']),\
                        int(event['Sec']))
            '''
            #uncomment if event info is in spickle file already  
            
            print stndat['stn']
            mag=stndat['Mw']
            dep=stndat['dep']
            ztor=5
            rake=stndat['rake']
            dip=stndat['dip']
            eqlat=stndat['eqlat']
            eqlon=stndat['eqlon']
            mech=stndat['mech']
            eqdate=stndat['eqtime']
            

            rjb = logspace(0,2.6,50)
            rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
            
            Tplot = [0.0, 0.2, 0.5, 2.0]
           
            acc = []
            stnrjb = []           
                        
            # loop thru periods
            for j, t in enumerate(Tplot):
                ax = plt.subplot(2, 2, j+1)
                Zea06r = []
                CB08r = []
                CY08r = []
                BA11r = []
                Aea13r = []
                AA13r = []
                
                for ii, r in enumerate(rrup):
                    jb = rjb[ii]
            
                    # get ground motion estimates from GMPEs
                    Zea06imt, CB08imt, CY08imt, BA11imt, Aea13imt, AA13imt \
                        = crustal_gsims(mag, dep, ztor, dip, rake, r, jb, vs30)
            
                    if t == 0.0:
                        Zea06r.append(Zea06imt['pga'][0])
                        CB08r.append(CB08imt['pga'][0])
                        CY08r.append(CY08imt['pga'][0])
                        BA11r.append(BA11imt['pga'][0])
                        Aea13r.append(Aea13imt['pga'][0])
                        AA13r.append(AA13imt['pga'])
                    else:
                        # interpolate log values to correct period
                        Zea06r.append(interp(t, Zea06imt['per'], Zea06imt['sa']))
                        CB08r.append(interp(t, CB08imt['per'], CB08imt['sa']))
                        CY08r.append(interp(t, CY08imt['per'], CY08imt['sa']))
                        BA11r.append(interp(t, BA11imt['per'], BA11imt['sa']))
                        Aea13r.append(interp(t, Aea13imt['per'], Aea13imt['sa']))
                        AA13r.append(interp(t, AA13imt['per'][::-1], AA13imt['sa'][::-1]))
                
                if colTrue == True:
                    plt.loglog(rjb, exp(Zea06r), lw=2., color=[cs[0][0],cs[0][1],cs[0][2]])
                    plt.loglog(rjb, exp(CB08r),  lw=2., color=[cs[1][0],cs[1][1],cs[1][2]])
                    plt.loglog(rjb, exp(CY08r),  lw=2., color=[cs[2][0],cs[2][1],cs[2][2]])
                    plt.loglog(rjb, exp(BA11r),  lw=2., color=[cs[3][0],cs[3][1],cs[3][2]])
                    plt.loglog(rjb, exp(Aea13r), lw=2., color=[cs[5][0],cs[5][1],cs[5][2]])
                    plt.loglog(rjb, exp(AA13r), lw=2., color=[cs[6][0],cs[6][1],cs[6][2]])
                else:
                    plt.loglog(rjb, exp(Zea06r), '-', lw=2.,  color='0.15')
                    plt.loglog(rjb, exp(CB08r),  '--', lw=2., color='0.15')
                    plt.loglog(rjb, exp(CY08r),  '-.', lw=2., color='0.15')
                    plt.loglog(rjb, exp(BA11r),  '-', lw=2.,  color='0.55')
                    plt.loglog(rjb, exp(Aea13r), '.', lw=2.,  color='0.55')
                    plt.loglog(rjb, exp(AA13r), '.', lw=2.,  color='0.55')
                
                # get response spectra for given period and plot
                if t == 0.0:
                    try:
                        accval=max([stndat['epga'], stndat['npga']])
                    except:
                        accval=stndat['zpga']
                else:
                    geomean = get_pkl_geomean(stndat)

                    # interpolate to period value
                    accval=interp(t, stndat['periods'], geomean[:,0])
                
                for s, stn in enumerate(stns):
                    if stn == stndat['stn']: 
                        repi, az, baz = distance(eqlat, eqlon, stnlat[s], stnlon[s])
                        # assume rjb = repi and eq is pt source                        
                 
                #else:
                 #   print 'Cannot find data for site:', stn
                
                # now plot             
                plt.loglog(repi, accval, 'k+', markersize=1.5*mag, markeredgewidth=1.75)
                plt.xlabel('Rjb (km)')
                plt.ylabel('Spectral Acceleration (g)')
                plt.xlim([3, 400])
                plt.ylim([1E-4, 1])
                plt.suptitle(''.join(['Mw',str(mag),' ',str(eqdate)]),fontsize=14)
                plt.title(titles[j])
                plt.grid(which='both', color='0.5')
                #plt.close()
                
                if j == 0:
                    plt.legend(['Zea06', 'CB08','CY08','BA11','Aea14','AA13'],loc=3)
                    leg = plt.gca().get_legend()
                    ltext  = leg.get_texts()
                    plt.setp(ltext, fontsize='xx-small')
                #save data into pickle file
                stnrjb.append(repi) 
                acc.append(accval)
            '''      
            stndat['Mw']=mag
            stndat['eqlat']=eqlat
            stndat['eqlon']=eqlon
            stndat['mech']=mech
            stndat['azimuth']=az
            stndat['dep']=dep
            stndat['dip']=dip
            stndat['rake']=rake
            stndat['baz']=baz
            stndat['dist']=list(set(stnrjb))
            stndat['accval']=acc
            stndat['Tplot']=Tplot
            stndat['eqtime']=eqdate
       
            psa_dat2.append(stndat)
            
        # save to picklefile2       pklfile = open(path.join(root, filename.replace('_2.','_3.')), 'wb')
        pklfile = open(path.join(root, filename), 'wb')        
        pickle.dump(psa_dat2, pklfile, -1)
        pklfile.close()
        '''       
        #plt.show()        
        plt.savefig(path.join(root,''.join(filename[:-6]+'attenuation_.pdf')), format='pdf', dpi=150)
             