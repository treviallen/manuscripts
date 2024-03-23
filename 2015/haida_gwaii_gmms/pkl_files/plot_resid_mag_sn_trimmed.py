def get_pkl_geomean(stndat):

    # read data
    try:
        # get geometric mean and convert to g
        geomean = exp((log(stndat['epsa']) + log(stndat['npsa'])) / 2.) # in g
    except:
        geomean = stndat['zpsa']  # if necessary
    return geomean

'''
start main
'''
import matplotlib.pyplot as plt
import pickle
from mapping_tools import distance
from numpy import array, arange, sqrt, exp, log, subtract, add, std
from calc_oq_gmpes import crustal_gsims
import warnings
import datetime as dt
from os import path, walk
from fnmatch import filter
from scipy import interp

warnings.filterwarnings("ignore")

plt.rcParams['pdf.fonttype'] = 42
CB08imt_rsum_a=[0]*100
Zea06imt_rsum_a=[0]*100
CY08imt_rsum_a=[0]*100
BA11imt_rsum_a=[0]*100
Aea14imt_rsum_a=[0]*100  ## _a M4-6    
AA13imt_rsum_a=[0]*100
count_a=0             

CB08imt_rsum_b=[0]*100
Zea06imt_rsum_b=[0]*100
CY08imt_rsum_b=[0]*100
BA11imt_rsum_b=[0]*100
Aea14imt_rsum_b=[0]*100  ## _b M5-6   
AA13imt_rsum_b=[0]*100
count_b=0

CB08imt_rsum_c=[0]*100
Zea06imt_rsum_c=[0]*100
CY08imt_rsum_c=[0]*100
BA11imt_rsum_c=[0]*100
Aea14imt_rsum_c=[0]*100  ## _c M6+    
AA13imt_rsum_c=[0]*100
count_c=0

'''
load sig/noise data
'''
pkl_file = '..\\stn_noise_models.pkl'
sn_dat = pickle.load(pkl_file)
   
for root, dirnames, filenames in walk('HG_miniseed'):
   for filename in filter(filenames,'*_2.pkl'):   
        pkl_file = open(path.join(root, filename),'rb')
        psadat = pickle.load(pkl_file)
        pkl_file.close()
        
        # set site details
        vs30 = [760.]
        ztor=5
        # set event details
        for i, stndat in enumerate(psadat): 
            eqtime=stndat['eqtime']
            az=stndat['azimuth']
            baz=stndat['baz']
            repi=stndat['dist'][0]
            mag=stndat['Mw']
            dep=stndat['dep']
            dip=stndat['dip']            
            rake=stndat['rake']  

            # assume rjb = repi and eq is pt source       
            rjb  = repi
            rrup = sqrt(repi**2 + dep**2)
            
            geomean = get_pkl_geomean(stndat)
            geomean=geomean[:,0]
            
            print stndat['stn'], mag, rrup
            if mag >5 and mag<=6 and rrup <100 and stndat['stn']!= 'LIB':
                count_b=count_b+1
                
                Zea06imt, CB08imt, CY08imt, BA11imt, Aea14imt, AA13imt \
                = crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
                
                CB08imt_b=interp(stndat['periods'],CB08imt['per'],CB08imt['sa'])
                CB08imt_rsum_b=add(CB08imt_rsum_b,(subtract(log(geomean),log(exp(CB08imt_b)))))
                CB08imt_rmean_b=CB08imt_rsum_b/count_b
                CB08imt_std_b=std(CB08imt_rmean_b)
                
                CY08imt_b=interp(stndat['periods'],CY08imt['per'],CY08imt['sa'])
                CY08imt_rsum_b=add(CY08imt_rsum_b,(subtract(log(geomean),log(exp(CY08imt_b)))))
                CY08imt_rmean_b=CY08imt_rsum_b/count_b
                CY08imt_std_b=std(CY08imt_rmean_b)
                
                Zea06imt_b=interp(stndat['periods'],Zea06imt['per'],Zea06imt['sa'])
                Zea06imt_rsum_b=add(Zea06imt_rsum_b,(subtract(log(geomean),log(exp(Zea06imt_b)))))
                Zea06imt_rmean_b=Zea06imt_rsum_b/count_b
                Zea06imt_std_b=std(Zea06imt_rmean_b)
                
                BA11imt_b=interp(stndat['periods'],BA11imt['per'],BA11imt['sa'])
                BA11imt_rsum_b=add(BA11imt_rsum_b,(subtract(log(geomean),log(exp(BA11imt_b)))))
                BA11imt_rmean_b=BA11imt_rsum_b/count_b
                BA11imt_std_b=std(BA11imt_rmean_b)
                
                Aea14imt_b=interp(stndat['periods'],Aea14imt['per'],Aea14imt['sa'])
                Aea14imt_rsum_b=add(Aea14imt_rsum_b,(subtract(log(geomean),log(exp(Aea14imt_b)))))                        
                Aea14imt_rmean_b=Aea14imt_rsum_b/count_b
                Aea14imt_std_b=std(Aea14imt_rmean_b)
                
                AA13imt_b=interp(stndat['periods'], AA13imt['per'][::-1], AA13imt['sa'][::-1])
                AA13imt_rsum_b=add(AA13imt_rsum_b,(subtract(log(geomean),log(exp(AA13imt_b)))))                        
                AA13imt_rmean_b=AA13imt_rsum_b/count_b
                AA13imt_std_b=std(AA13imt_rmean_b)
            else:
                if mag<=5 and rrup <100 and stndat['stn']!= 'LIB':
                    count_a=count_a+1
                    
                    Zea06imt, CB08imt, CY08imt, BA11imt, Aea14imt, AA13imt \
                    = crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
                
                    CB08imt_a=interp(stndat['periods'],CB08imt['per'],CB08imt['sa'])
                    CB08imt_rsum_a=add(CB08imt_rsum_a,(subtract(log(geomean),log(exp(CB08imt_a)))))
                    CB08imt_rmean_a=CB08imt_rsum_a/count_a
                    CB08imt_std_a=std(CB08imt_rmean_a)
                    
                    CY08imt_a=interp(stndat['periods'],CY08imt['per'],CY08imt['sa'])
                    CY08imt_rsum_a=add(CY08imt_rsum_a,(subtract(log(geomean),log(exp(CY08imt_a)))))
                    CY08imt_rmean_a=CY08imt_rsum_a/count_a
                    CY08imt_std_a=std(CY08imt_rmean_a)
                    
                    Zea06imt_a=interp(stndat['periods'],Zea06imt['per'],Zea06imt['sa'])
                    Zea06imt_rsum_a=add(Zea06imt_rsum_a,(subtract(log(geomean),log(exp(Zea06imt_a)))))
                    Zea06imt_rmean_a=Zea06imt_rsum_a/count_a
                    Zea06imt_std_a=std(Zea06imt_rmean_a)
                    
                    BA11imt_a=interp(stndat['periods'],BA11imt['per'],BA11imt['sa'])
                    BA11imt_rsum_a=add(BA11imt_rsum_a,(subtract(log(geomean),log(exp(BA11imt_a)))))
                    BA11imt_rmean_a=BA11imt_rsum_a/count_a
                    BA11imt_std_a=std(BA11imt_rmean_a)
                    
                    Aea14imt_a=interp(stndat['periods'],Aea14imt['per'],Aea14imt['sa'])
                    Aea14imt_rsum_a=add(Aea14imt_rsum_a,(subtract(log(geomean),log(exp(Aea14imt_a)))))                        
                    Aea14imt_rmean_a=Aea14imt_rsum_a/count_a
                    Aea14imt_std_a=std(Aea14imt_rmean_a)   
                
                    AA13imt_a=interp(stndat['periods'], AA13imt['per'][::-1], AA13imt['sa'][::-1])
                    AA13imt_rsum_a=add(AA13imt_rsum_a,(subtract(log(geomean),log(exp(AA13imt_a)))))                        
                    AA13imt_rmean_a=AA13imt_rsum_a/count_a
                    AA13imt_std_a=std(AA13imt_rmean_a)            
            
                else:
                    if mag>6 and rrup <100 and stndat['stn']!= 'LIB':
                        count_c=count_c+1
                        
                        Zea06imt, CB08imt, CY08imt, BA11imt, Aea14imt, AA13imt \
                        = crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
                    
                        CB08imt_c=interp(stndat['periods'],CB08imt['per'],CB08imt['sa'])
                        CB08imt_rsum_c=add(CB08imt_rsum_c,(subtract(log(geomean),log(exp(CB08imt_c)))))
                        CB08imt_rmean_c=CB08imt_rsum_c/count_c
                        CB08imt_std_c=std(CB08imt_rmean_c)
                        
                        CY08imt_c=interp(stndat['periods'],CY08imt['per'],CY08imt['sa'])
                        CY08imt_rsum_c=add(CY08imt_rsum_c,(subtract(log(geomean),log(exp(CY08imt_c)))))
                        CY08imt_rmean_c=CY08imt_rsum_c/count_c
                        CY08imt_std_c=std(CY08imt_rmean_c)
                        
                        Zea06imt_c=interp(stndat['periods'],Zea06imt['per'],Zea06imt['sa'])
                        Zea06imt_rsum_c=add(Zea06imt_rsum_c,(subtract(log(geomean),log(exp(Zea06imt_c)))))
                        Zea06imt_rmean_c=Zea06imt_rsum_c/count_c
                        Zea06imt_std_c=std(Zea06imt_rmean_c)
                        
                        BA11imt_c=interp(stndat['periods'],BA11imt['per'],BA11imt['sa'])
                        BA11imt_rsum_c=add(BA11imt_rsum_c,(subtract(log(geomean),log(exp(BA11imt_c)))))
                        BA11imt_rmean_c=BA11imt_rsum_c/count_c
                        BA11imt_std_c=std(BA11imt_rmean_c)
                        
                        Aea14imt_c=interp(stndat['periods'],Aea14imt['per'],Aea14imt['sa'])
                        Aea14imt_rsum_c=add(Aea14imt_rsum_c,(subtract(log(geomean),log(exp(Aea14imt_c)))))                        
                        Aea14imt_rmean_c=Aea14imt_rsum_c/count_c
                        Aea14imt_std_c=std(Aea14imt_rmean_c) 
                        
                        AA13imt_c=interp(stndat['periods'], AA13imt['per'][::-1], AA13imt['sa'][::-1])
                        AA13imt_rsum_c=add(AA13imt_rsum_c,(subtract(log(geomean),log(exp(AA13imt_c)))))                        
                        AA13imt_rmean_c=AA13imt_rsum_c/count_c
                        AA13imt_std_c=std(AA13imt_rmean_c)

fig = plt.figure(figsize=(10, 10)) 
for i, model in enumerate(['Zea06imt', 'CB08imt', 'CY08imt', 'BA11imt', 'Aea14imt', 'AA13imt']):
    rmean_b=eval(''.join(model+'_rmean_b'))
    rstd_b=eval(''.join(model+'_std_b'))
    rmean_a=eval(''.join(model+'_rmean_a'))
    rstd_a=eval(''.join(model+'_std_a'))
    rmean_c=eval(''.join(model+'_rmean_c'))
    rstd_c=eval(''.join(model+'_std_c'))
    ax=plt.subplot(3,2,i)
    plt.grid(which='none', color='0.5')
    plt.xlim([0.025, 5])
    plt.ylim([-3.0, 3])
    plt.xlabel('Period (s)',fontsize=8)
    plt.ylabel('residuals', fontsize=8) 
    plt.title(model,fontsize=10)
    plt.setp(ax.get_xticklabels(),fontsize=8)
    plt.setp(ax.get_yticklabels(),fontsize=8)
    plt.semilogx(stndat['periods'],[0]*len(stndat['periods']),'k--')
    
    plt.semilogx(stndat['periods'],rmean_a,'r',label=''.join('Mw<5 '+str([count_a])))
    plt.semilogx(stndat['periods'],subtract(rmean_a,rstd_a),'r:',stndat['periods'],\
    add(rmean_a,rstd_a),'r:')    
    
    plt.semilogx(stndat['periods'],rmean_b,'b',label=''.join('5<Mw<6 '+str([count_b])))
    plt.semilogx(stndat['periods'],subtract(rmean_b,rstd_b),'b:',stndat['periods'],\
    add(rmean_b,rstd_b),'b:')
    
    plt.semilogx(stndat['periods'],rmean_c,'g',label=''.join('Mw>6 '+str([count_c])))
    plt.semilogx(stndat['periods'],subtract(rmean_c,rstd_c),'g:',stndat['periods'],\
    add(rmean_c,rstd_c),'g:')    
    
    plt.suptitle('log(Observed)-log(Predicted) $\pm\\sigma$ (<100km)',size=12)
    plt.subplots_adjust(top=0.9,hspace=0.4)
                        
    if i ==1: 
        plt.legend(loc=0,fontsize=8)
                      
plt.savefig('HG_miniseed/plots/Mw_res.pdf', format='pdf', dpi=150)                    
plt.show()
#plt.close()