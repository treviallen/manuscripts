#import pickle

#wadat = pickle.load(open('wa_recs.pkl','rb'))
#bkndat = pickle.load(open('wa_eqSource.pkl','rb'))

from obspy import read, Trace, Stream, UTCDateTime
from obspy.taup import TauPyModel, taup_time
from mapping_tools import distance, km2deg, reckon
from response import get_response_info, paz_response, deconvolve_instrument
from spectral_analysis import calc_fft, prep_psa_simple, calc_response_spectra, get_cor_velocity
from misc_tools import listdir_extension
from data_fmt_tools import return_sta_data
from datetime import datetime, timedelta
from os import path, getcwd, remove
from numpy import asarray, arange, array, log10, mean, percentile, where, isnan, interp
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
mpl.style.use('classic')


from data_fmt_tools import get_stn_dataless_seed, get_station_distance, \
                           remove_low_sample_data


###############################################################################
# set params
###############################################################################

tstep = 5 # sec
maptimes = arange(-60, 12*60+1, tstep)
freqmin = 0.01 # Hz
freqmax = 0.02


##########################################################################################
# parse eq epicentres
##########################################################################################

def parse_usgs_events(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
    lines = open(usgscsv).readlines()[1:]
    
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'datetime': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict

usgscsv = '20180225_png.csv'
# parse catalogue
ev = parse_usgs_events(usgscsv)[0]

# set event eqlo = ev['lon']
eqla = ev['lat']
eqmag = ev['mag']
eqdp = ev['dep']
eqdt = ev['datetime']

###############################################################################
# get travel times
###############################################################################
print('Getting travel times...')
# set arrivals
model = TauPyModel(model="iasp91")
kmrng = arange(0, 3500, 100)  
ptimes = []
stimes = []
for kmr in kmrng:
    arrivals = model.get_travel_times(source_depth_in_km=ev['dep'], distance_in_degree=km2deg(kmr), phase_list=['P', 'S'])
    #print(arrivals)
    
    gotp = False
    gots = False
    for a in arrivals:
        if a.name.upper() == 'P' and gotp == False:
            ptimes.append(a.time)
            gotp = True
        if a.name.upper() == 'S' and gots == False:
            stimes.append(a.time)
            gots = True
            
    if len(arrivals) == 0:
        ptimes.append(0.)
        stimes.append(0.)

ptimes = array(ptimes)
stimes = array(stimes)

###############################################################################
# parse mseed
###############################################################################
"""
mseedfiles = listdir_extension('mseed', 'mseed')
#mseedfiles = listdir_extension('ausarray_mseed', 'mseed')

# loop through mseed files
sta_dict = []
for mseedfile in mseedfiles:
    print(mseedfile)
    # read mseed
    st = read(path.join('mseed', mseedfile))
    #st = read(path.join('ausarray_mseed', mseedfile))
    
    # remove junk channels
    st = remove_low_sample_data(st)
    
    # merge
    st = st.merge(method=0, fill_value='interpolate')
    
    # loop through traces and get Z component
    for tr in st:
        if tr.stats.channel.endswith('Z'):
            # filter
            tr_new = tr.resample(1.0)
            tr_new = tr_new.detrend(type='constant')
            tr_new = tr_new.taper(0.05, type='hann', max_length=None, side='both')
            
            nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
                  = get_response_info(tr.stats.station, eqdt, tr.stats.channel)
            # eqdt.datetime
            
            #print(pazfile, tr.stats.station)
            # get fft of trace
            freq, wavfft = calc_fft(tr_new.data, tr.stats.sampling_rate)
            
            # get response for given frequencies
            real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                                gain, inst_ty)
            # deconvolve response
            corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)
            
            # make new instrument corrected velocity trace
            pgv, ivel = get_cor_velocity(corfftr, corffti, freq, inst_ty)
            
            tr_new.data = ivel.real
            #print(tr.data)
            
            tr_new = tr_new.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
            # correct response
            
            times = tr_new.times('utcdatetime')
            sampled_vel = []
            # loop through time steps and get peak velocity
            for mt in maptimes:
                sample_time = UTCDateTime(ev['datetime']) + mt
                
                index = times.searchsorted(sample_time)
                try:
                    sampled_vel.append(tr_new.data[index])
                except:
                    sampled_vel.append(0)
    
    sta_data = return_sta_data(st[0].stats.station)
    
    sta_dat = {'sampled_vel':array(sampled_vel), 'sta':st[0].stats.station, 'stlo': sta_data['stlo'], 'stla': sta_data['stla']}
    
    sta_dict.append(sta_dat)
    
    #arrivals = model.get_travel_times(source_depth_in_km=arrival_dep, distance_in_degree=rngdeg)

print('Saving pkl file...')
pklfile = open("sta_dict.pkl", "wb" )
pickle.dump(sta_dict, pklfile) #, protocol=-1)
pklfile.close()
"""
print('Loading pkl file...')
sta_dict = pickle.load(open("sta_dict.pkl", "rb" ))

###############################################################################
# set map
###############################################################################

urcrnrlat = 1.0
llcrnrlat = -28.
urcrnrlon = 152.
llcrnrlon = 126
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

#plt.tick_params(labelsize=16)


# make a map for each time step
props = dict(boxstyle='round', facecolor='w', alpha=1)
azimuths = arange(0, 360, 1)
#norm = mpl.colors.Normalize(vmin=-6, vmax=6)
norm = mpl.colors.Normalize(vmin=-1E-6, vmax=1E-6)
for i, mt in enumerate(maptimes):
    
    fig = plt.figure(1, figsize=(18,10))
    ax = fig.add_subplot(111)

    # initialise map
    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='c',area_thresh=1000.)

    # draw coastlines, state and country boundaries, edge of map.
    #m.shadedrelief()
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.fillcontinents(color='w',lake_color='0.9')
    m.drawmapboundary(fill_color='0.9')
    m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

    # build plotting data
    stla = []
    stlo = []
    vel = []
    for sta in sta_dict:
        stla.append(sta['stla'])
        stlo.append(sta['stlo'])
        vel.append(sta['sampled_vel'][i])
        
    print('max', max(vel), 'min', min(vel))
    
    negidx = where(array(vel) < 0)[0]
    log_vel = log10(abs(array(vel)))
    log_vel[negidx] *= -1
    idx = where(isnan(log_vel) == False)[0]
    vel = array(vel)
    #x, y = m(array(stlo)[idx], array(stla)[idx])
    #m.scatter(x, y, c=log_vel[idx], marker='o', s=30, cmap='seismic', norm=norm, alpha=1.0, zorder=1000)
    
    x, y = m(array(stlo), array(stla))
    m.scatter(x, y, c=vel, marker='o', s=40, cmap='seismic', norm=norm, alpha=1.0, zorder=1000)
    
    # add p/s wave trains
    if mt >= 0.0:
       # add epicentre star
       x, y = m(ev['lon'], ev['lat'])
       m.plot(x, y, 'r*', ms=25, label='Epicentre')
       
       # add p-wave
       p_radius = interp(mt, ptimes, kmrng)
       px = []
       py = []
       for az in azimuths:
           xy = reckon(ev['lat'], ev['lon'], p_radius, az)
           px.append(xy[0])
           py.append(xy[1])
       
       x, y = m(px, py)
       m.plot(x, y, '-', c='royalblue', lw=1.5, label='P-Phase')
       
       # add s-wave
       s_radius = interp(mt, stimes, kmrng)
       px = []
       py = []
       for az in azimuths:
           xy = reckon(ev['lat'], ev['lon'], s_radius, az)
           px.append(xy[0])
           py.append(xy[1])
       
       x, y = m(px, py)
       m.plot(x, y, '-', c='darkorange', lw=1.5, label='S-Phase')
    
        # add legend
        plt.legend(loc=2, fontsize=14, numpoints=1)
        
    # add label
    sample_time = UTCDateTime(ev['datetime']) + mt
    xpos, ypos = m(152.5, -27)
    plt.text(xpos, ypos, sample_time.isoformat()[0:-7], ha='right', va='bottom', fontsize=16, bbox=props)
    
    
    print('mapped_velocity_'+str(i+1)+'.png')
    if i < 9:
        plt.savefig('png/mapped_velocity_00'+str(i+1)+'.png', fmt='png', bbox_inches='tight')
    elif i < 99:
        plt.savefig('png/mapped_velocity_0'+str(i+1)+'.png', fmt='png', bbox_inches='tight')
    else:
        plt.savefig('png/mapped_velocity_'+str(i+1)+'.png', fmt='png', bbox_inches='tight')
    plt.clf()

###############################################################################
# make animation
###############################################################################
import matplotlib.image as mgimg
from matplotlib import animation

fig = plt.figure()

# initiate an empty  list of "plotted" images 
pngfiles = listdir_extension('png', 'png')
myimages = []

#loops through available png:s
for p in pngfiles:

    ## Read in picture
    fname = path.join('png', p) 
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img)

    # append AxesImage object to the list
    myimages.append([imgplot])

## create an instance of animation
my_anim = animation.ArtistAnimation(fig, myimages, interval=1000, blit=True, repeat_delay=1000)

## NB: The 'save' method here belongs to the object you created above
my_anim.save("animation.mp4")

## Showtime!
plt.show()
    
