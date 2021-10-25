import pickle
import shapefile
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import matplotlib as mpl
from mapping_tools import get_field_data

mpl.style.use('classic')

def return_mag_dist(recs, polygons):
    maxdist = 1750.
    for j, poly in enumerate(polygons):
        
        if zone_code[j] == 'BS':
            mag = []
            rhyp = []
            dep = []
            vs30 = []
                    
            for rec in recs: 
                if rec['rhyp'] <= maxdist:
                    pt = Point(rec['eqlo'], rec['eqla'])
                    if pt.within(poly):
                        mag.append(rec['mag'])
                        rhyp.append(rec['rhyp'])
                        dep.append(rec['dep'])
                        vs30.append(rec['vs30'])
    
    return mag, rhyp, dep, vs30

##########################################################################################
# load shapefile
##########################################################################################

shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
zone_code = get_field_data(sf, 'CODE', 'str')
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

##########################################################################################
# parse amp SA files
##########################################################################################
# reads sa data files and returns period (T) and acceleration (SA) vectors
print('Loading pkl file...')
recs = pickle.load(open("stdict_ampfact.pkl", "rb" ))

amag, arhyp, adep, avs30 = return_mag_dist(recs, polygons)

##########################################################################################
# parse base SA files
##########################################################################################
# reads sa data files and returns period (T) and acceleration (SA) vectors
print('Loading pkl file...')
recs = pickle.load(open("stdict.pkl", "rb" ))

bmag, brhyp, bdep, bvs30 = return_mag_dist(recs, polygons)

##########################################################################################
# plot dist
##########################################################################################
letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

fig = plt.figure(1, figsize=(24,7))

ax = plt.subplot(131)
plt.plot(brhyp, bmag, 'o', mfc='dodgerblue', mec='none', ms=7, label='Base & Site Model', zorder=100)
plt.plot(arhyp, amag, 'o', mec='orangered', mfc='none', mew=2, ms=7.5, label='Site Model Only', zorder=1)


plt.xlabel('Hypocentral Distance (km)', fontsize=17)
plt.ylabel('Moment Magnitude', fontsize=17)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(which='both')

plt.legend(loc=2, fontsize=14, numpoints=3)

# add letter
ylims = ax.get_ylim()
xlims = ax.get_xlim()
plt.text(xlims[0]-0.04*xlims[1], ylims[1]+(ylims[1]-ylims[0])*0.05, '(a)', va='bottom', ha ='right', fontsize=20)

##########################################################################################
# plot dep
##########################################################################################

ax = plt.subplot(132)
plt.plot(bdep, bmag, 'o', mfc='dodgerblue', mec='none', ms=8, label='Base & Site Model', zorder=100)
#plt.plot(adep, amag, 'o', mec='orangered', mfc='none', mew=2, ms=7.5, label='Site Model Only', zorder=1)

plt.xlabel('Hypocentral Depth (km)', fontsize=17)
plt.ylabel('Moment Magnitude', fontsize=17)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(which='both')

ylims = ax.get_ylim()
xlims = ax.get_xlim()
plt.text(xlims[0]-0.04*xlims[1], ylims[1]+(ylims[1]-ylims[0])*0.05, '(b)', va='bottom', ha ='right', fontsize=20)

##########################################################################################
# plot vs30
##########################################################################################

ax = plt.subplot(133)
plt.plot(bvs30, bmag, 'o', mfc='dodgerblue', mec='none', ms=7, label='Base & Site Model', zorder=100)
plt.plot(avs30, amag, 'o', mec='orangered', mfc='none', mew=2, ms=7.5, label='Site Model Only', zorder=1)

plt.xlabel('$\mathregular{V_{S30}}$ (m/s)', fontsize=17)
plt.ylabel('Moment Magnitude', fontsize=17)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(which='both')

ylims = ax.get_ylim()
xlims = ax.get_xlim()
plt.text(xlims[0]-0.04*xlims[1], ylims[1]+(ylims[1]-ylims[0])*0.05, '(c)', va='bottom', ha ='right', fontsize=20)

##########################################################################################
# finish
##########################################################################################

plt.savefig('figures/2021_banda_mag_dist.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig_4.eps', fmt='eps', dpi=300, bbox_inches='tight')
plt.show()
