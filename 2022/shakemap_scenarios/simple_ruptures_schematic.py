import shapefile
#from mapping_tools import get_field_data, distance, reckon, get_line_parallels
from fault_tools import mag2wid_L10, mag2rupwid_WC94, mag2wid_L14
from numpy import radians, cos, sin, tan, mean, array
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

mw = 6.7
dip = 35.
dep = 10.

rupwid = mag2wid_L14(mw, 'scrrs')

rngkm = rupwid * cos(radians(dip))
rupdep = rupwid * sin(radians(dip))
dep_diff = dep - rupdep

if rupdep > dep:
    
    posazpts, negazpts = get_line_parallels(shape.points, rngkm)
    

# shuffle down such that hypo dep captured in rup dep
else:
  # do +ve
    hdist2tor = dep_diff / tan(radians(dip))        
                
    hdist2bor = dep / tan(radians(dip))        
           
# now plot
fig = plt.figure(1, figsize=(8,4))
ax = fig.add_subplot(111)
ax.set_aspect('equal')

# plt h levels
plt.plot([-5, hdist2bor+5], [0, 0], '-', c='0.5')
plt.plot([-5, hdist2bor+5], [-dep_diff, -dep_diff], '--', c='0.5')
plt.plot([-5, hdist2bor+5], [-dep, -dep], '--', c='0.5')
plt.plot([15, 15], [-dep_diff, -dep], '--', c='0.5')

# plt fault
plt.plot([0, hdist2bor], [0, -dep], 'k--', lw=2., label='Fault Surface Projection')
plt.plot([hdist2tor, hdist2bor], [-dep_diff, -dep], 'r-', lw=2., label='Modelled Rupture Width')
plt.legend(loc=1, fontsize=10)

# write text
plt.text(-5, 0.1, 'Ground Surface', va='bottom', fontsize=12)
plt.text(-5, -dep_diff+0.1, '$\mathregular{z_{tor}}$', va='bottom', fontsize=14)
plt.text(-5, -dep+0.1, '$\mathregular{z_{bor}}$', va='bottom', fontsize=14)
plt.text(15.1, ((-dep + -dep_diff)/2.), r'$\Delta$'+'z', va='center', fontsize=14)

plt.xlim([-7, hdist2bor+7])
plt.ylim([-dep-2, 3])
plt.xlabel('Fault-Normal Distance (km)', fontsize=15)
plt.ylabel('Depth (km)', fontsize=15)

plt.savefig('rupture_schematic.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()