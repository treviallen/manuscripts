from sys import argv
from numpy import array, nanmean, nanstd
import matplotlib.pyplot as plt
from calculate_magnitudes import calc_WGW96, calc_A10
import matplotlib as mpl
mpl.style.use('classic')

mlfile = 'moe_stn_mags_sansTAS_54.csv'
#mlfile = '/Users/tallen/Documents/Code/process_waves/ml/201207200911.ml'
#mlfile = argv[1]

lines = open(mlfile).readlines()[1:]

stn = []
logA = []
rhyp = []
r35 = []
hb87 = []
mlm92 = []
wgw96 = []
a10 = []

for line in lines:
    dat = line.strip().split(',')
    stn.append(dat[1])
    #logA.append(float(dat[2]))
    rhyp.append(float(dat[2]))
    r35.append(float(dat[3]))
    hb87.append(float(dat[6]))
    mlm92.append(float(dat[8]))
    a10.append(float(dat[10]))
    wgw96.append(float(dat[9]))
        
# plot ML with distance
fig = plt.figure(1, figsize=(12.75, 5))
ax = plt.subplot(111)
plt.plot(rhyp, r35, 'ro', ms=10)
plt.plot(rhyp, hb87, 'b^', ms=10)
plt.plot(rhyp, mlm92, 'gs', ms=10)
plt.plot(rhyp, wgw96, 'd', color='orange', ms=10)
plt.plot(rhyp, a10, 'cv', ms=10)

plt.legend(['R35', 'HB87','MLM92','WGW96','A10'],loc=2, fontsize=12, numpoints=1)
plt.xlabel('Hypocentral Distance (km)')
plt.ylabel('Local Magnitude')

pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0,  pos1.width * 0.93, pos1.height] 
ax.set_position(pos2) # set a new position

print('R35  ', nanmean(r35), nanstd(r35))
print('HB87 ', nanmean(hb87), nanstd(hb87))
print('MLM92', nanmean(mlm92), nanstd(mlm92))
print('WGW96', nanmean(wgw96), nanstd(wgw96))
print('A10  ', nanmean(a10), nanstd(a10))

plt.plot([0, 500],[nanmean(r35), nanmean(r35)], 'r--', lw=2)
plt.plot([0, 500],[nanmean(hb87), nanmean(hb87)], 'b--', lw=2)
plt.plot([0, 500],[nanmean(mlm92), nanmean(mlm92)], 'g--', lw=2)
plt.plot([0, 500],[nanmean(wgw96), nanmean(wgw96)], '--', color='orange', lw=2)
plt.plot([0, 500],[nanmean(a10), nanmean(a10)], 'c--', lw=2)

# annotate nanmean + nanstd
plt.text(503, nanmean(r35), 'R35: '+ str('%0.2f' % nanmean(r35)) + ' (' + str('%0.2f' % nanstd(r35))+')', ha='left', va='center', fontsize=11)
plt.text(503, nanmean(hb87), 'HB87: '+ str('%0.2f' % nanmean(hb87)) + ' (' + str('%0.2f' % nanstd(hb87))+')', ha='left', va='center', fontsize=11)
plt.text(503, nanmean(mlm92)+0.01, 'MLM92: '+ str('%0.2f' % nanmean(mlm92)) + ' (' + str('%0.2f' % nanstd(mlm92))+')', ha='left', va='center', fontsize=11)
plt.text(503, nanmean(wgw96), 'WGW96: '+ str('%0.2f' % nanmean(wgw96)) + ' (' + str('%0.2f' % nanstd(wgw96))+')', ha='left', va='center', fontsize=11)
plt.text(503, nanmean(a10)-0.01, 'A10: '+ str('%0.2f' % nanmean(a10)) + ' (' + str('%0.2f' % nanstd(a10))+')', ha='left', va='center', fontsize=11)

outpng = mlfile+'.png'
plt.savefig(outpng, format='png', bbox_inches='tight')
plt.show()
