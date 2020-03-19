from sys import argv
from numpy import array, mean, std
import matplotlib.pyplot as plt
from calculate_magnitudes import calc_WGW96, calc_A10

mlfile = '/Users/tallen/Documents/Code/process_waves/ml/201206191053.ml'
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
    dat = line.strip().split('\t')
    stn.append(dat[0])
    logA.append(float(dat[4]))
    rhyp.append(float(dat[2]))
    r35.append(float(dat[5]))
    hb87.append(float(dat[6]))
    mlm92.append(float(dat[7]))
    if len(dat) == 9:
        a10.append(float(dat[8]))
    else:
        a10.append(calc_A10(0, logA[-1], rhyp[-1]))
    wgw96.append(calc_WGW96(0, logA[-1], rhyp[-1]))
        
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

print 'R35  ', mean(r35), std(r35)
print 'HB87 ', mean(hb87), std(hb87)
print 'MLM92', mean(mlm92), std(mlm92)
print 'WGW96', mean(wgw96), std(wgw96)
print 'A10  ', mean(a10), std(a10)

plt.plot([0, 500],[mean(r35), mean(r35)], 'r--', lw=2)
plt.plot([0, 500],[mean(hb87), mean(hb87)], 'b--', lw=2)
plt.plot([0, 500],[mean(mlm92), mean(mlm92)], 'g--', lw=2)
plt.plot([0, 500],[mean(wgw96), mean(wgw96)], '--', color='orange', lw=2)
plt.plot([0, 500],[mean(a10), mean(a10)], 'c--', lw=2)

# annotate mean + std
plt.text(503, mean(r35)+0.01, 'R35: '+ str('%0.2f' % mean(r35)) + ' (' + str('%0.2f' % std(r35))+')', ha='left', va='center', fontsize=11)
plt.text(503, mean(hb87), 'HB87: '+ str('%0.2f' % mean(hb87)) + ' (' + str('%0.2f' % std(hb87))+')', ha='left', va='center', fontsize=11)
plt.text(503, mean(mlm92)-0.05, 'MLM92: '+ str('%0.2f' % mean(mlm92)) + ' (' + str('%0.2f' % std(mlm92))+')', ha='left', va='center', fontsize=11)
plt.text(503, mean(wgw96), 'WGW96: '+ str('%0.2f' % mean(wgw96)) + ' (' + str('%0.2f' % std(wgw96))+')', ha='left', va='center', fontsize=11)
plt.text(503, mean(a10), 'A10: '+ str('%0.2f' % mean(a10)) + ' (' + str('%0.2f' % std(a10))+')', ha='left', va='center', fontsize=11)

outpng = mlfile+'.png'
plt.savefig(outpng, format='png', bbox_inches='tight')
plt.show()
