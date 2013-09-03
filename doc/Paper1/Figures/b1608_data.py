from pylab import *

rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

data = loadtxt('../../../src/DNest3/1608_all_mags.txt')

data[:,0] -= data[:,0].min()

figure(figsize=(16, 8))
colors = ['b', 'r', 'g', 'y']
shapes = ['o', 's', 'p', '*']
sizes  = [6, 5, 6, 7]
names  = ['A', 'B', 'C', 'D']
for i in xrange(0, 4):
	which = data[:,3] == i
	plot(data[which, 0], data[which, 1], shapes[i], color=colors[i], markersize=sizes[i], label='Image {I}'.format(I=names[i]))
legend(loc='upper left', numpoints=1)
axis([-50, 1250, -4, 0])
title('B1608+656 Data')

xlabel('$t$ (days)')
ylabel('$Y$ (magnitudes)')

savefig('b1608_data.pdf', bbox_inches='tight')
import os
os.system('pdftops -eps b1608_data.pdf')
show()

