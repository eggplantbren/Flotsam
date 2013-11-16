from pylab import *

t = 1.05**arange(0., 100.)
output = loadtxt('output.txt')

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	plot(t, output[i, :], 'bo-')
	draw()

ioff()
show()

