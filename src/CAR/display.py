from pylab import *

t = arange(0., 100.)
t[98] = 105.
t[99] = 110.
output = loadtxt('output.txt')

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	plot(t, output[i, :], 'bo-')
	draw()

ioff()
show()

