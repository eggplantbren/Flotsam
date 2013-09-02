from pylab import *

def ar1(n=10000, a=0.999, b=1.):
	y = zeros(n)
	for i in xrange(1, n):
		y[i] = a*y[i-1] + b*randn()
	return y

seed(123)

qso = ar1()/20.

# Observation times
y1 = qso[8100:9100:10].copy()
y2 = qso[8000:9000:10].copy()
N = y1.size

# Offset
y2 += 2.

# Invent timestamps
t = arange(0, N)

# Add noise
y1 += 0.3*randn(N)
y2 += 0.3*randn(N)

# Make Flotsam data
data = empty((2*N, 4))

data[0:N, 0] = t
data[0:N, 1] = y1
data[0:N, 2] = 0.3
data[0:N, 3] = 0
data[N:, 0] = t
data[N:, 1] = y2
data[N:, 2] = 0.3
data[N:, 3] = 1

errorbar(data[0:N, 0], data[0:N, 1], yerr=data[0:N, 2], fmt='bo')
errorbar(data[N:, 0], data[N:, 1], yerr=data[N:, 2], fmt='ro')
show()

savetxt('easy_data.txt', data)
