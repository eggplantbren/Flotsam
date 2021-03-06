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

# Add microlensing to make hard data
# Load the rays
import cPickle as pickle
print('Loading...')
rays = pickle.load(open('rays.pickle', 'rb'))
print('Done.')

# FIRST IMAGE
x0 = -1.
y0 = -1.
vx = 0.01
vy = 0.01
r = 0.05

x = x0 + vx*arange(0, N)
y = y0 + vy*arange(0, N)

plot(rays[0::37,0], rays[0::37,1], 'b.', markersize=1)
plot(x, y, 'ro-', markersize=10)
axis([-4., 4., -4., 4.])
show()

mag = zeros(N)
ion()
hold(False)
for i in xrange(0, N):
	mag[i] = 2.5*log10(sum((rays[:,0] - x[i])**2 + (rays[:,1] - x[i])**2 <= r**2))
	plot(mag[0:(i+1)])
	draw()
ioff()
show()

data[0:N, 1] += mag - mean(mag)

# SECOND IMAGE
x0 = 1.
y0 = 1.
vx = 0.01
vy = -0.01
r = 0.05

x = x0 + vx*arange(0, N)
y = y0 + vy*arange(0, N)

hold(True)
plot(rays[0::37,0], rays[0::37,1], 'b.', markersize=1)
plot(x, y, 'ro-', markersize=10)
axis([-4., 4., -4., 4.])
show()

mag = zeros(N)
ion()
hold(False)
for i in xrange(0, N):
	mag[i] = 2.5*log10(sum((rays[:,0] - x[i])**2 + (rays[:,1] - x[i])**2 <= r**2))
	plot(mag[0:(i+1)])
	draw()
ioff()
show()
data[N:, 1] += mag - mean(mag)

savetxt('hard_data.txt', data)

hold(True)
errorbar(data[0:N, 0], data[0:N, 1], yerr=data[0:N, 2], fmt='bo')
errorbar(data[N:, 0], data[N:, 1], yerr=data[N:, 2], fmt='ro')
show()
