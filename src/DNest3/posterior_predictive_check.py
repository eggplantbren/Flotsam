import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import scipy.linalg as la
import time

plt.rc("font", size=18, family="serif", serif="Computer Sans")
plt.rc("text", usetex=True)

data = np.loadtxt('leon.txt')
models = np.atleast_2d(np.loadtxt('mean_covariance.txt'))

num_points = int((-1 + np.sqrt(1. + 4.*models.shape[1]))/2.)

plt.ion()
for i in xrange(0, models.shape[0]):
	model = models[i, :]
	mu = np.matrix(model[0:num_points]).T
	C  = model[num_points:].reshape((num_points, num_points))
	L = np.matrix((la.cholesky(C)).T)

	for k in xrange(0, 3):
		y = mu + L*np.matrix(rng.randn(num_points, 1))

		plt.hold(False)
		plt.plot(data[:,0], data[:,1], 'b.', markersize=1, label='Actual Data')
		plt.hold(True)
		plt.plot(data[:,0], y, 'r.', markersize=1, label='Simulated Data')
		plt.title('{count}/{total}'.format(count=(3*i + k + 1), total=3*models.shape[0]))
		plt.legend(numpoints=1, loc='lower left')
		plt.xlabel('Time')
		plt.ylabel('Magnitude')
		plt.draw()
		time.sleep(1)

plt.ioff()
plt.show()

