import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import scipy.linalg as la

data = np.loadtxt('leon.txt')
models = np.atleast_2d(np.loadtxt('mean_covariance.txt'))

num_points = int((-1 + np.sqrt(1. + 4.*models.shape[1]))/2.)

plt.ion()
for i in xrange(0, models.shape[0]):
	model = models[i, :]
	mu = model[0:num_points]
	C  = model[num_points:].reshape((num_points, num_points))
	L = np.matrix((la.cholesky(C)).T)
	y = mu + L*np.matrix(rng.randn(num_points, 1))

	plt.hold(False)
	plt.plot(data[:,0], data[:,1], 'b.', markersize=1)
	plt.hold(True)
	plt.plot(data[:,0], y, 'r.', markersize=1)
	plt.draw()

plt.ioff()
plt.show()

