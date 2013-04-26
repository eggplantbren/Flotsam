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

which = rng.randint(9)
plt.figure()
for i in xrange(0, 9):
	model = models[i, :]
	mu = np.matrix(model[0:num_points]).T
	C  = model[num_points:].reshape((num_points, num_points))
	L = np.matrix((la.cholesky(C)).T)

	y = mu + L*np.matrix(rng.randn(num_points, 1))

	plt.subplot(3,3,i+1)
	if i == which:
		plt.plot(data[:,0], data[:,1], 'b.', markersize=1)
	else:
		plt.plot(data[:,0], y, 'b.', markersize=1)
	plt.xlim([data[:,0].min(), data[:,0].max()])
	plt.gca().set_xticklabels([])
	plt.gca().set_yticklabels([])
#	plt.title('{count}/{total}'.format(count=(3*i + k + 1), total=3*models.shape[0]))
#	plt.legend(numpoints=1, loc='lower left')
#	plt.xlabel('Time')
#	plt.ylabel('Magnitude')

plt.subplot(3,3,8)
plt.xlabel(r'$t$ (days)')
plt.subplot(3,3,4)
plt.ylabel(r'$y$ (magnitudes)')

plt.tight_layout()
plt.savefig('posterior_predictive_check.pdf', bbox_inches='tight')
plt.show()

