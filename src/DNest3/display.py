import numpy as np
import matplotlib.pyplot as plt

sample = np.atleast_2d(np.loadtxt('sample.txt'))
sample_info = np.atleast_2d(np.loadtxt('sample_info.txt'))

for i in xrange(0, sample.shape[1]):
	plt.subplot(2,1,1)
	plt.plot(sample[:,i])
	plt.title(i)
	plt.subplot(2,1,2)
	plt.plot(sample_info[:,0], sample[:,i], 'k.', markersize=2)
	plt.title('Convergence Plot')
	plt.show()

