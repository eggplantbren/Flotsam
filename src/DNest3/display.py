import numpy as np
import matplotlib.pyplot as plt

sample = np.atleast_2d(np.loadtxt('sample.txt'))

for i in xrange(0, sample.shape[1]):
	plt.plot(sample[:,0])
	plt.title(i)
	plt.show()

