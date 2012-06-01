import numpy as np
import matplotlib.pyplot as plt
import time

sample = np.loadtxt('flotsam.out')
for i in xrange(0, sample.shape[1]):
	plt.hist(sample[:,i], 20)
	plt.title(i)
	plt.show()

