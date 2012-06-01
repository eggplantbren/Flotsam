import numpy as np
import matplotlib.pyplot as plt

class Data:
	"""
	An object of this class is a dataset in Flotsam time delay inference
	format. Columns: (time, magnitude, sig_magnitude, image ID)
	"""

	def __init__(self):
		"""
		Sets loaded to False, does not construct data arrays
		"""
		self.loaded = False
		pass

	def load(self, filename):
		"""
		Loads the data from the given file.
		"""
		data = np.loadtxt(filename)
		self.t = data[:,0]
		self.y = data[:,1]
		self.sig = data[:,2]
		self.id = data[:,3].astype('int')
		assert self.id.min() == 0
		self.computeSummaries()
		self.loaded = True

	def plot(self):
		for i in xrange(self.id.min(), self.id.max() + 1):
			which = np.nonzero(self.id == i)[0]
			plt.errorbar(self.t[which], self.y[which],\
					yerr=self.sig[which], fmt='.')
		plt.show()

	def computeSummaries(self):
		"""
		Compute useful summary statistics
		"""
		self.numImages = self.id.max() + 1
		self.tRange = self.t.max() - self.t.min()

		# Summary statistics for each image
		self.yMean = np.empty(self.numImages)
		self.ySqMean = np.empty(self.numImages)
		self.yStDev = np.empty(self.numImages)
		for i in xrange(0, self.numImages):
			which = np.nonzero(self.id == i)[0]
			w = self.sig[which]**(-2)
			w = w/w.sum()
			self.yMean[i] = (w*self.y[which]).sum()
			self.ySqMean[i] = (w*self.y[which]**2).sum()
			self.yStDev[i] = np.sqrt(self.ySqMean[i]\
							-self.yMean[i]**2)

if __name__ == '__main__':
	"""
	This is just for simple testing purposes.
	"""
	
	data = Data()
	data.load('j1131.txt')
	print(data.tRange, data.yMean, data.yStDev)
	data.plot()

