import numpy as np

class Data
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
		self.loaded = True

	# Useful summary statistics
	@property
	def tRange(self):
		return self.t.max() - self.t.min()

	@property
	def weights(self):
		w = self.sig**(-2)
		w = w/w.sum()
		return w

	@property
	def yMean(self):
		return (self.weights*self.y).sum()

	@property
	def ySqMean(self):
		return (self.weights*self.y**2).sum()

	@property
	def yStDev(self):
		return np.sqrt(self.ySqMean - self.yMean**2)

if __name__ == '__main__':
	"""
	This is just for simple testing purposes.
	"""
	
	
