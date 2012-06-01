import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from Data import Data

data = Data()
data.load('j1131.txt')

class Limits:
	"""
	A singleton class that just holds prior bounds
	"""
	mMin = data.yMean - 10.*data.yStDev
	mMax = data.yMean + 10.*data.yStDev
	mRange = mMax - mMin

	tauMin = -data.tRange
	tauMax =  data.tRange
	tauRange = tauMax - tauMin

class TDModel:
	"""
	An object of this class is a point in the Flotsam parameter space
	for inferring time delays in the presence of microlensing.
	"""
	def __init__(self):
		pass

	def fromPrior(self):
		# Mean magnitudes
		self.m = Limits.mMin + Limits.mRange*rng.rand(data.numImages)

		# Time delays
		self.tau = Limits.tauMin\
				+ Limits.tauRange*rng.rand(data.numImages)
		self.tau[0] = 0.

	@property
	def logPrior(self):
		logP = 0.
		if np.any(np.logical_or(self.m < Limits.mMin,\
				self.m > limits.mMax)):
			logP = -np.inf
		if np.any(np.logical_or(self.tau < Limits.tauMin,\
				self.tau > limits.tauMax)):
			logP = -np.inf
		return logP

	@property
	def logLikelihood(self):
		logL = 0.
		return logL

	@property
	def vector(self):
		"""
		Stack all parameters into a vector
		"""
		return np.hstack([self.m, self.tau])

if __name__ == '__main__':
	model = TDModel()
	model.fromPrior()
	print(model.m)

