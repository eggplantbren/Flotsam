import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from Data import Data

class Limits:
	"""
	A singleton class that just holds prior bounds
	"""
	@staticmethod
	def initialise(data):
		Limits.mMin = data.yMean - 10.*data.yStDev
		Limits.mMax = data.yMean + 10.*data.yStDev
		Limits.mRange = Limits.mMax - Limits.mMin

		Limits.tauMin = -data.tRange
		Limits.tauMax =  data.tRange
		Limits.tauRange = Limits.tauMax - Limits.tauMin

class TDModel:
	"""
	An object of this class is a point in the Flotsam parameter space
	for inferring time delays in the presence of microlensing.
	"""
	def __init__(self, numImages):
		self.numImages = numImages

	def fromPrior(self):
		# Mean magnitudes
		self.m = Limits.mMin + Limits.mRange*rng.rand(self.numImages)

		# Time delays
		self.tau = Limits.tauMin\
				+ Limits.tauRange*rng.rand(self.numImages)
		self.tau[0] = 0.

	@property
	def logPrior(self):
		logP = 0.
		if np.any(np.logical_or(self.m < Limits.mMin,\
				self.m > Limits.mMax)):
			logP = -np.inf
		if np.any(np.logical_or(self.tau < Limits.tauMin,\
				self.tau > Limits.tauMax)):
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

	def fromVector(self, vec):
		"""
		Set all parameters from the vector
		"""
		self.m = vec[0:self.numImages]
		self.tau = vec[self.numImages:2*self.numImages]
		

if __name__ == '__main__':
	data = Data()
	data.load('j1131.txt')
	Limits.initialise(data)

	model = TDModel(data.numImages)
	model.fromPrior()
	print(model.m)

