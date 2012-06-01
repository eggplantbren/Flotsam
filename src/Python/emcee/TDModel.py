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
		Limits.magMin = data.yMean - 10.*data.yStDev
		Limits.magMax = data.yMean + 10.*data.yStDev
		Limits.magRange = Limits.magMax - Limits.magMin

		Limits.tauMin = -data.tRange
		Limits.tauMax =  data.tRange
		Limits.tauRange = Limits.tauMax - Limits.tauMin

		Limits.logSig_mlMin = np.log(1E-2*data.yStDev)
		Limits.logSig_mlMax = np.log(1E2*data.yStDev)
		Limits.logSig_mlRange = Limits.logSig_mlMax\
						- Limits.logSig_mlMin

class TDModel:
	"""
	An object of this class is a point in the Flotsam parameter space
	for inferring time delays in the presence of microlensing.
	"""
	def __init__(self, numImages):
		self.numImages = numImages

	def fromPrior(self):
		# Mean magnitudes
		self.mag = Limits.magMin + Limits.magRange*rng.rand(self.numImages)

		# Time delays
		self.tau = Limits.tauMin\
				+ Limits.tauRange*rng.rand(self.numImages)
		self.tau[0] = 0.

		# Microlensing amplitudes
		self.logSig_ml = Limits.logSig_mlMin + Limits.logSig_mlRange\
					*rng.rand(self.numImages)

	@property
	def logPrior(self):
		logP = 0.
		if np.any(np.logical_or(self.mag < Limits.magMin,\
				self.mag > Limits.magMax)):
			logP = -np.inf
		if np.any(np.logical_or(self.tau < Limits.tauMin,\
				self.tau > Limits.tauMax)):
			logP = -np.inf
		if np.any(np.logical_or(self.logSig_ml < Limits.logSig_mlMin,\
				self.logSig_ml > Limits.logSig_mlMax)):
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
		return np.hstack([self.mag, self.tau, self.logSig_ml])

	def fromVector(self, vec):
		"""
		Set all parameters from the vector
		"""
		self.mag = vec[0:self.numImages]
		self.tau = vec[self.numImages:2*self.numImages]
		self.logSig_ml = vec[2*self.numImages:3*self.numImages]

if __name__ == '__main__':
	data = Data()
	data.load('j1131.txt')
	Limits.initialise(data)

	model = TDModel(data.numImages)
	model.fromPrior()

