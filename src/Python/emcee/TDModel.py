import numpy as np
import numpy.random as rng
import scipy.linalg as la
import matplotlib.pyplot as plt
from Data import Data

class Limits:
	"""
	A singleton class that just holds prior bounds
	"""
	@staticmethod
	def initialise(data):
		# Mean magnitudes
		Limits.magMin = data.yMean - 10.*data.yStDev
		Limits.magMax = data.yMean + 10.*data.yStDev
		Limits.magRange = Limits.magMax - Limits.magMin

		# Time delays
		Limits.tauMin = -data.tRange
		Limits.tauMax =  data.tRange
		Limits.tauRange = Limits.tauMax - Limits.tauMin

		# Microlensing
		Limits.logSig_mlMin = np.log(1E-2*data.yStDev)
		Limits.logSig_mlMax = np.log(1E2*data.yStDev)
		Limits.logSig_mlRange = Limits.logSig_mlMax\
						- Limits.logSig_mlMin
		Limits.alpha_mlMin = 1.
		Limits.alpha_mlMax = 2.
		Limits.alpha_mlRange = Limits.alpha_mlMax\
						- Limits.alpha_mlMin

		# QSO Variability
		Limits.logSig_qsoMin = np.log(1E-2*data.yStDev).mean()
		Limits.logSig_qsoMax = np.log(1E2*data.yStDev).mean()
		Limits.logSig_qsoRange = Limits.logSig_qsoMax\
						- Limits.logSig_qsoMin
		Limits.logTau_qsoMin = np.log(1E-2*data.tRange)
		Limits.logTau_qsoMax = np.log(1E3*data.tRange)
		Limits.logTau_qsoRange = Limits.logTau_qsoMax\
						- Limits.logTau_qsoMin

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

		# Microlensing smoothness
		self.alpha_ml = Limits.alpha_mlMin + Limits.alpha_mlRange\
					*rng.rand()

		# QSO Variability
		self.logSig_qso = Limits.logSig_qsoMin +\
					Limits.logSig_qsoRange*rng.rand()
		self.logTau_qso = Limits.logTau_qsoMin +\
					Limits.logTau_qsoRange*rng.rand()

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
		if self.alpha_ml < Limits.alpha_mlMin or self.alpha_ml >\
					Limits.alpha_mlMax:
			logP = -np.inf
		if self.logSig_qso < Limits.logSig_qsoMin or\
			self.logSig_qso > Limits.logSig_qsoMax:
			logP = -np.inf
		if self.logTau_qso < Limits.logTau_qsoMin or\
			self.logTau_qso > Limits.logTau_qsoMax:
			logP = -np.inf
		return logP

	def logLikelihood(self, data):
		assert self.numImages == data.numImages

		which = [np.nonzero(data.id == i)[0]\
				for i in xrange(0, self.numImages)]

		# Mean vector
		m = np.empty(data.t.size)
		for i in xrange(0, self.numImages):
			m[which[i]] = self.mag[i]

		# Covariance matrix
		C = np.zeros((data.t.size, data.t.size))
		ids = np.meshgrid(data.id, data.id)

		equal = ids[0] == ids[1]

##		try:
##			L = la.cholesky(C)
##		except:
##			return -np.inf
#		y = data.y - m
#		logDeterminant = 2.0*np.sum(np.log(np.diag(L)))
#		solution = la.cho_solve((L, True), y)
#		exponent = np.dot(y, solution)
#		logL = -0.5*data.t.size*np.log(2.0*np.pi) - 0.5*logDeterminant - 0.5*exponent
		return 0.

	@property
	def vector(self):
		"""
		Stack all parameters into a vector
		"""
		return np.hstack([self.mag, self.tau, self.logSig_ml,\
					self.alpha_ml, self.logSig_qso,\
					self.logTau_qso])

	def fromVector(self, vec):
		"""
		Set all parameters from the vector
		"""
		self.mag = vec[0:self.numImages]
		self.tau = vec[self.numImages:2*self.numImages]
		self.logSig_ml = vec[2*self.numImages:3*self.numImages]
		self.alpha_ml = vec[3*self.numImages]
		self.logSig_qso = vec[3*self.numImages + 1]
		self.logTau_qso = vec[3*self.numImages + 2]

if __name__ == '__main__':
	data = Data()
	data.load('j1131.txt')
	Limits.initialise(data)

	model = TDModel(data.numImages)
	model.fromPrior()

