import numpy as np
import emcee
from Data import Data
from TDModel import Limits, TDModel

# Define the Density
class Flotsam(object):
	def __init__(self, data):
		self.data = data
		self.model = TDModel(data.numImages)

	def __call__(self, params):
		self.model.fromVector(params)
		return self.model.logPrior + self.model.logLikelihood(data)

data = Data()
data.load('j1131.txt')
Limits.initialise(data)

nwalkers = 100

# Make an initial guess for the positions.
m = TDModel(data.numImages)
m.fromPrior()

params = np.empty((nwalkers, m.vector.size))
for i in xrange(0, nwalkers):
	m.fromPrior()
	params[i, :] = m.vector

# Instantiate the class
flotsam = Flotsam(data)

# The sampler object
sampler = emcee.EnsembleSampler(nwalkers, m.vector.size, flotsam, threads=1)

# Sample, outputting to a file
f = open("flotsam.out", "w")
for pos, prob, rstate in sampler.sample(params, iterations=1000):
	# Write the current position to a file, one line per walker
	f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
	f.write("\n")
	f.flush()
	print('Acceptance Fraction = %.3f'%sampler.acceptance_fraction.mean())

f.close()

