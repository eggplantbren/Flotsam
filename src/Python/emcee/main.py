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
		return model.logPrior + model.logLikelihood

data = Data()
data.load('j1131.txt')
Limits.initialise(data)

nwalkers = 100

# Make an initial guess for the positions.
model = TDModel(data.numImages)
model.fromPrior()
params = np.empty((nwalkers, model.vector.size))
for i in xrange(0, nwalkers):
	model.fromPrior()
	params[i, :] = model.vector

# Instantiate the class
flotsam = Flotsam(data)

# The sampler object
sampler = emcee.EnsembleSampler(nwalkers, model.vector.size, flotsam, threads=10)

# Sample, outputting to a file
f = open("flotsam.out", "w")
for pos, prob, rstate in sampler.sample(params, iterations=2000):
	# Write the current position to a file, one line per walker
	f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
	f.write("\n")
	f.flush()
f.close()

