import numpy as np
import matplotlib.pyplot as plt

cut = 0.
sample = np.atleast_2d(np.loadtxt('sample.txt'))
sample_info = np.atleast_2d(np.loadtxt('sample_info.txt'))
cut = int(cut*sample.shape[0])
sample = sample[cut:]
sample_info = sample_info[cut:]

for i in xrange(0, sample.shape[1]):
	plt.subplot(2,1,1)
	plt.plot(sample[:,i])
	plt.title(i)
	plt.subplot(2,1,2)
	plt.plot(sample_info[:,0], sample[:,i], 'k.', markersize=2)
	plt.title('Convergence Plot')
	plt.show()

plt.rc("font", size=18, family="serif", serif="Computer Sans")
plt.rc("text", usetex=True)
posterior_sample = np.atleast_2d(np.loadtxt('posterior_sample.txt'))

plt.hist(posterior_sample[:,3], 30, normed=True, label='Posterior Samples')
plt.xlabel('Time Delay $\\tau$ (days)')
plt.ylabel('Posterior Probability')
plt.title('Time Delay = {mean:.2f} $\\pm$ {sd:.2f} days'\
		.format(mean=posterior_sample[:,3].mean(),
		sd=posterior_sample[:,3].std()))
plt.savefig('hist.pdf', bbox_inches='tight')
plt.show()

#import triangle
#index = np.array([11, 12, 3])
#labels = [r'$f_{\rm bad}$',
#		r'$K$', r'$\tau$']

#ndim, nsamples = index.size, posterior_sample.shape[0]
#triangle.corner(posterior_sample[:, index], labels=labels, quantiles=[0.16, 0.5, 0.84], bins=30, smooth=False)
#plt.show()

#plt.plot(posterior_sample[:,11], posterior_sample[:,12],
#		'bo', markersize=5, alpha=0.25, label='Posterior Samples')
#plt.xlabel(r'$f_{\rm bad}$')
#plt.ylabel(r'Error bar boost $K$')
#plt.xlim([-0.02, 1.02])
#plt.ylim([ -0.1, 12.1])
#plt.axhline(1., linestyle='--', color='k')
#plt.axvline(0., linestyle='--', color='k')
#plt.savefig('error_bars.pdf', bbox_inches='tight')
#plt.show()

