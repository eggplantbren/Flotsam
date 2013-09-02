import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import scipy

"""
This little code makes a magnification map, so I can generate
microlensing light curves for battle testing Flotsam.
"""

# Reproducibility
rng.seed(123)

# Number of rays per dimension
N = 4000

# The rectangle in the lens plane
x_min, x_max = -10., 10.
y_min, y_max = -10., 10.
dx = (x_max - x_min)/N
dy = (y_max - y_min)/N

# A grid of rays
x = np.linspace(x_min + 0.5*dx, x_max - 0.5*dx, N)
y = np.linspace(y_min + 0.5*dy, y_max - 0.5*dy, N)
[x, y] = np.meshgrid(x, y)
y = y[::-1, :]

# Convergence, external shear, number of stars
# Note: "gamma" DOES NOTHING
rho, gamma, num_stars = 0.5, 0., 500

# Generate the stars, in a circle
radius = np.min([0.5*(x_max - x_min), 0.5*(y_max - y_min)])

# Mass of each star
mass = np.pi*radius**2*rho/num_stars

x_stars = np.empty(num_stars)
y_stars = np.empty(num_stars)
for i in xrange(0, num_stars):
	x_stars[i] = x_min + (x_max - x_min)*rng.rand()
	y_stars[i] = y_min + (y_max - y_min)*rng.rand()
	rsq = x_stars[i]**2 + y_stars[i]**2
	while rsq > radius**2:
		x_stars[i] = x_min + (x_max - x_min)*rng.rand()
		y_stars[i] = y_min + (y_max - y_min)*rng.rand()
		rsq = x_stars[i]**2 + y_stars[i]**2

plt.plot(x_stars, y_stars, 'r*')
plt.title('The Stars')
plt.axis('scaled')
plt.show()


# Compute source plane positions
xs = x.copy()
ys = y.copy()

plt.ion()
plt.hold(False)
for i in xrange(0, num_stars):
	rsq = (x - x_stars[i])**2 + (y - y_stars[i])**2
	xs -= mass*(x - x_stars[i])/(np.pi*rsq)
	ys -= mass*(y - y_stars[i])/(np.pi*rsq)

	print('{i}/{n}'.format(i=(i+1), n=num_stars))

	if (i+1)%50 == 0:
		counts = scipy.histogram2d(xs.flatten(), ys.flatten(), bins=250,\
					range=[[x_min, x_max], [y_min, y_max]])
		plt.imshow(counts[0])
		plt.title('{a}/{b}'.format(a=(i+1), b=num_stars))
		plt.draw()

rays = np.empty((xs.size, 2))
rays[:,0] = xs.flatten()
rays[:,1] = ys.flatten()

# Save the rays
import cPickle as pickle
print('Saving...')
pickle.dump(rays, open('rays.pickle', 'wb'))
print('Done.')

plt.ioff()
plt.show()

