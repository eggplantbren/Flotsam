from pylab import *

rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

data = loadtxt('../../../src/Microlensing/easy_data.txt')
N = data.shape[0]//2

figure(figsize=(20, 6))
hold(True)

subplot(1,3,1)
errorbar(data[0:N, 0], data[0:N, 1], yerr=data[0:N, 2], fmt='bo', label='Image 1')
errorbar(data[N:, 0], data[N:, 1], yerr=data[N:, 2], fmt='ro', label='Image 2')
xlabel('$t$ (days)')
ylabel('$Y$ (magnitudes)')
title('Easy Data')
axis([-3, 103, -3, 4])
legend(loc='upper left', numpoints=1)

easy = data
data = loadtxt('../../../src/Microlensing/hard_data.txt')

subplot(1,3,2)
errorbar(data[0:N, 0], data[0:N, 1], yerr=data[0:N, 2], fmt='bo', label='Image 1')
errorbar(data[N:, 0], data[N:, 1], yerr=data[N:, 2], fmt='ro', label='Image 2')
xlabel('$t$ (days)')
ylabel('$Y$ (magnitudes)')
title('Hard Data')
axis([-3, 103, -3, 4])
legend(loc='upper left', numpoints=1)

subplot(1,3,3)
plot(data[0:N, 0], data[0:N, 1] - easy[0:N, 1], 'b', linewidth=2, label='Image 1')
plot(data[N:, 0], data[N:, 1] - easy[N:, 1], 'r', linewidth=2, label='Image 2')
xlabel('$t$ (days)')
ylabel('$\\mu(t)$ (magnitudes)')
title('Microlensing in Hard Data')
axis([-3, 103, -3, 4])
legend(loc='upper left', numpoints=1)

savefig('sim_data.pdf', bbox_inches='tight')
import os
os.system('pdftops -eps sim_data.pdf')
show()

