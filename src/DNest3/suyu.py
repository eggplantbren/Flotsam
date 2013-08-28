from pylab import *

#This system is special
#in that the three relative time delays between the four images were measured accurately with errors of only a few
#percent: tAB = 31.5
#+2.0
#-1.0 days, tCB = 36.0
#+1.5
#-1.5 days,
#and tDB = 77.0
#+2.0
#-1.0 days (Fassnacht et al. 1999, 2002).

# IMAGE ID 0 = A
# IMAGE ID 1 = B
# IMAGE ID 2 = C
# IMAGE ID 3 = D

posterior_sample = loadtxt('posterior_sample.txt')

tAB = posterior_sample[:,4] - posterior_sample[:,5]
tCB = posterior_sample[:,6] - posterior_sample[:,5]
tDB = posterior_sample[:,7] - posterior_sample[:,5]

hist(tAB, 50, alpha=0.2, color='b')
hist(tCB, 50, alpha=0.2, color='r')
hist(tDB, 50, alpha=0.2, color='g')

# Find 68% credible intervals
def summary(theta):
	s = sort(theta)
	a = s[int(0.16*len(theta))]
	b = s[int(0.5*len(theta))]
	c = s[int(0.84*len(theta))]
	return (b - a, b, c - b)

summary1 = summary(tAB)
summary2 = summary(tCB)
summary3 = summary(tDB)

axhline(0, linewidth=2, color='k')

# Plot summaries
errorbar([summary1[1]], [-0.3], xerr=[[summary1[0]], [summary1[2]]], fmt='bo', linewidth=2)
errorbar([summary2[1]], [-0.2], xerr=[[summary2[0]], [summary2[2]]], fmt='ro', linewidth=2)
errorbar([summary3[1]], [-0.1], xerr=[[summary3[0]], [summary3[2]]], fmt='go', linewidth=2)

# Plot Fassnacht's time delays
errorbar([31.5], [-0.5], xerr=[[1.], [2.]], fmt='bo', linewidth=2)
errorbar([36.], [-0.5], xerr=[[1.5], [1.5]], fmt='ro', linewidth=2)
errorbar([77.], [-0.5], xerr=[[1.], [2.]], fmt='go', linewidth=2)

ylim(-1)
xlabel('Time Delay $\\tau$ (days)')
show()

print(mean(tAB), std(tAB))
print(mean(tCB), std(tCB))
print(mean(tDB), std(tDB))

print(summary1)
print(summary2)
print(summary3)




