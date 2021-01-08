import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import lfilter

f = open("caursarea_2.txt", "r")
g = open("caurspressure_2.txt", "r")

a = []
px = []
py = []

sigma = 52.98
epsilon = 275.7

for line in f:
	l = line.strip().split()
	a.append(float(l[0]))

for line in g:
	l = line.strip().split(",")
	px.append(float(l[0]))
	py.append(float(l[1]))

a = [i / (17056) for i in a]

fig = plt.figure(figsize=(6,8))
ax = fig.add_subplot(2,1,2)

n = 100  # the larger n is, the smoother curve will be
b = [1.0 / n] * n
aa = 1
ppxx = lfilter(b,aa,px)
ppyy = lfilter(b,aa,py)


ax.plot(a, ppxx, 'r')
ax.plot(a, ppyy, 'b', linestyle='--', dashes=(3,3))
ax.legend([r'$\Pi_{||}$', r'$\Pi_{\perp}$'])
ax.set_ylabel("Moduli (mN/m)", weight='bold')
ax.set_xlabel("Area per nanoparticle " + r'$(\AA^2 / nanoparticle)$', weight='bold')

#perform the derivatives

def derivative(x, y):
	x_out = []
	y_out = []

	for i in range(1, len(x)-1):
		t1 = y[i-1] * (2 * x[i] - x[i] - x[i+1]) / ( (x[i-1]-x[i])*(x[i-1]-x[i+1]) )
		t2 = y[i] * (2 * x[i] - x[i-1] - x[i+1]) / ( (x[i]-x[i-1])*(x[i]-x[i+1]) )
		t3 = y[i+1] * (2 * x[i] - x[i-1] - x[i]) / ( (x[i+1]-x[i-1])*(x[i+1]-x[i]) )
		x_out.append(x[i])
		y_out.append(t1 + t2 + t3)

	return np.array(x_out), np.array(y_out)

parx, pary = derivative(a, ppyy)
perpx, perpy = derivative(a, ppxx)


K = -.5 * parx * (pary + perpy) 
G = -.5 * parx * (pary - perpy)

skip = 100
parx = parx[skip:]
K = K[skip:]
G = G[skip:]

KK = lfilter(b,aa,K)
GG = lfilter(b,aa,G)

ax = fig.add_subplot(2,1,1)
ax.plot(parx, KK, 'k')
ax.plot(parx, GG, 'g', linestyle='--', dashes=(3,3))
ax.legend(["K", "G"])
ax.set_ylabel("Moduli (mN/m)", weight='bold')
#plt.title("Simulated Isotherm and Moduli", weight='bold', fontsize=16)
plt.savefig("isotherm", dpi=300)









