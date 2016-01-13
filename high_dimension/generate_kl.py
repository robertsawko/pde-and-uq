from numpy import logspace, linspace, pi, sin
from numpy.random import randn
# from kl_expansion import kl_covariance
from kl_expansion import kl_expansion_fourier_projected
import pickle
from matplotlib.pyplot import figure
from scipy.signal import welch
from generate_shinozuka import ftransform_analytical


m = 2000
eta = randn(m)
eta0 = randn(1)
a = pi
N = 4000
b = 100
x = linspace(-a, a, N)
kl = kl_expansion_fourier_projected(x, m=m, n=20, a=a, c=b, eta=eta)
y3 = 0.5 * (1.0 + eta0 + sin(x)) + kl

with open('u0.pickle', 'wb') as f:
    pickle.dump([x, y3], f)

fig = figure()
ax = fig.gca()
l3, = ax.plot(x, y3)

fig.show()
figc = figure()
f, pxx = welch(kl, fs=1 / (x[1] - x[0]), nperseg=N // 10)
axc = figc.gca()
axc.loglog(
    2 * pi * f, pxx / 4 / pi,
    label='Empirical m={0:n}'.format(m))
omega = logspace(-1, 4, 10000)
a = ftransform_analytical(omega, b=b)
axc.loglog(omega, a, label='Analytical')
figc.show()
