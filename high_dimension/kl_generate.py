import numpy as np
# from kl_expansion import kl_covariance
from kl_expansion import kl_expansion_fourier_projected
from kl_expansion import kl_expansion
from scipy.interpolate import interp1d
import pickle


m = 200
eta = np.random.randn(m)
eta0 = np.random.randn(1)
a = np.pi
x = np.linspace(-a, a, 4000)
#y1 = 0.5*(1.0 + eta0 + np.sin(x)) \
    #+ kl_expansion(x, m=m, a=a, c=1/0.01, eta=eta)
#y2 = 0.5*(1.0 + eta0 + np.sin(x)) \
    #+ kl_expansion_fourier_projected(x, m=m, n=10, a=a, c=1/0.01, eta=eta)
y3 = 0.5*(1.0 + eta0 + np.sin(x)) \
    + kl_expansion_fourier_projected(x, m=m, n=20, a=a, c=1/0.01, eta=eta)

# u0 = interp1d(x, y3, kind='cubic')


with open('u0.pickle', 'w') as f:
    pickle.dump([x, y3], f)


import matplotlib as mpl
# Define a backend
# MUST be invoked before pyplot import
mpl.use("pgf")
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,2))
ax = fig.gca()
#l1, = ax.plot(x, y1, alpha=0.5)
#l2, = ax.plot(x, y2)
l3, = ax.plot(x, y3)
plt.legend(
    [l3],
    [
     #   "Original KL",
     #   "Fourier projected KL n={0}".format(10),
        "Fourier projected KL n={0}".format(20)
    ],
    loc="upper left"
)

fig.savefig(
    "u0.pdf",
    transparent=True, bbox_inches='tight', pad_inches=0)

