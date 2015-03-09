from proteus import *
from proteus.default_p import *
from proteus.TransportCoefficients import TC_base
import numpy as np

"""
1D, burgers equation with Karhunen-Loeve initial condition
"""

nd = 1
L = [2 * np.pi]


class ForcedBurgersEqn(TC_base):
    def __init__(self, nu=0.0, sigma=0.05, rofx=lambda x, t: 0):
        TC_base.__init__(
            self, nc=1,
            mass={0: {0: 'linear'}},
            advection={0: {0: 'nonlinear'}},
            diffusion={0: {0: {0: 'constant'}}},
            potential={0: {0: 'nonlinear'}},
            reaction={0: {0: 'constant'}},
            hamiltonian={},
            variableNames=['u']
            )
        self.nu = nu
        self.rofx = rofx
        self.sigma = sigma

    def evaluate(self, t, c):
        u = c[('u', 0)]
        c[('m', 0)][:] = u
        c[('dm', 0, 0)][:] = 1
        c[('f', 0)][..., 0] = 0.5*u**2
        c[('df', 0, 0)][..., 0] = u
        c[('a', 0, 0)][..., 0, 0] = self.nu
        # x = c['x'][:, 0][:, 0]
        c[('r', 0)][:] = -self.sigma * self.rofx(c['x'][:], t)


class sineIC:
    def __init__(self, eta):
        self.eta = eta

    def uOfXT(self, x, t):
        return self.eta + np.sin(x[0])


def forcing_term(x, t, xi):
    return xi*np.sin(t)


def emptyBC(x):
    return None

xi = np.pi / 10 * np.random.randn()
eta = np.random.randn()

coefficients = ForcedBurgersEqn(
    nu=0,
    rofx=lambda x, t: forcing_term(x, t, xi),
    sigma=1
)

initialConditions = {0: sineIC(eta)}

dirichletConditions = {0: emptyBC}
advectiveFluxBoundaryconditions = {0: emptyBC}
diffusiveFluxBoundaryconditions = {0: {0: emptyBC}}


def getPBC(x):
    eps = 1.0e-8
    if x[0] < eps or x[0] >= L[0] - eps:
        return np.array([0.0, 0.0, 0.0])

periodicDirichletConditions = {0: getPBC}
T = 1
