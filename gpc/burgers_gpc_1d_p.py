from proteus import *
from proteus.default_p import *
from proteus.TransportCoefficients import TC_base
import numpy as np
import pickle

"""
1D, burgers equation with Karhunen-Loeve initial condition
"""

nd = 1
L = [2 * np.pi]
nc = 9

with open('e_div_gamma.pickle') as f:
    e_div_gamma = pickle.load(f)


class gPCBurgersEqn(TC_base):
    def __init__(self, nu=0.0, nc=3):
        advection = {}
        diffusion = {}
        mass = {}
        potential = {}
        for k in range(nc):
            mass[k] = {k: 'linear'}
            diffusion[k] = {k: {k: 'constant'}}
            potential[k] = {k: {k: 'u'}}
            advection[k] = dict((j, 'nonlinear') for j in range(nc))
        TC_base.__init__(
            self, nc=nc,
            mass=mass,
            advection=advection,
            diffusion=diffusion,
            potential=potential,
            reaction={},
            hamiltonian={},
            variableNames=['v{0}'.format(i) for i in range(nc)]
            )
        self.nu = nu
        self.nc = nc

    def evaluate(self, t, c):
        for k in range(self.nc):
            vk = c[('u', k)]
            c[('m', k)][:] = vk
            c[('dm', k, k)][:] = 1
            c[('a', k, k)][..., 0, 0] = self.nu

            c[('f', k)][..., 0] = 0.0
            for i in range(self.nc):
                vi = c[('u', i)]
                for j in range(self.nc):
                    vj = c[('u', j)]
                    c[('f', k)][..., 0] = c[('f', k)][..., 0] + \
                        0.5 * vi * vj * e_div_gamma[i, j, k]
                    c[('df', k, i)][..., 0] = vj * e_div_gamma[i, j, k]


coefficients = gPCBurgersEqn(nu=1e-6, nc=nc)


def emptyBC(x):
    return None

dirichletConditions = dict((i, emptyBC) for i in range(nc))
advectiveFluxBoundaryconditions = dict((i, emptyBC) for i in range(nc))
diffusiveFluxBoundaryConditions = dict((i, {}) for i in range(nc))


def getPBC(x):
    eps = 1.0e-8
    if x[0] < eps or x[0] >= L[0] - eps:
        return np.array([0.0, 0.0, 0.0])

periodicDirichletConditions = dict((i, getPBC) for i in range(nc))


class vkIC:
    def __init__(self, k):
        self.k = k

    def uOfXT(self, x, t):
        if self.k is 0:
            return np.sin(x[0])
        elif self.k is 1:
            return 1.0
        else:
            return 0.0

initialConditions = dict((k, vkIC(k=k)) for k in range(nc))
fluxBoundaryConditions = dict((i, 'outFlow') for i in range(nc))

T = 1.0
