from proteus import *
from proteus.default_p import *
from proteus.TransportCoefficients import TC_base
import numpy as np

"""
1D, burgers equation with Karhunen-Loeve initial condition
"""

nd = 1
L = [2 * np.pi]

#with open('u0.pickle') as f:
    #[x, y] = pickle.load(f)


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


def forcing_term(x, t):
    m = 5  # xi.shape[0]
    f = 1
    xi = np.array(
        [
            [-1.03279051,  1.12144128],
            [0.63810904,  0.33041598],
            [-0.48385153,  0.89359614],
            [-1.27955579,  2.11849731],
            [-0.45502265,  0.75728598]
        ])

    for k in np.arange(1, m + 1):
        f = f + (-1)**k * \
            (xi[k - 1, 0] * np.sin(2 * k * x[..., 0]) \
            + xi[k - 1, 1] * np.cos(3 * k * x[..., 0])) \
            * exp(-np.sin(2*k*t))
    return f

xi = np.random.randn(5, 2)
coefficients = ForcedBurgersEqn(
    nu=1e-6, rofx=lambda x, t: forcing_term(x, t, xi=xi))


def emptyBC(x, flag):
    return None

dirichletConditions = {0: emptyBC}
advectiveFluxBoundaryconditions = {0: emptyBC}
diffusiveFluxBoundaryconditions = {0: {0: emptyBC}}


def getPBC(x, flag):
    eps = 1.0e-8
    if x[0] < eps or x[0] >= L[0] - eps:
        return np.array([0.0, 0.0, 0.0])

periodicDirichletConditions = {0: getPBC}


class KarhunenLoeveIC:
    def __init__(self, u0):
        self.u0 = u0 

    def uOfXT(self, x, t):
        return self.u0(x[0])

# initialConditions = {0: KarhunenLoeveIC(x, y)}
fluxBoundaryConditions = {0: 'outFlow'}

T = 2.0
