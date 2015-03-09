from proteus import *
from proteus.default_p import *
import numpy as np
from proteus.TransportCoefficients import TC_base

# from proteus import Domain

# Numer of dimensions
nd = 3

La, Lb = 5.0, 5.0
# Domain lengths
L = (2 * np.pi, 2 * La, 2 * Lb)

# Case name
# name = "kinetic"
# domain = Domain.RectangularDomain(L=(2 * np.pi, 8), x=[0, -4], name=name)


class KineticEqn(TC_base):
    def __init__(self, nu=0.0, sigma=0.05, fofx=lambda x, t: 0):
        nc = 2
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        mass[0] = {0: 'constant'}
        mass[1] = {1: 'linear'}
        reaction[0] = {0: 'constant'}
        reaction[1] = {0: 'linear'}
        # This is telling how many arguments you've got! so essentially the
        # line belows says "hamiltonian for equation 0 is a linear function of
        # gradient of 0th variable and linear function of a gradient of 1st
        # variable."
        hamiltonian[0] = {0: 'linear', 1: 'linear'}
        hamiltonian[1] = {1: 'linear'}
        TC_base.__init__(self,
                         nc=nc,
                         mass=mass,
                         diffusion=diffusion,
                         potential=potential,
                         reaction=reaction,
                         hamiltonian=hamiltonian,
                         variableNames=['I', 'p'])
        self.nu = nu
        self.fofx = fofx
        self.sigma = sigma

    def evaluate(self, t, c):
        I = c[('u', 0)]
        p = c[('u', 1)]
        x = c['x'][..., 0]
        a = c['x'][..., 1] - La
        b = c['x'][..., 2] - Lb
        gradI = c['grad(u)', 0]
        gradp = c['grad(u)', 1]

        c[('m', 0)][:] = 0

        c[('m', 1)][:] = p
        c[('dm', 1, 1)][:] = 1

        c[('H', 0)][:] = gradI[..., 1] - gradp[..., 0]
        # derivatives with respect  to gradient on hamiltonian terms
        c[('dH', 0, 0)][..., 0] = 0
        c[('dH', 0, 0)][..., 1] = 1.0
        c[('dH', 0, 1)][..., 0] = -1.0
        c[('dH', 0, 1)][..., 1] = 0

        c[('H', 1)][:] =\
            a * gradp[..., 0] + self.sigma * self.fofx(x, t) * gradp[..., 1]
        c[('dH', 1, 1)][..., 0] = a
        c[('dH', 1, 1)][..., 1] = self.sigma * self.fofx(x, t)

        c[('r', 0)][:] = 0

        c[('r', 1)][:] = I
        c[('dr', 1, 0)][:] = 1


def sine_forcing(x, t):
    b = x[2] - Lb
    f = b * np.sin(t)
    return f

coefficients = KineticEqn(
    nu=1e-6, sigma=1, fofx=lambda x, t: sine_forcing(x, t))


class pdf_IC:
    def __init__(self):
        self.C = 5.0 / np.pi**2

    def uOfXT(self, x, t):
        a = x[1] - 5
        b = x[2] - 5
        # print a
        return self.C * np.exp(-(a - np.sin(x[0]))**2 / 2) * np.exp(-50.0 * b**2 / np.pi**2)


class I_IC:
    def __init__(self):
        pass
        self.pi2 = np.pi**2
        self.C = -5.0 / self.pi2

    def uOfXT(self, x, t):
        a = x[1] - La
        b = x[2] - Lb
        return self.C * np.cos(x[0])\
            * np.exp(
                - (np.pi**2 * np.sin(x[0])**2
                    - 2 * self.pi2 * np.sin(x[0]) * a
                    + self.pi2 * a**2 + 100 * b**2)
                / (2 * self.pi2))


initialConditions = {
    0: I_IC(),
    1: pdf_IC()}


def emptyBC(x, flag):
    return None


def dirichletZeroBC(x, flag):
    a = x[1] - 5
    b = x[2] - 5
    if np.abs(a) == 5 or np.abs(b) == 5:
        return lambda x, t: 0.0
    else:
        return None

dirichletConditions = dict((i, dirichletZeroBC) for i in range(2))
advectiveFluxBoundaryConditions = dict((i, emptyBC) for i in range(2))
diffusiveFluxBoundaryConditions = dict((i, {}) for i in range(2))


def getPBC(x, flag):
    eps = 1.0e-8
    if x[0] < eps or x[0] >= L[0] - eps:
        return np.array([0.0, x[1], x[2]])

periodicDirichletConditions = dict((i, getPBC) for i in range(2))
# fluxBoundaryConditions = dict((i, 'outFlow') for i in range(2))

T = 1.0
