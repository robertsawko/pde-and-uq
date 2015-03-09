from proteus import *
from proteus.default_p import *
from proteus.TransportCoefficients import TC_base
import numpy as np
"""
1D, burgers equation Riemann ic
"""

## \page Tests Test Problems 
# \ref burgers_riem_1d_p.py "Burgers equation Riemann problem"
#

# 030508 1d scalar nonlinear hyperbolic pde with convex flux
# classical example but low level test problem.

# #\ingroup test
# \file burgers_riem_1d_p.py
# @{
#
# \brief Burgers equation Riemann problem.
#
# The equation is
#
# \f[
# u_j + \pd{u^2/2}{x} = 0
# \f]
#
# \todo finish burgers_riem_1d_p.py doc

nd = 1
L = [1]


class BurgersEqn(TC_base):
    def __init__(self, nu=0.0):
        TC_base.__init__(
            self, nc=1,
            mass={0: {0: 'linear'}},
            advection={0: {0: 'nonlinear'}},
            diffusion={0: {0: {0: 'constant'}}},
            potential={0: {0: 'nonlinear'}},
            reaction={},
            hamiltonian={},
            variableNames=['u']
            )
        self.nu = nu

    def evaluate(self, t, c):
        u = c[('u', 0)]
        c[('m', 0)][:] = u
        c[('dm', 0, 0)][:] = 1
        c[('f', 0)][..., 0] = 0.5*u**2
        c[('df', 0, 0)][..., 0] = u
        c[('a', 0, 0)][..., 0, 0] = self.nu


coefficients = BurgersEqn(nu=1e-6)


def emptyBC(x):
    return None

dirichletConditions = {0: emptyBC}
advectiveFluxBoundaryconditions = {0: emptyBC}
diffusiveFluxBoundaryconditions = {0: {0: emptyBC}}


def getPBC(x):
    eps = 1.0e-8
    if x[0] < eps or x[0] >= L[0] - eps:
        return np.array([0.0, 0.0, 0.0])

periodicDirichletConditions = {0: getPBC}


class RiemIC:
    def __init__(self):
        self.uLeft = 2.0
        self.uRight = 1.0

    def uOfXT(self, x, t):
        if x[0] >= 0.33 and x[0] <= 0.66:
            return self.uLeft
        else:
            return self.uRight


initialConditions = {0: RiemIC()}
fluxBoundaryConditions = {0: 'outFlow'}

T = 2.0
