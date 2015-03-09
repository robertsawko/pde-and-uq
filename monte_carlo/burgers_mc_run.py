from proteus.iproteus import *
Profiling.logLevel = 5
Profiling.verbose = True
from numpy.random import randn
import numpy as np

import burgers_1d_p as physics
import burgers_1d_dgp2_lim_n as numerics
pList = [physics]
nList = [numerics]
so = default_so
so.sList = [default_s]
opts.cacheArchive = True
so.archiveFlag = Archiver.ArchiveFlags.EVERY_SEQUENCE_STEP

N = 50000


class sineIC:
    def __init__(self, eta):
        self.eta = eta

    def uOfXT(self, x, t):
        return self.eta + np.sin(x[0])


def forcing_term(x, t, xi):
    return xi*np.sin(t)


# p1 basis
for i in range(N):
    so.name = pList[0].name = "mc_run{0:05d}".format(i)
    xi = np.pi / 10 * randn()
    eta = randn()
    pList[0].coefficients = pList[0].ForcedBurgersEqn(
        nu=1e-6,
        rofx=lambda x, t: forcing_term(x, t, xi),
        sigma=1
    )
    pList[0].initialConditions = {0: sineIC(eta)}

    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
