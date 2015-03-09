from proteus.iproteus import *
import numpy as np
import pickle
Profiling.logLevel = 5
Profiling.verbose = True

lcList = [6, 0.01]
mList = [10, 200]

import burgers_kl_1d_p as physics
import burgers_kl_1d_dgp2_lim_n as numerics
pList = [physics]
nList = [numerics]
so = default_so
so.sList = [default_s]
opts.cacheArchive = True
so.archiveFlag = Archiver.ArchiveFlags.EVERY_SEQUENCE_STEP

# This was generate with np.random.randn(5, 2) for repeatability
with open('eta.pickle') as f:
    etaf, etaIC = pickle.load(f)

for ind, lc in enumerate(lcList):
    m = mList[ind]
    pList[0].initialConditions \
        = {0: pList[0].KarhunenLoeveIC(m, 25, etaIC[0:m + 1], lc)}
    pList[0].coefficients = pList[0].ForcedBurgersEqn(nu=1e-9)
    so.name = pList[0].name = "lc{0}_noforcing".format(lc)
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
    pList[0].coefficients = pList[0].ForcedBurgersEqn(
        nu=1e-6,
        rofx=lambda x, t: pList[0].forcing_term(x, t, etaf)
    )
    so.name = pList[0].name = "lc{0}_forcing".format(lc)
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
