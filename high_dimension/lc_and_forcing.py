from proteus.iproteus import *
import numpy as np
import pickle
from scipy.interpolate import interp1d
import burgers.physics as ph
import burgers.numerics as nu

Profiling.verbose = True

# lcList = [6, 0.01]
# mList = [10, 200]
lcList = [0.01]
mList = [200]

with open('u0.pickle', 'rb') as f:
    x, y = pickle.load(f)

pList = [ph]
nList = [nu]
so = default_so
so.sList = [default_s]
# opts.cacheArchive = True
opts.logLevel = 1
so.archiveFlag = Archiver.ArchiveFlags.EVERY_USER_STEP

print "Interpolating...",
u0 = interp1d(x, y, kind='nearest')
print "done."

for ind, lc in enumerate(lcList):
    m = mList[ind]
    pList[0].initialConditions \
        = {0: pList[0].KarhunenLoeveIC(u0)}
    pList[0].coefficients = pList[0].ForcedBurgersEqn(
        nu=1e-9,
        sigma=0.0,
        rofx=lambda x, t: 0)
    so.name = pList[0].name = "lc{0}_noforcing".format(lc)
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
    pList[0].coefficients = pList[0].ForcedBurgersEqn(
        nu=1e-9,
        sigma=0.05,
        # rofx=lambda x, t: pList[0].forcing_term(x, t, xi=xi)
        rofx=lambda x, t: pList[0].forcing_term(x, t)
    )
    so.name = pList[0].name = "lc{0}_forcing".format(lc)
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
