from proteus.iproteus import *
import numpy as np
Profiling.logLevel = 5
Profiling.verbose = True

mesh_sizes = np.array([50, 100, 200, 400, 1000])

import burgers_gpc_1d_p as physics
import burgers_gpc_1d_dgp2_lim_n as numerics
pList = [physics]
nList = [numerics]
so = default_so
so.sList = [default_s]
opts.cacheArchive = True
so.archiveFlag = Archiver.ArchiveFlags.EVERY_USER_STEP

for mesh_size in mesh_sizes:
    so.name = pList[0].name = "gpc_dgp2_{0:04d}".format(mesh_size)
    pList[0].coefficients = pList[0].gPCBurgersEqn(nu=1e-4, nc=5)
    nList[0].nn = mesh_size
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns
