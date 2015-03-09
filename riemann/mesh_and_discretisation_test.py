from proteus.iproteus import *
import numpy as np
Profiling.logLevel = 5
Profiling.verbose = True

mesh_sizes = np.array([10, 50, 100, 200])

import burgers_riem_1d_p as physics
import burgers_riem_1d_dgp1_lim_n as numerics
pList = [physics]
nList = [numerics]
so = default_so
so.sList = [default_s]
opts.cacheArchive = True
so.archiveFlag = Archiver.ArchiveFlags.EVERY_SEQUENCE_STEP

for mesh_size in mesh_sizes:
    so.name = pList[0].name = "burgers_dgp1_{0:03d}-".format(mesh_size)
    nList[0].nn = mesh_size
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns

nList[0].timeOrder = 3
nList[0].nStagesTime = nList[0].timeOrder
nList[0].runCFL = 0.1
nList[0].limiterType = TimeIntegration.DGlimiterP2Lagrange1d
nList[0].femSpaces = {0: FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis}
nList[0].elementQuadrature = Quadrature.SimplexGaussQuadrature(pList[0].nd, 4)
nList[0].elementBoundaryQuadrature =\
    Quadrature.SimplexGaussQuadrature(pList[0].nd-1, 4)
nList[0].numericalFluxType = NumericalFlux.RusanovNumericalFlux_Diagonal
for mesh_size in mesh_sizes:
    so.name = pList[0].name = "burgers_dgp2_{0:03d}-".format(mesh_size)
    nList[0].nn = mesh_size
    ns = NumericalSolution.NS_base(so, pList, nList, so.sList, opts)
    ns.calculateSolution(so.name)
    del ns