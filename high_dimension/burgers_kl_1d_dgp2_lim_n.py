from proteus import *
from proteus.default_n import *
from burgers_kl_1d_p import *
from proteus.Archiver import ArchiveFlags
import numpy as np
archiveFlag = ArchiveFlags.EVERY_USER_STEP

tnList = np.linspace(0, T, 1000)

timeOrder = 3
nStagesTime = timeOrder

DT = None
runCFL = 0.1

limiterType = TimeIntegration.DGlimiterP2Lagrange1d  # None

timeIntegration = SSPRKPIintegration
stepController = Min_dt_RKcontroller
nDTout = 10

femSpaces = {0: DG_AffineQuadraticOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd, 4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1, 4)

nn = 10001
nLevels = 1

subgridError = None
massLumping = False

numericalFluxType = RusanovNumericalFlux_Diagonal

shockCapturing = None
multilevelNonlinearSolver = NLNI

usingSSPRKNewton = True
levelNonlinearSolver = SSPRKNewton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

periodicDirichletConditions = {0: getPBC}
parallelPeriodic = False
parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
