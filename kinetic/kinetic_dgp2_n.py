from proteus import *
from proteus.default_n import *
from kinetic_p import *
from proteus.Archiver import ArchiveFlags
archiveFlag = ArchiveFlags.EVERY_USER_STEP

import proteus.default_n as n

nSteps = 1000
dt = T/float(nSteps+1)
tnList = [i*dt for i in range(nSteps+1)]

timeOrder = 3
nStagesTime = timeOrder

DT = None
runCFL = 0.1

limiterType = TimeIntegration.DGlimiterP2Lagrange2d  # None

# timeIntegration = TimeIntegration.SSPRKPIintegration
# stepController = StepControl.Min_dt_RKcontroller
timeIntegration = TimeIntegration.BackwardEuler
stepController = StepControl.FixedStep
nDTout = 100

femSpaces = dict(
    (i, FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis)
    for i in range(2))

elementQuadrature = Quadrature.SimplexGaussQuadrature(nd, 4)
elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(nd-1, 4)

nnx = 41
nny = 56
nnz = 2
nLevels = 1

subgridError = None
massLumping = False

numericalFluxType = NumericalFlux.HamiltonJacobi_DiagonalLesaintRaviart

shockCapturing = None
multilevelNonlinearSolver = NLNI

usingSSPRKNewton = False
levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res

conservativeFlux = None

checkMass = True

periodicDirichletConditions = dict((i, getPBC) for i in range(2))
parallelPeriodic = False
parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
