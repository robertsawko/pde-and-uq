from proteus import *
from proteus.default_n import *
from burgers_gpc_1d_p import *
from proteus.Archiver import ArchiveFlags
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP

timeOrder = 3
nStagesTime = timeOrder

DT = None
runCFL = 0.1

limiterType = TimeIntegration.DGlimiterP2Lagrange1d  # None

timeIntegration = SSPRKPIintegration
# stepController = FixedStep
nSteps = 1000
# dt = T/float(nSteps+1)
stepController = Min_dt_RKcontroller
nDTout = 10
tnList = [float(n)/nSteps*T for n in range(nSteps + 1)]

tnList.insert(1, 1e-6)

femSpaces = dict(
    (i, FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis)
    for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd, 4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1, 4)

nn = 201
nLevels = 1

subgridError = None
massLumping = False

numericalFluxType = RusanovNumericalFlux_Diagonal

shockCapturing = None
multilevelNonlinearSolver = NLNI
# maxNonlinearIts = 100

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
l_atol_res = 0.001*nl_atol_res

conservativeFlux = None

checkMass = True

periodicDirichletConditions = dict((i, getPBC) for i in range(nc))
parallelPeriodic = False
parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
