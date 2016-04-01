import proteus as pr
from proteus.default_n import *
from burgers_riem_1d_p import *
from proteus.Archiver import ArchiveFlags
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
tnList = [0, T]

DT = None
runCFL = 0.3
limiterType = pr.TimeIntegration.DGlimiterP1Lagrange1d

timeIntegration = pr.TimeIntegration.SSPRKPIintegration
stepController = pr.StepControl.Min_dt_RKcontroller
nDTout = 10

femSpaces = {0: pr.FemTools.DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = pr.Quadrature.SimplexGaussQuadrature(nd, 3)
elementBoundaryQuadrature = pr.Quadrature.SimplexGaussQuadrature(nd-1, 3)

nn = 101
subgridError = None
massLumping = False

numericalFluxType = pr.NumericalFlux.RusanovNumericalFlux_Diagonal


shockCapturing = None

multilevelNonlinearSolver = NLNI

levelNonlinearSolver = pr.NonlinearSolvers.SSPRKNewton  # Newton

nonlinearSmoother = None

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = pr.LinearAlgebraTools.SparseMatrix

multilevelLinearSolver = KSP_petsc4py

levelLinearSolver = KSP_petsc4py

linearSmoother = None

linTolFac = 0.001

conservativeFlux = None

periodicDirichletConditions = {0: getPBC}
# parallel settings
parallelPeriodic = True
parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
