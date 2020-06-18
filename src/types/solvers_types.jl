#import KrylovMethods
import IterativeSolvers
using IncompleteLU
import LinearMaps


abstract type AbstractSolver end
abstract type AbstractLinearSolver <: AbstractSolver end
abstract type AbstractNonLinearSolver <: AbstractSolver end

abstract type AbstractDirectSolver <: AbstractLinearSolver end
abstract type AbstractIterativeSolver <: AbstractLinearSolver end

##### PRECONDITIONERS #####
abstract type AbstractPreconditioner end

struct ILU <: AbstractPreconditioner
	drop_tol::Float64
end
(p::ILU)(J) = ilu(J, τ = p.drop_tol)

# Default backslash solver #
struct BackslashSolver <: AbstractDirectSolver end # single core direct solver

(l::BackslashSolver)(K,RHS,model) = K\RHS

# MUMPS multithreaded sparse direct solver #
struct MUMPS <: AbstractDirectSolver end  # TODO

# GMRES from IterativeSolvers package

@kwdef mutable struct GMRESIterativeSolvers{T, Tl, Tr} <: AbstractIterativeSolver
	"Tolerance for solver"
	tol::T = 1e-4

	"Number of restarts"
	restart::Int64 = 200

	"Maximum number of iterations"
	maxiter::Int64 = 100

	"Display information during iterations"
	verbose::Bool = false

	"Record information"
	log::Bool = true

	"Start with zero guess"
	initially_zero::Bool = true

	"Left preconditioner"
	Pl::Tl = IterativeSolvers.Identity()

	"Right preconditioner"
	Pr::Tr = IterativeSolvers.Identity()
end

function (l::GMRESIterativeSolvers{T, Tl, Tr})(J, rhs, model) where {T, Tl, Tr}
	#J_map = v -> apply(J, v)
	#dim = ndofs(model.dofhandler)
	#Jmap = LinearMaps.LinearMap{T}(J_map, dim, dim ; ismutating = false)
	(l.Pl isa AbstractPreconditioner) && (l.Pl = l.Pl(J))
	res = IterativeSolvers.gmres(J, rhs; tol = l.tol, log = l.log, verbose = l.verbose, restart = l.restart, maxiter = l.maxiter, initially_zero = l.initially_zero, Pl = l.Pl, Pr = l.Pr)
	(res[2].iters >= l.maxiter) && (@warn "IterativeSolvers.gmres iterated maxIter = $(res[2].iters) times without achieving the desired tolerance.\n")
	return res[1]
end




@kwdef struct ConjugateGradients <: AbstractIterativeSolver
    package::Symbol
    kwargs::NamedTuple = NamedTuple()
end

@kwdef struct MINRES <: AbstractIterativeSolver
    package::Symbol
    kwargs::NamedTuple = NamedTuple()
end




@kwdef struct NewtonRaphson{T<:AbstractLinearSolver} <: AbstractNonLinearSolver
    atol::Float64 = 1e-4 # absolute residual norm tolerance
    max_iter_number::Int = 10 # maximum non linear iterations
    max_div_iter::Int = 4 # maximum diverging iterations before restart
    autodiff::Bool = false # Has the jacobian to be computed numerically
    iter_on_nonassociated_plasticity::Bool = false
    linear_solver::T = BackslashSolver()
end
