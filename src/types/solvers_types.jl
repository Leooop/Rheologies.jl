

abstract type AbstractSolver end
abstract type AbstractLinearSolver <: AbstractSolver end
abstract type AbstractNonLinearSolver <: AbstractSolver end

abstract type AbstractDirectSolver <: AbstractLinearSolver end
abstract type AbstractIterativeSolver <: AbstractLinearSolver end

##### PRECONDITIONERS #####
# abstract type AbstractPreconditioner end
#
# struct ILU <: AbstractPreconditioner
# 	drop_tol::Float64
# end
# (p::ILU)(J) = ilu(J, τ = p.drop_tol)

# Default backslash solver #
struct BackslashSolver <: AbstractDirectSolver end # single core direct solver

(l::BackslashSolver)(K,RHS,model) = K\RHS

# MUMPS multithreaded sparse direct solver #
struct MUMPS <: AbstractDirectSolver end  # TODO

# GMRES from IterativeSolvers package

@kwdef mutable struct GMRESIterativeSolver{T<:Real, Tl, Tr} <: AbstractIterativeSolver
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

function (l::GMRESIterativeSolver{T, Tl, Tr})(J, rhs, model) where {T, Tl, Tr}
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




@kwdef struct NewtonRaphson{T<:AbstractLinearSolver,R1<:Real,R2<:Real} <: AbstractNonLinearSolver
    atol::R1 = 1e-4 # absolute residual norm tolerance
    max_iter_number::Int = 10 # maximum non linear iterations
	max_iter_atol::R2 = 0.0 # absolute residual norm tolerance at max iteration
    max_div_iter::Int = 4 # maximum diverging iterations before restart
    autodiff::Bool = false # Has the jacobian to be computed numerically
    iter_on_nonassociated_plasticity::Bool = false
    linear_solver::T = BackslashSolver()
end

struct NLsolveNewton{LS,R1<:Real,R2<:Real} <: AbstractNonLinearSolver
    xtol::R1
	ftol::R2
	iterations::Int
	store_trace::Bool
	show_trace::Bool
	extended_trace::Bool
	linesearch::LS

	# Custom inner constructor used to import the relevant packages and perform some checks
	function NLsolveNewton{LS,R1,R2}(xtol,ftol,iterations,store_trace,show_trace,extended_trace,linesearch) where {LS,R1<:Real,R2<:Real}
		# load modules
		@info("loading NLsolve package")
        @eval(Rheologies, import NLsolve)
		@info("loading LineSearches package")
        @eval(Rheologies, import LineSearches)

		# Some checks
		@assert xtol >= 0
		@assert ftol >= 0
		return new{LS,R1,R2}(xtol,ftol,iterations,store_trace,show_trace,extended_trace,linesearch)
	end
end

# general outer constructor without type parameters
NLsolveNewton(xtol::R1,ftol::R2,iterations,store_trace,show_trace,extended_trace,linesearch::LS) where {LS,R1<:Real,R2<:Real} = NLsolveNewton{LS,R1,R2}(xtol,ftol,iterations,store_trace,show_trace,extended_trace,linesearch)
# kwargs method
NLsolveNewton(;xtol = 0.0,ftol = 1e-8,iterations = 1000,store_trace = false,show_trace = false,extended_trace = false,linesearch = LineSearches.Static) = NLsolveNewton(xtol,ftol,iterations,store_trace,show_trace,extended_trace,linesearch)
