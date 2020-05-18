
abstract type Solver end
abstract type LinearSolver <: Solver end
abstract type NonLinearSolver <: Solver end

abstract type DirectSolver <: LinearSolver end

struct BackslashSolver <: DirectSolver end # single core direct solver
struct MUMPS <: DirectSolver end # multicore sparse direct solver

@kwdef struct ConjugateGradients <: LinearSolver
    package::Symbol
    kwargs::NamedTuple = NamedTuple()
end

@kwdef struct MINRES <: LinearSolver
    package::Symbol
    kwargs::NamedTuple = NamedTuple()
end


@kwdef struct NewtonRaphson <: NonLinearSolver
    atol::Float64 = 1e-4 # absolute residual norm tolerance
    max_iter_number::Int = 10 # maximum non linear iterations
    max_div_iter::Int = 4 # maximum diverging iterations before restart
    autodiff::Bool = false # Has the jacobian to be computed numerically
    iter_on_nonassociated_plasticity::Bool = false
    linear_solver::LinearSolver = BackslashSolver()
end
