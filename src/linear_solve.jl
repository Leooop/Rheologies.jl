
### LINEAR SOLVE FUNCTIONS ####
# Direct solvers
function linear_solve!(u,K,RHS,::BackslashSolver)
    u .= K \ RHS
    return nothing
end
linear_solve(K,RHS,::BackslashSolver) = K \ RHS

function linear_solve!(u,K,RHS,::MUMPS)
    u .= MUMPSjInv.solveMUMPS(K, RHS, 1, 1)
    return nothing
end


# indirect solvers :
function linear_solve!(δu,K,res,linear_solver::ConjugateGradients)
    if linear_solver.package == :KrylovMethods
        #!isdefined(KrylovMethods) && @eval(import KrylovMethods)
        δu, flag, relres, iter, resvec = KrylovMethods.cg(K, res, linear_solver.kwargs...)
        @assert flag == 0
        return nothing
    elseif linear_solver.package == :IterativeSolvers
        #!isdefined(IterativeSolvers) && @eval(import IterativeSolvers)
        #Preconditioners.UpdatePreconditioner!(p, K)
        IterativeSolvers.cg!(δu, K, res; linear_solver.kwargs...)
        return nothing
    else
        @error "linear solver package $(linear_solver.package) is not implemented."
    end
end

function linear_solve!(δu,K,res,linear_solver::MINRES)
    if linear_solver.package == :IterativeSolvers
        IterativeSolvers.minres!(δu, K, res; linear_solver.kwargs...)
        return nothing
    else
        @error "linear solver package $(linear_solver.package) is not interfaced."
    end
end
