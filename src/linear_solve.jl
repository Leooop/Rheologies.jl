### LINEAR SOLVE FUNCTIONS ####
function linear_solve!(u,K,RHS,::BackslashSolver)
    u .= Symmetric(K) \ RHS
    return nothing
end

linear_solve(K,RHS,::BackslashSolver) = Symmetric(K) \ RHS


function linear_solve!(δu,K,res,linear_solver::ConjugateGradient)
    if linear_solver.package == :KrylovMethods
        #!isdefined(KrylovMethods) && @eval(import KrylovMethods)
        δu, flag, relres, iter, resvec = KrylovMethods.cg(Symmetric(K), res, maxIter = linear_solver.max_iter)
        @assert flag == 0
        return nothing
    elseif linear_solver.package == :IterativeSolvers
        #!isdefined(IterativeSolvers) && @eval(import IterativeSolvers)
        IterativeSolvers.cg!(δu, Symmetric(K), res; maxiter = linear_solver.max_iter)
        return nothing
    else
        @error "linear solver package $(linear_solver.package) is not implemented."
    end
end
