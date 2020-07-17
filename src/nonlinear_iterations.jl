nonlinear_solve!(u, u_prev, δu, model::Model{dim,2,D,V,E,P,S}, restart_flag) where {dim,D,V,E,P,S<:AbstractLinearSolver} = @error "Please choose a non linear solver"

function nonlinear_solve!(u::Vector, u_prev::Vector, δu::Vector, model::Model{DIM,2,TD,TV,TE,TP,TS}, restart_flag) where {DIM,TD,TV,TE,TP,TS<:AbstractNonLinearSolver}

    # Unpack some model fields
    grid, dh, dbc, mp, states, res, K, clock, solver = model.grid, model.dofhandler, model.dirichlet_bc, model.material_properties, model.material_state, model.RHS, model.K, model.clock, model.solver

    println("Δt in nlsolve : ", clock.Δt)

    # number of base functions per element
    nbasefuncs = getnbasefunctions.(model.cellvalues_tuple)


    # flags for restarting nonlinear iterations
    divergence_flag = false # flag when diverging residual norm
    max_iter_flag = false # flag when reaching max number of NL iterations

    # some initializations
    divergence_count = 0 # used to trigger reset of non linear iterations if residual norm diverges
    norm_prev = 0.0

    # apply boundary dirichlet boundary conditions to u :
    update!(dbc, clock.current_time) # evaluates the D-bndc at time t
    apply!(u, dbc)  # set the prescribed values in the solution vector


    newton_itr = -1 # initialize non linear iteration count
    while true; newton_itr += 1
        @timeit "Newton iter" begin
            #@timeit "assemble"
            if newton_itr == 0
                tt = @elapsed doassemble_elast!(res, K, model, nbasefuncs, u, u_prev)
            else
                tt = @elapsed doassemble_res!(res, model, nbasefuncs, u, u_prev)
            end
            println("assemble time : ", tt)
            # compute residual norm
            norm_res = norm(res[JuAFEM.free_dofs(dbc)])

            # print current Newton iteration
            print("\nIteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_res))\n")

            #### Max D TEST :
            if model.material_properties[1].damage isa Damage
                dofs_D = get_field_dofs(:D,model)
                D = u[dofs_D]
                D_prev = u_prev[dofs_D]
                max_ΔD = maximum(D.-D_prev)
                @assert all(D.-D_prev .>= 0)
                print("\t max ΔD = ", max_ΔD, "\n")
            end
            ####

            #### EXIT CHECKS ####
            if (norm_res < solver.atol)
                break
            elseif (solver.max_iter_atol > 0) & (newton_itr == solver.max_iter_number) & (norm_res < solver.max_iter_atol)
                (norm_res >= solver.atol) && @warn("last newton iteration, accepted residual is ",norm_res)
                break # if
            end

            # if too many successive divergence, reset the nonlinear iterations from last time iteration solution u_converged using a lower timestep
            if newton_itr == 0
                norm_prev = norm_res
            else
                if norm_res >= norm_prev
                    divergence_count += 1
                else
                    divergence_count = 0
                end
                norm_prev = norm_res
            end

            divergence_flag = (divergence_count == solver.max_div_iter)
            max_iter_flag = (newton_itr == solver.max_iter_number)

            if divergence_flag | max_iter_flag
                divergence_flag && (@info "$(solver.max_div_iter) successive non linear iterations were diverging. Restart with decreased timestep")
                max_iter_flag   && (@info "reached maximum number of non linear iterations : $(solver.max_iter_number). Restart with decreased timestep")
                restart_flag = true
                return restart_flag
            end

            if newton_itr == 0
                @timeit "apply_dbc" apply_zero!(K, res, dbc)
                @timeit "linear_solve" δu .= solver.linear_solver(K,-res,model)

                ##### TEST #####
                dofs_D = get_field_dofs(:D,model)
                println("δD elast extrema = ", extrema(δu[dofs_D]))
                ################
            else
                dofs_D = get_field_dofs(:D,model)
                ### Compute Jacobian using finite difference ###
                f!(res,u) = doassemble_res!(res, model, nbasefuncs, u, u_prev)
                #f_trailing_args!(res,u,model,nbasefuncs,u_prev) = doassemble_res!(res, model, nbasefuncs, u, u_prev)
                f(u) = (res = similar(u) ; doassemble_res!(res, model, nbasefuncs, u, u_prev) ; return res)
                print(" hand made Jac time : ")
                @time J = compute_jacobian(f!,u,1e-6)
                println("hand made Jacobian L2-norm conditioning = ", cond(Array(J)))

                uu_jac =
                DD_jac = 
                # print(" FiniteDiff Jac time : ")
                # @time begin
                #     # @eval using FiniteDifferences #, FiniteDiff, SparsityDetection, SparseDiffTools
                #     # input = rand(length(u))
                #     # output = similar(input)
                #     # sparsity_pattern = jacobian_sparsity(f!,output,input)
                #     # sparsejac = Float64.(sparse(sparsity_pattern))
                #     # colors = matrix_colors(sparsejac)
                #     # J2 = FiniteDiff.JacobianCache(u,colorvec=colors,sparsity=sparsejac)
                #     # J2 = FiniteDiff.finite_difference_jacobian(f,u)
                #     # fdm = FiniteDifferences.central_fdm(2,1)
                #     # J2 = FiniteDifferences.jacobian(fdm, f, u)[1]
                #
                #     # println("FiniteDiff Jacobian L2-norm conditioning = ", cond(J2))
                #     # println("max difference 2 jacs = ", maximum(abs.(J.-J2)))
                #     # J = sparse(J2)
                # end


                ### Apply zeros accordingly to dbcs ###
                @timeit "apply_dbc" apply_zero!(J, res, dbc)

                println("presence of nans : ", any(isnan.(J)))
                println("maximum((J'.-J)./J) : ", maximum([J[i,j].-J[j,i]./J[i,j] for i in length(res), j in length(res) if !isnan(J[i,j])]))
                println("positive definiteness of sym(J) : ", isposdef(Symmetric(J)))
                #@timeit "linear_solve" δu .= Symmetric(K) \ -res

                ## Preconditioning :
                # using Preconditioners
                # P = DiagonalPreconditioner(J) # similar to the following
                P = Diagonal(J)
                invP = inv(P)
                println("Preconditioned Jacobian L2-norm conditioning = ", cond(Array(invP*J)))
                @timeit "linear_solve" δu .= solver.linear_solver(invP*J,Array(-invP*res),model)
                #@timeit "linear_solve" δu .= solver.linear_solver(J,-res,model)

                ##### TEST #####
                println("δD damaged extrema = ", extrema(δu[dofs_D]))
                # TODO try to set negative δD to zero, remove it after test
                # for i in dofs_D
                #     if δu[i] < 0
                #         δu[i] = 0.0
                #     end
                # end
                # println("δD damaged extrema corrected = ", extrema(δu[dofs_D]))
                ################
            end
            # @timeit "linear_solve" begin
            #     ## TODO better
            #     Pl = Preconditioners.AMGPreconditioner(Symmetric(K))
            #     IterativeSolvers.cg!(δu, Symmetric(K), res, Pl = Pl, solver.linear_solver.kwargs...)
            # end
            # displacement correction
            u .+= δu

            # Make sure that the elastic solve didn't affect Damage
            # if newton_itr == 0
            #     dofs_D = get_field_dofs(:D,model)
            #     u[dofs_D] .= u_prev[dofs_D]
            # end
        end
        ###### TEST
        vtk_grid("TEST_u-D_iter$(newton_itr)", model.dofhandler) do vtkfile
            vtk_point_data(vtkfile, model.dofhandler, u)
        end
        ######
    end
end



function nonlinear_solve!(u::Vector, uprev::Vector, δu::Vector, model::Model{DIM,1,TD,TV,TE,TP,TS}, restart_flag) where {DIM,TD,TV,TE,TP,TS<:AbstractNonLinearSolver}

    # Unpack some model fields
    grid, dh, dbc, mp, states, K, res, clock, solver = model.grid, model.dofhandler, model.dirichlet_bc, model.material_properties, model.material_state, model.K, model.RHS, model.clock, model.solver

    println("Δt in nlsolve : ", clock.Δt)

    # number of base functions per element
    nbasefuncs = getnbasefunctions(model.cellvalues_tuple[1])

    # flags for restarting nonlinear iterations
    divergence_flag = false # flag when diverging residual norm
    max_iter_flag = false # flag when reaching max number of NL iterations

    # some initializations
    divergence_count = 0 # used to trigger reset of non linear iterations if residual norm diverges
    norm_prev = 0.0

    # apply boundary dirichlet boundary conditions to u :
    update!(dbc, clock.current_time) # evaluates the D-bndc at time t
    apply!(u, dbc)  # set the prescribed values in the solution vector


    newton_itr = -1 # initialize non linear iteration count
    while true; newton_itr += 1
        @timeit "Newton iter" begin
            #@timeit "assemble"
            tt = @elapsed doassemble!(model::Model{DIM,1,TD,TV,TE,TP}, nbasefuncs, u ; noplast = (newton_itr == 0))
            println("assemble time : ", tt)
            # compute residual norm
            norm_res = norm(res[JuAFEM.free_dofs(dbc)])

            # print current Newton iteration :
            if newton_itr > 0
                print("\nIteration: $newton_itr \t residual: $(@sprintf("%.8f", norm_res))\n")
            else
                print("First elastic iteration: $newton_itr \t residual: $(@sprintf("%.8f", norm_res))\n")
            end

            #### Max D TEST :
            if TD <: Damage
                max_ΔD = 0
                for cell_states in states
                    cell_temp_D = mean([state.temp_D for state in cell_states])
                    cell_D = mean([state.D for state in cell_states])
                    max_ΔD = max(max_ΔD,cell_temp_D - cell_D)
                end
                print("\t max ΔD = ", max_ΔD, "\n")
            end
            ####

            #### EXIT CHECKS ####
            if (norm_res < solver.atol)
                break
            elseif (solver.max_iter_atol > 0) & (newton_itr == solver.max_iter_number) & (norm_res < solver.max_iter_atol)
                (norm_res >= solver.atol) && @warn("last newton iteration, accepted residual is ",norm_res)
                break # if
            end

            # if too many successive divergence, reset the nonlinear iterations from last time iteration solution u_converged using a lower timestep
            if newton_itr == 0
                norm_prev = norm_res
            else
                if norm_res >= norm_prev
                    divergence_count += 1
                else
                    divergence_count = 0
                end
                norm_prev = norm_res
            end

            divergence_flag = (divergence_count == solver.max_div_iter)
            max_iter_flag = (newton_itr == solver.max_iter_number)

            if divergence_flag | max_iter_flag
                divergence_flag && (@info "$(solver.max_div_iter) successive non linear iterations were diverging. Restart with decreased timestep")
                max_iter_flag   && (@info "reached maximum number of non linear iterations : $(solver.max_iter_number). Restart with decreased timestep")
                restart_flag = true
                return restart_flag
            end


            ### Linear Solve for δu ###
            @timeit "apply_dbc" apply_zero!(K, res, dbc)
            #@timeit "linear_solve" δu .= Symmetric(K) \ -res
            @timeit "linear_solve" δu .= solver.linear_solver(K,-res,model)
            # @timeit "linear_solve" begin
            #     ## TODO better
            #     Pl = Preconditioners.AMGPreconditioner(Symmetric(K))
            #     IterativeSolvers.cg!(δu, Symmetric(K), res, Pl = Pl, solver.linear_solver.kwargs...)
            # end
            # displacement correction
            u .+= δu

        end
    end
end

# function output_values(model::Model{dim,1,D,V,E,P}) where {dim,D,V,E,P<:Plasticity}
#     grid, dh, dbc, mp, states, K, res, solver = model.grid, model.dofhandler, model.dirichlet_bc, model.material_properties, model.material_state, model.K, model.RHS, model.solver
#
#     τ_values = zeros(getncells(grid))
#     F_values = zeros(getncells(grid))
#     p_values = zeros(getncells(grid))
#     σxx_values = zeros(getncells(grid))
#     σyy_values = zeros(getncells(grid))
#     σzz_values = zeros(getncells(grid))
#     ϵxx_values = zeros(getncells(grid))
#     ϵyy_values = zeros(getncells(grid))
#     ϵzz_values = zeros(getncells(grid))
#     ϵkk_values = zeros(getncells(grid))
#     γ_values = zeros(getncells(grid))
#     norm_ep_values = zeros(getncells(grid))
#     I2D = SymmetricTensor{2,3}(δ) # second order identity tensor
#     for (el, cell_states) in enumerate(states)
#         for state in cell_states
#             #el == 1 && show(state.σ)
#             p = -vol(state.temp_σ)[1]
#             s = state.temp_σ + p*I2D
#             p_values[el] += p
#             ϕ = mp[el].plasticity.ϕ
#             C = mp[el].plasticity.C
#             τ_yield = sind(ϕ)*p + cosd(ϕ)*C
#             F_values[el] += get_τ(s,mp[el].plasticity) - τ_yield
#             τ_values[el] += get_τ(s,mp[el].plasticity)
#             norm_ep_values[el] += norm(state.temp_ϵᵖ)
#             σxx_values[el] += state.temp_σ[1,1]
#             σyy_values[el] += state.temp_σ[2,2]
#             σzz_values[el] += state.temp_σ[3,3]
#             ϵxx_values[el] += state.temp_ϵ[1,1]
#             ϵyy_values[el] += state.temp_ϵ[2,2]
#             ϵzz_values[el] += state.temp_ϵ[3,3]
#             ϵkk_values[el] += tr(state.temp_ϵ)
#             dev_ϵ_el = dev(state.temp_ϵ)
#             γ_values[el] += sqrt(2*dev_ϵ_el ⊡ dev_ϵ_el)
#         end
#         p_values[el] /= length(cell_states)
#         τ_values[el] /= length(cell_states)
#         norm_ep_values[el] /= length(cell_states)
#         F_values[el] /= length(cell_states)
#         σxx_values[el] /= length(cell_states)
#         σyy_values[el] /= length(cell_states)
#         σzz_values[el] /= length(cell_states)
#         ϵxx_values[el] /= length(cell_states)
#         ϵyy_values[el] /= length(cell_states)
#         ϵzz_values[el] /= length(cell_states)
#         ϵkk_values[el] /= length(cell_states)
#         γ_values[el] /= length(cell_states)
#     end
#     return (p = p_values, τ = τ_values, norm_ep = norm_ep_values, σxx = σxx_values, σyy = σyy_values, σzz = σzz_values, ϵxx = ϵxx_values, ϵyy = ϵyy_values, ϵzz = ϵzz_values, ϵkk = ϵkk_values, γ = γ_values)
# end
#
# function output(vals, u, dh, dbc, clock)
#     vtk_grid("Rheology_$(clock.iter)", dh) do vtkfile
#          vtk_point_data(vtkfile, dh, u) # displacement field
#          vtk_point_data(vtkfile, dbc) # dirichlet boundary conditions
#          vtk_cell_data(vtkfile, vals.p, "pressure")
#          vtk_cell_data(vtkfile, vals.τ, "tau")
#          vtk_cell_data(vtkfile, vals.norm_ep, "norm_ep")
#          vtk_cell_data(vtkfile, vals.σxx, "sigma xx")
#          vtk_cell_data(vtkfile, vals.σyy, "sigma yy")
#          vtk_cell_data(vtkfile, vals.σzz, "sigma zz")
#          vtk_cell_data(vtkfile, vals.ϵxx, "epsilon xx")
#          vtk_cell_data(vtkfile, vals.ϵyy, "epsilon yy")
#          vtk_cell_data(vtkfile, vals.ϵzz, "epsilon zz")
#          vtk_cell_data(vtkfile, vals.ϵkk, "volumetric strain")
#          vtk_cell_data(vtkfile, vals.γ, "dev strain invariant")
#          println("converged newton_itr saved")
#     end
# end
