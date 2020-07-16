
function solve(model::Model; initial_values = nothing, output_writer = nothing, log = false)
    # TimerOutput handling
    reset_timer!()
    (log == false) && disable_timer!()


    # Reinitialize clock and output_writer # NO NEED, FORCE LOADING ALL INPUT FILE BEFORE SOLVING
    #reinit!(model.clock)
    #reinit!(output_writer,model)

    # iterate over time
    u = iterate(model,output_writer,initial_values)

    # print timings of code annotated with @timeit macro
    (log == true) && print_timer(title = "Analysis with $(getncells(model.grid)) elements", linechars = :ascii)
    println(model.material_state[1][1].temp_σ)
    return model, u
end


function iterate(model::Model{2,N,Nothing,Nothing,E,Nothing}, output_writer, initial_values = nothing) where {N,E<:Elasticity}
    @info "Rheology is purely elastic. Displacement field is solved once at the end of the requested time interval"
    c = model.clock
    c.current_time = c.tspan[2]
    # number of base functions per element
    nbasefuncs = getnbasefunctions.(model.cellvalues_tuple)

    # assembly of
    @timeit "assemble" doassemble!(model,nbasefuncs)

    #Apply Dirichlet boundary conditions
    update!(model.dirichlet_bc, c.current_time)
    apply!(model.K, model.RHS, model.dirichlet_bc)

    # Solve
    #u = Symmetric(model.K) \ model.RHS
    @timeit "linear_solve" u = linear_solve(model.K, model.RHS, model.solver)

    # update material states :
    fill_state_from_u!(model,u) # fill temporary stress and strains
    update_material_state!(model) # set temporary values to converged type fields

    # output
    @timeit "export" write_output!(model, u, output_writer)

    return u
end

function iterate(model::Model{2,2,D,V,E,P}, output_writer, initial_values = nothing) where {N,D,V,E,P}

    # Unpack some model fields
    dh, dbc, cv, clock = model.dofhandler, model.dirichlet_bc, model.cellvalues_tuple, model.clock

    # Pre-allocate solution vectors, etc.
    n_dofs = ndofs(dh)  # total number of dofs
    u  = zeros(n_dofs)
    if initial_values != nothing
        get_initial_solution_vector!(u,dh,initial_values)
    end
    u_converged = copy(u) # backup solution vector
    δu = similar(u)

    # Tuple of the number of shape functions per element and per field
    nbasefuncs = getnbasefunctions.(model.cellvalues_tuple)


    while clock.current_time <= clock.tspan[2]
        @timeit "time iteration" begin

            Δt_max_damage = get_damage_constrained_Δt(model,u,0.3)
            println("Δt_max_damage = ",Δt_max_damage)
            timestep!(clock,Δt_max_damage) # update clock

            print("\n TIME ITERATION $(clock.iter)\n",
            " current simulation time = $(clock.current_time):\n",
            " timestep = $(clock.Δt)\n")

            restart_flag = false

            # Apply dirichlet bc and iteratively solve for u :
            # @timeit "nonlinear solve"
            restart_flag = nonlinear_solve!(u,u_converged,δu,model,restart_flag)

            if restart_flag == true
                u .= u_converged
                undo_timestep!(clock)
                clock.Δt *= clock.Δt_fact_down # decreased timestep
            else # converged
                update_material_state!(model) # update converged state values
                @timeit "export" write_output!(model, u, output_writer) # output

                clock.Δt *= clock.Δt_fact_up # increase timestep
                println(clock.Δt)
                u_converged .= u
            end


            # ##### NLsolve VERSION :#####
            # perform_elastic_solve!(u,model)
            #
            # ###### TEST
            # # vtk_grid("TEST_u-D_after_elastic_solve_it_$(clock.iter)", model.dofhandler) do vtkfile
            # #     vtk_point_data(vtkfile, model.dofhandler, u)
            # # end
            # #write_output!(model, u, output_writer) # output
            # #println("temp_σyy = ", model.material_state[15][1].temp_σ[2,2])
            # #println("temp_ϵyy = ", model.material_state[15][1].temp_ϵ[2,2])
            # ######
            #
            # f!(res,u) = doassemble_res!(res, model, nbasefuncs, u, u_converged)
            # tt = @elapsed (sol = NLsolve.nlsolve(f!, u, method = :newton))
            # println("nlsolve time = ", tt)
            #
            # if NLsolve.converged(sol) == false
            #     println(sol)
            #     u .= u_converged
            #     undo_timestep!(clock)
            #     clock.Δt *= clock.Δt_fact_down # decreased timestep
            # else # converged
            #     u .= sol.zero # converged solution
            #     println(sol)
            #     update_material_state!(model) # update converged state values
            #     @timeit "export" write_output!(model, u, output_writer) # output
            #
            #     clock.Δt *= clock.Δt_fact_up # increase timestep
            #     println(clock.Δt)
            #     u_converged .= u
            # end

            ############################


            clock.current_time == clock.tspan[2] && break # end time loop if requested end time is reached

        end
    end
    return u
end

# TODO check the performances penalty magnitude of such an approach : one solution could
# be to build a container containing most EP_model fields and to reuse it for each subsequent elastic solve
"The current strategy is to convert our damaged model into an equivalent undamaged model and to use the
existing linear elastic solves implemented for undamaged rheologies"
function perform_elastic_solve!(u,model)

    tt = @elapsed begin
    # model parameters
    u_interp = typeof(model.dofhandler.field_interpolations[1]).parameters[3]
    el_geom = typeof(model.grid.cells[1])
    variables = PrimitiveVariables{1}((:u,), (u_interp,), el_geom)
    dh = create_dofhandler(model.grid, variables)
    dbc = modify_dirichlet_bc(dh, model.dirichlet_bc)
    mp = convert_to_undamaged_material_properties(model.grid,model.material_properties)
    K = create_sparsity_pattern(dh)
    RHS = zeros(ndofs(dh))

    EP_model = Model( model.grid,
                      dh,
                      dbc,
                      model.neumann_bc,
                      model.body_forces,
                      (model.cellvalues_tuple[1],),
                      model.facevalues,
                      mp,
                      model.material_state,
                      K,
                      RHS,
                      model.clock,
                      model.solver,
                      model.multithreading )
    end
    println("EP_model building = ", tt)


    dbc, K, res, solver, clock = EP_model.dirichlet_bc, EP_model.K, EP_model.RHS, EP_model.solver, EP_model.clock
    nbasefuncs = getnbasefunctions(EP_model.cellvalues_tuple[1])
    δu = similar(res)

    update!(dbc, clock.current_time) #update newly formed dirichlet bc in EP_model
    dofs_u = get_field_dofs(:u, model) # dofs associated with displacement field in u (u-D formulation)

    tt = @elapsed doassemble!(EP_model, nbasefuncs, u[dofs_u] ; noplast = true)
    println("assemble time : ", tt)
    # compute residual norm
    norm_res = norm(res[JuAFEM.free_dofs(dbc)])
    print("First elastic iteration   \t residual: $(@sprintf("%.8f", norm_res))\n")

    ### Linear Solve for δu ###
    apply_zero!(K, res, dbc)
    δu .= solver.linear_solver(K,-res,EP_model)
    # displacement correction

    @assert length(dofs_u) == length(δu)
    u[dofs_u] .+= δu

    # update initial model temporary state
    rcell, rqp = rand(1:getncells(model.grid)), rand(1:getnquadpoints(model.cellvalues_tuple[1]))
    copy_temp_state!(model,EP_model.material_state)
    @assert EP_model.material_state[rcell][rqp].temp_σ == model.material_state[rcell][rqp].temp_σ

    return nothing
end


function iterate(model::Model{2,1,D,V,E,P}, output_writer, initial_values = nothing) where {N,D,V,E,P}

    # Unpack some model fields
    dh, dbc, cv, clock = model.dofhandler, model.dirichlet_bc, model.cellvalues_tuple, model.clock

    # Pre-allocate solution vectors, etc.
    n_dofs = ndofs(dh)  # total number of dofs
    u  = zeros(n_dofs)  # solution vector
    u_converged = copy(u) # backup solution vector
    δu = zeros(n_dofs)  # displacement correction

    while clock.current_time <= clock.tspan[2]
        @timeit "time iteration" begin

            timestep!(clock) # update clock

            print("\n TIME ITERATION $(clock.iter)\n",
            " current simulation time = $(clock.current_time):\n",
            " timestep = $(clock.Δt)\n")

            restart_flag = false

            # Apply dirichlet bc and iteratively solve for u :
            # @timeit "nonlinear solve"
            restart_flag = nonlinear_solve!(u,u_converged,δu,model,restart_flag)

            if restart_flag == true
                u .= u_converged
                undo_timestep!(clock)
                clock.Δt *= clock.Δt_fact_down # decreased timestep
            else # converged
                update_material_state!(model) # update converged state values
                @timeit "export" write_output!(model, u, output_writer) # output

                clock.Δt *= clock.Δt_fact_up # increase timestep
                println(clock.Δt)
                u_converged .= u
            end

            clock.current_time == clock.tspan[2] && break # end time loop if requested end time is reached

        end
    end
    return u
end
# iterate(model::Model{2,1,Nothing,Nothing,E,P}, output_writer) where {E<:Elasticity,P<:Plasticity} = iterate(model::Model{2,1,Nothing,V,E,P}, output_writer) where {V,E,P}
# iterate(model::Model{2,1,D,Nothing,E,Nothing}, output_writer) where {D<:Damage,E<:Elasticity} = iterate(model::Model{2,1,Nothing,V,E,P}, output_writer) where {V,E,P}
# iterate(model::Model{2,1,D,Nothing,E,P}, output_writer) where {D<:Damage,E<:Elasticity,P<:Plasticity} = iterate(model::Model{2,1,Nothing,V,E,P}, output_writer) where {V,E,P}
# iterate(model::Model{2,1,Nothing,V,E,Nothing}, output_writer) where {V<:Viscosity,E<:Elasticity} = iterate(model::Model{2,1,Nothing,V,E,P}, output_writer) where {V,E,P}


function fill_state_from_u!(model,u)
    @inbounds for (i, cell) in enumerate(CellIterator(model.dofhandler))
        cell_states, r = model.material_state[i], model.material_properties[i]
        eldofs = celldofs(cell)
        nu = getnbasefunctions(model.cellvalues_tuple[1])
        ue = u[eldofs][1:nu]
        fill_cell_state!(model,cell,cell_states,r,ue)
    end
end

function fill_cell_state!(model,cell,cell_states,r,ue)
    cvu = model.cellvalues_tuple[1]
    n_basefuncs = getnbasefunctions(cvu)
    reinit!(cvu, cell)
    Dᵉ = r.elasticity.Dᵉ
    @inbounds for q_point in 1:getnquadpoints(cvu)
        state = cell_states[q_point]
        # For each integration point, compute stress and material stiffness
        ϵ2D = function_symmetric_gradient(cvu, q_point, ue)
        # Total strain recomputed each time because of the newton correction
        state.temp_ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        state.temp_σ = Dᵉ ⊡ state.temp_ϵ
    end
end

export_vtk(model::Model{2,Nothing,V,E,P},u,filename) where {V,E,P} = nothing
function export_vtk(model::Model{2,Nothing,V,E,Nothing},u,filename) where {V,E}
    p_values = zeros(getncells(model.grid))
    σxx_values = zeros(getncells(model.grid))
    σyy_values = zeros(getncells(model.grid))
    σzz_values = zeros(getncells(model.grid))
    σxy_values = zeros(getncells(model.grid))
    ϵxx_values = zeros(getncells(model.grid))
    ϵyy_values = zeros(getncells(model.grid))
    ϵzz_values = zeros(getncells(model.grid))
    ϵxy_values = zeros(getncells(model.grid))
    ϵkk_values = zeros(getncells(model.grid))
    for (el, cell_states) in enumerate(model.material_state)
        for state in cell_states
            p_values[el] += - 1/3 * tr(state.σ)
            σxx_values[el] += state.σ[1,1]
            σyy_values[el] += state.σ[2,2]
            σzz_values[el] += state.σ[3,3]
            σxy_values[el] += state.σ[1,2]
            ϵxx_values[el] += state.ϵ[1,1]
            ϵyy_values[el] += state.ϵ[2,2]
            ϵzz_values[el] += state.ϵ[3,3]
            ϵxy_values[el] += state.ϵ[1,2]
            ϵkk_values[el] += tr(state.ϵ)
        end
        nqp = length(cell_states)
        p_values[el] /= nqp
        σxx_values[el] /= nqp
        σyy_values[el] /= nqp
        σzz_values[el] /= nqp
        σxy_values[el] /= nqp
        ϵxx_values[el] /= nqp
        ϵyy_values[el] /= nqp
        ϵzz_values[el] /= nqp
        ϵxy_values[el] /= nqp
        ϵkk_values[el] /= nqp
    end
    vtk_grid(filename, model.dofhandler) do vtkfile
        vtk_point_data(vtkfile, model.dofhandler, u)
        vtk_point_data(vtkfile, model.dirichlet_bc)
        vtk_cell_data(vtkfile, p_values, "pressure")
        vtk_cell_data(vtkfile, σxx_values, "sigma xx")
        vtk_cell_data(vtkfile, σyy_values, "sigma yy")
        vtk_cell_data(vtkfile, σzz_values, "sigma zz")
        vtk_cell_data(vtkfile, σxy_values, "sigma xy")
        vtk_cell_data(vtkfile, ϵxx_values, "epsilon xx")
        vtk_cell_data(vtkfile, ϵyy_values, "epsilon yy")
        vtk_cell_data(vtkfile, ϵzz_values, "epsilon zz")
        vtk_cell_data(vtkfile, ϵxy_values, "epsilon xy")
        vtk_cell_data(vtkfile, ϵkk_values, "volumetric strain")
    end
end
