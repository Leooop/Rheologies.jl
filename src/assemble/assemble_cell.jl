
## ELASTIC u ##
function assemble_cell!(Ke, fe, model::Model{dim,1,Nothing,Nothing,E,Nothing}, cell, cv, nu) where {dim,E}

    reinit!(cv, cell)
    # get stiffness tensor of the cell
    cell_id = cell.current_cellid.x
    Dᵉ = model.material_properties[cell_id].elasticity.Dᵉ
    bodyforces = model.body_forces[cell_id].components

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:nu
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            δu = shape_value(cv, q_point, i)
            fe[i] += (δu ⋅ bodyforces) * dΩ # sign is ok
            for j in 1:i
                Δϵ2D = shape_symmetric_gradient(cv, q_point, j)
                Δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,Δϵ2D))
                Ke[i, j] += δϵ ⊡ Dᵉ ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
    apply_Neumann_bc!(fe, model, cell, nu)
end

## ELASTIC u,p ##
function assemble_cell!(Ke, fe, model::Model{dim,2,Nothing,Nothing,E,Nothing}, cell, ɛdev, cvu, cvp, nu, np, u▄, p▄) where {dim,E}

    reinit!(cvu, cell)
    reinit!(cvp, cell)

    cell_id = cell.current_cellid.x
    r = model.material_properties[cell_id]
    bodyforces = model.body_forces[cell_id].components
    Gmod, Kmod = r.elasticity.G, r.elasticity.K # get elast moduli

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, q_point)

        for i in 1:nu
            ɛdev[i] = dev(shape_symmetric_gradient(cvu, q_point, i))
        end

        for i in 1:nu
            #divδu = shape_divergence(cellvalues_u, q_point, i)
            δu = shape_value(cvu, q_point, i)
            fe[BlockIndex(u▄, i)] += (δu ⋅ bodyforces) * dΩ
            for j in 1:i
                Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ
            end
        end

        for i in 1:np
            δp = shape_value(cvp, q_point, i)
            for j in 1:nu
                divδu = shape_divergence(cvu, q_point, j)
                Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
            end
            for j in 1:i
                p = shape_value(cvp, q_point, j)
                Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/Kmod * δp * p  * dΩ
            end

        end
    end


    symmetrize_lower!(Ke)

    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    apply_Neumann_bc!(fe, model, cell, nu)

end


## DAMAGED-ELASTO-VISCO-PLASTIC ##
# import Base.getindex
# Base.getindex(barray::PseudoBlockArray{Float64,1,Array{Float64,1},Tuple{BlockedUnitRange{Array{Int64,1}}}}, blockindex::BlockIndex{1}) =
# getindex(barray::PseudoBlockArray{T,N,R,BS} where BS<:Tuple{Vararg{AbstractUnitRange{Int64},N}} where R<:AbstractArray{T,N}, blockindex::BlockIndex{N}) where {T, N}

function assemble_res_cell!(re, model::Model{dim,2,TD,Nothing,TE,TP}, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev) where {dim,TD,TE,TP}

    reinit!(cvu, cell)
    reinit!(cvD, cell)

    cell_id = cell.current_cellid.x
    r = model.material_properties[cell_id]
    bodyforces = model.body_forces[cell_id].components


    ue = nodal_vars_el[1:nu]
    De = exp.(nodal_vars_el[nu+1:end])
    De_prev = exp.(nodal_vars_el_prev[nu+1:end])

    fill!(re, 0)
    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point in 1:getnquadpoints(cvu)

        # get qp state
        # sol_type = eltype(nodal_vars)
        # if sol_type <: ForwardDiff.Dual
        #     state_dual = PlasticMaterialState{sol_type}()
        # else
        state = model.material_state[cell_id][q_point]
        # end

        # get Gauss differential term
        dΩ = getdetJdV(cvu, q_point)

        # get strain
        ϵ2D = function_symmetric_gradient(cvu, q_point, ue)
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))

        # get damage variable
        D = function_value(cvD, q_point, De)
        D_prev = function_value(cvD, q_point, De_prev)

        ###### TEST
        if D <= r.damage.D₀
            println("cell_id = ",cell_id)
            println("ϵ = ",ϵ)
            println("D = ",D)
            println("all_D_qp = ", De)
        end
        ######
        if D < 1
            # get stress
            σ = compute_σij(r,D,ϵ)

            #### TEST
            if (cell_id == 1) & (q_point == 1) & !(eltype(σ) <: ForwardDiff.Dual)
                println("σ at cell_id & qp = 1 : ", σ)
                println("σ_prev at cell_id & qp = 1 : ", state.temp_σ)
            end
            if D <= r.damage.D₀
                println("σ = ",σ)
            end
            ####

            # get damage growth rate
            KI = compute_KI(r,σ,D)

            ###### TEST
            if isnan(KI)
                println("cell_id = ",cell_id) # 1
                println("ϵ = ",ϵ)
                println("σ = ",σ) # All NaN !!!
                println("D = ",D) # OK
            end
            ######

            dDdt = compute_subcrit_damage_rate(r, KI, D)

            for i in 1:nu
                #
                δϵ2D = shape_symmetric_gradient(cvu, q_point, i)
                δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
                δu = shape_value(cvu, q_point, i)
                # increment residual
                #re[BlockIndex(u▄, i)] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
                re[i] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
            end

            for i in 1:nD
                #re[BlockIndex(D▄, i)] += ((D - D_prev)/model.clock.Δt - dDdt) * dΩ
                #(D <= D_prev) && (regularization_term = abs(D-D_prev))

                re[i+nu] += ((D - D_prev)/model.clock.Δt - dDdt ) * dΩ #
            end

            # update qp state
            sol_type = eltype(nodal_vars_el)
            if !(sol_type <: ForwardDiff.Dual) # allows AD without updating state
                set_temp_state!(r,model.clock,state,σ,ϵ)
            end

        else # loss of cohesion, drucker-Prager rheology
            rp = Rheology(nothing,nothing,r.elasticity,r.plasticity)

            sol_type = eltype(nodal_vars_el)
            if sol_type <: ForwardDiff.Dual
                state_dual = PlasticMaterialState{sol_type}() # initializes with zero values
                add_state!(state_dual,state) # add current state
                σ, _ = compute_stress_tangent(ϵ, rp, state_dual, model.clock, noplast = false)
            else
                σ, _ = compute_stress_tangent(ϵ, rp, state, model.clock, noplast = false)
            end

            dDdt = 0.0
            for i in 1:nu
                δϵ2D = shape_symmetric_gradient(cvu, q_point, i)
                δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
                δu = shape_value(cvu, q_point, i)
                # increment residual
                #re[BlockIndex(u▄, i)] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
                re[i] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
            end

            for i in 1:nD
                #re[BlockIndex(D▄, i)] += ((D - D_prev)/model.clock.Δt - dDdt) * dΩ
                #(D <= D_prev) && (regularization_term = abs(D-D_prev))

                re[i+nu] += ((D - D_prev)/model.clock.Δt - dDdt ) * dΩ #
            end
        end
    end

    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    apply_Neumann_bc!(re, model, cell, nu ; inc_sign = -)

end

# TODO : put config out of cell loop, use DiffResults to get func eval + jacobian in a single pass.
function assemble_cell_AD!(re, Ke, model::Model{dim,2,TD,Nothing,TE,TP}, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev) where {dim,TD,TE,TP}
    assemble_res_cell!(re, model, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev)
    f!(re,nodal_vars_el) = assemble_res_cell!(re, model, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev)
    re2 = similar(re)
    jacobian_config = ForwardDiff.JacobianConfig(f!, re2, nodal_vars_el)#, chunk = Chunk(nodal_vars))
    ForwardDiff.jacobian!(Ke, f!, re2, nodal_vars_el, jacobian_config, Val{true}())
end

function assemble_cell_elast!(re, Ke, model::Model{DIM,2,TD,V,E,P}, cell, cvu, cvD, nu, nD, u▄, D▄, ue, De, De_prev) where {DIM,TD,V,E,P}

    reinit!(cvu, cell)
    reinit!(cvD, cell)

    cell_id = cell.current_cellid.x
    r = model.material_properties[cell_id]
    bodyforces = model.body_forces[cell_id].components

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point in 1:getnquadpoints(cvu)

        # get qp state
        state = model.material_state[cell_id][q_point]

        # get Gauss differential term
        dΩ = getdetJdV(cvu, q_point)

        # get strain
        ϵ2D = function_symmetric_gradient(cvu, q_point, ue)
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        σ, C = compute_stress_tangent(ϵ, r, state, model.clock, noplast = true)

        for i in 1:nu
            #
            δϵ2D = shape_symmetric_gradient(cvu, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            δu = shape_value(cvu, q_point, i)
            # increment residual
            #re[BlockIndex(u▄, i)] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
            re[i] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
            δϵC = δϵ ⊡ C
            for j in 1:i
                Δϵ2D = shape_symmetric_gradient(cvu, q_point, j)
                Δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,Δϵ2D))
                Ke[BlockIndex((u▄, u▄), (i, j))] += δϵC ⊡ Δϵ * dΩ
            end
        end

        for i in 1:nD
            #re[BlockIndex(D▄, i)] += ((D - D_prev)/model.clock.Δt - dDdt) * dΩ
            re[i+nu] = 1e3 #* dΩ
            Ke[BlockIndex((D▄, D▄), (i, i))] = 1.0 #* dΩ
        end
    end
    symmetrize_lower!(Ke) #when j in 1:i
    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    apply_Neumann_bc!(re, model, cell, nu ; inc_sign = -)

end

function assemble_cell!(Ke, re, model::Model{dim,1,TD,V,E,P}, cell, cv, n_basefuncs, ue, noplast = false) where {dim,TD,V,E,P}
    # unpack
    cell_id = cell.current_cellid.x
    r, states = model.material_properties[cell_id], model.material_state[cell_id]
    bodyforces = model.body_forces[cell_id].components
    #norm_bodyforces = norm(bodyforces)
    #unit_bodyforces = bodyforces./norm_bodyforces

    # reinit cellvalues
    reinit!(cv, cell)

    @inbounds for q_point in 1:getnquadpoints(cv)
        # For each integration point, compute stress and material stiffness
        ϵ2D = function_symmetric_gradient(cv, q_point, ue)
        # Total strain recomputed each time because of the newton correction
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        # @timeit "tangent computation"
        σ, C = compute_stress_tangent(ϵ, r, states[q_point], model.clock, noplast = noplast)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs #TODO change loops order
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            δu = shape_value(cv,q_point,i) # keep it 2D to fit body_forces dim
            re[i] += ((δϵ ⊡ σ) - (δu⋅bodyforces)) * dΩ # add internal force and external body forces to residual. ⊡ : double contraction, ⋅ : single contraction
            δϵC = δϵ ⊡ C
            for j in 1:i
                Δϵ2D = shape_symmetric_gradient(cv, q_point, j)
                Δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,Δϵ2D))
                Ke[i, j] += δϵC ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke) #when j in 1:i

    # Add traction as a negative contribution to the element residual `re`:
    apply_Neumann_bc!(re, model, cell, n_basefuncs ; inc_sign = -)
end

### individual assembly of res function and derivative K used for linesearch

function assemble_cell_res!(re, model::Model{dim,1,TD,V,E,P}, cell, cv, n_basefuncs, ue, noplast = false) where {dim,TD,V,E,P}
    # unpack
    cell_id = cell.current_cellid.x
    r, states = model.material_properties[cell_id], model.material_state[cell_id]
    bodyforces = model.body_forces[cell_id].components
    #norm_bodyforces = norm(bodyforces)
    #unit_bodyforces = bodyforces./norm_bodyforces

    # reinit cellvalues
    reinit!(cv, cell)

    @inbounds for q_point in 1:getnquadpoints(cv)
        # For each integration point, compute stress and material stiffness
        ϵ2D = function_symmetric_gradient(cv, q_point, ue)
        # Total strain recomputed each time because of the newton correction
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        # @timeit "tangent computation"
        σ, C = compute_stress_tangent(ϵ, r, states[q_point], model.clock, noplast = noplast)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs #TODO change loops order
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            δu = shape_value(cv,q_point,i) # keep it 2D to fit body_forces dim
            re[i] += ((δϵ ⊡ σ) - (δu⋅bodyforces)) * dΩ # add internal force and external body forces to residual. ⊡ : double contraction, ⋅ : single contraction
        end
    end

    # Add traction as a negative contribution to the element residual `re`:
    apply_Neumann_bc!(re, model, cell, n_basefuncs)
end

function assemble_cell_K!(Ke, model::Model{dim,1,TD,V,E,P}, cell, cv, n_basefuncs, ue, noplast = false) where {dim,TD,V,E,P}
    # unpack
    cell_id = cell.current_cellid.x
    r, states = model.material_properties[cell_id], model.material_state[cell_id]
    #norm_bodyforces = norm(bodyforces)
    #unit_bodyforces = bodyforces./norm_bodyforces

    # reinit cellvalues
    reinit!(cv, cell)

    @inbounds for q_point in 1:getnquadpoints(cv)
        # For each integration point, compute stress and material stiffness
        ϵ2D = function_symmetric_gradient(cv, q_point, ue)
        # Total strain recomputed each time because of the newton correction
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        # @timeit "tangent computation"
        σ, C = compute_stress_tangent(ϵ, r, states[q_point], model.clock, noplast = noplast)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs #TODO change loops order
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            for j in 1:i
                Δϵ2D = shape_symmetric_gradient(cv, q_point, j)
                Δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,Δϵ2D))
                Ke[i, j] += δϵ ⊡ C ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke) #when j in 1:i
end
