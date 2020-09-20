
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
            #divδu = shape_divergence(cvu, q_point, i)
            δu = shape_value(cvu, q_point, i)
            δɛ = shape_symmetric_gradient(cvu, q_point, i)
            fe[i] += (δu ⋅ bodyforces) * dΩ
            for j in 1:i
                # Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ # initial one
                Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * δϵ ⊡ ɛdev[j] * dΩ #TODO check thid or revert change
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


## Functions for viscosity based drucker-prager plasticity :

function assemble_cell!(Ke, fe, model::Model{dim,2,Nothing,TV,TE,TP}, cell, ɛdev, cvu, cvp, nu, np, u▄, p▄, nodal_vars_prev) where {dim,TV,TE,TP<:ViscousDruckerPrager}

    reinit!(cvu, cell)
    reinit!(cvp, cell)

    cell_id = cell.current_cellid.x
    r = model.material_properties[cell_id]
    bodyforces = model.body_forces[cell_id].components
    Gmod, Kmod = r.elasticity.G, r.elasticity.K # get elast moduli
    η = r.viscosity.η
    Δt = model.clock.Δt
    u_prev = nodal_vars_prev
    states = model.material_state[cell_id]

    if r.elasticity.ν != 0.5
        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point in 1:getnquadpoints(cvu)

            state = states[q_point]
            ϵ_prev = function_symmetric_gradient(cvu, q_point,u_prev[1:nu])

            #ϵ_prev = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ_prev_2D))
            dΩ = getdetJdV(cvu, q_point)

            for i in 1:nu
                ɛdev[i] = dev(shape_symmetric_gradient(cvu, q_point, i))
            end

            for i in 1:nu
                #divδu = shape_divergence(cvu, q_point, i)
                δu = shape_value(cvu, q_point, i)
                δN = shape_symmetric_gradient(cvu, q_point, i)
                fe[i] += (ɛdev[i] ⊡ ((2η/Δt)*ϵ_prev) + δu ⋅ bodyforces) * dΩ
                # fe[i] += (δN ⊡ ((2η/Δt)*ϵ_prev) + δu ⋅ bodyforces) * dΩ
                for j in 1:i
                    ϵ = shape_symmetric_gradient(cvu, q_point, j)
                    Ke[BlockIndex((u▄, u▄), (i, j))] +=  ɛdev[i] ⊡ (2Gmod * ɛdev[j] + (2η/Δt) * ϵ) * dΩ # initial one
                    # Ke[BlockIndex((u▄, u▄), (i, j))] +=  δN ⊡ (2Gmod * ɛdev[j] + (2η/Δt) * ϵ) * dΩ
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
    else
        # code Sister pde
        Z = Gmod*Δt/(η + Gmod*Δt)
        @inbounds for q_point in 1:getnquadpoints(cvu)

            state = states[q_point]
            ϵ_prev_2D = function_symmetric_gradient(cvu, q_point,u_prev[1:nu])
            ϵ_prev = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ_prev_2D))
            σ_prev = state.σ

            if cell_id == 1
                println("ϵ_prev : ", ϵ_prev)
                println("σ_prev : ", σ_prev)
            end

            dΩ = getdetJdV(cvu, q_point)

            # ∇N = ɛdev
            # for i in 1:nu
            #     ∇N[i] = shape_symmetric_gradient(cvu, q_point, i)
            # end

            for i in 1:nu
                #divδu = shape_divergence(cvu, q_point, i)
                δu = shape_value(cvu, q_point, i)
                δN_2D = shape_gradient(cvu, q_point, i)
                δN = Tensor{2,3}((i,j)->get_3D_func(i,j,δN_2D))
                fe[i] += (δN ⊡ ( (2Z*η/Δt)*ϵ_prev - (1 - Z)*σ_prev ) + δu ⋅ bodyforces) * dΩ
                # fe[i] += (δN ⊡ ((2η/Δt)*ϵ_prev) + δu ⋅ bodyforces) * dΩ
                for j in 1:i # 1:i initialy with symmetrize_lower
                    ϵ_2D = shape_symmetric_gradient(cvu, q_point, j)
                    ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ_2D))
                    Ke[BlockIndex((u▄, u▄), (i, j))] +=  δN ⊡ ((2Z*η/Δt) * ϵ) * dΩ # initial one
                    # Ke[BlockIndex((u▄, u▄), (i, j))] +=  δN ⊡ (2Gmod * ɛdev[j] + (2η/Δt) * ϵ) * dΩ
                end
            end

            for i in 1:np
                δp = shape_value(cvp, q_point, i)
                for j in 1:nu
                    divδu = shape_divergence(cvu, q_point, j)
                    Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
                end
                # for j in 1:i
                #     p = shape_value(cvp, q_point, j)
                #     Ke[BlockIndex((p▄, p▄), (i, j))] += 0.0
                # end

            end
        end
    end


    symmetrize_lower!(Ke)

    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    apply_Neumann_bc!(fe, model, cell, nu)

end


# TODO : put config out of cell loop, use DiffResults to get func eval + jacobian in a single pass.
function assemble_cell_AD!(re, Ke, model::Model{dim,2,TD,TV,TE,TP}, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev) where {dim,TD,TV<:Viscosity,TE,TP<:ViscousDruckerPrager}
    assemble_res_cell!(re, model, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev)
    f!(re,nodal_vars_el) = assemble_res_cell!(re, model, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev)
    re2 = similar(re)
    jacobian_config = ForwardDiff.JacobianConfig(f!, re2, nodal_vars_el)#, chunk = Chunk(nodal_vars))
    ForwardDiff.jacobian!(Ke, f!, re2, nodal_vars_el, jacobian_config, Val{true}())
end

## DAMAGED-ELASTO-VISCO-PLASTIC ##
# import Base.getindex
# Base.getindex(barray::PseudoBlockArray{Float64,1,Array{Float64,1},Tuple{BlockedUnitRange{Array{Int64,1}}}}, blockindex::BlockIndex{1}) =
# getindex(barray::PseudoBlockArray{T,N,R,BS} where BS<:Tuple{Vararg{AbstractUnitRange{Int64},N}} where R<:AbstractArray{T,N}, blockindex::BlockIndex{N}) where {T, N}
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

        # with homogeneous elasticity
        # r_nodamage = Rheology(nothing,r.viscosity,r.elasticity,r.plasticity) # TODO find a better way to dispatch on compute_stress_tangent
        # σ, C = compute_stress_tangent(ϵ, r_nodamage, state, model.clock, noplast = true)

        # with damaged elasticity :
        D_prev = function_value(cvD, q_point, De_prev)
        C = compute_damaged_stiffness_tensor(r,ϵ,D_prev)
        σ = C ⊡ ϵ

        if (cell_id == 1) & (q_point == 1)
            println("σ = ", σ)
            println("ϵ = ", ϵ)
        end
        # update temp state :
        state.temp_σ = σ
        state.temp_ϵ = ϵ
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
            re[i+nu] = 1e9 #* dΩ
            Ke[BlockIndex((D▄, D▄), (i, i))] = 1.0 #* dΩ
        end
    end
    symmetrize_lower!(Ke) #when j in 1:i
    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    apply_Neumann_bc!(re, model, cell, nu ; inc_sign = -)

end

function assemble_res_cell!(re, model::Model{dim,2,TD,Nothing,TE,TP}, cell, cvu, cvD, nu, nD, u▄, D▄, nodal_vars_el, nodal_vars_el_prev) where {dim,TD,TE,TP}

    reinit!(cvu, cell)
    reinit!(cvD, cell)

    cell_id = cell.current_cellid.x
    r = model.material_properties[cell_id]
    bodyforces = model.body_forces[cell_id].components

    ue = nodal_vars_el[1:nu]
    De = exp.(nodal_vars_el[nu+1:end])
    De_prev = exp.(nodal_vars_el_prev[nu+1:end])
    dDdt = 0.0

    fill!(re, 0)

    if (cell_id == 1) & !(eltype(De) <: ForwardDiff.Dual)
        println("De at cell id 1 : ")
        display(De)
    end
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
        if D <= 1 # TEST TODO put back D < 1
            # get stress
            σ = compute_σij(r,D,ϵ)

            # get damage growth rate
            KI = compute_KI(r,σ,D)

            ###### TEST
            if isnan(KI)
                println("KI is NaN")
                println("\tcell_id = ",cell_id) # 1
                println("\tϵ = ",ϵ)
                println("\tσ = ",σ) # All NaN !!!
                println("\tD = ",D) # OK
            end
            ######

            dDdt = compute_subcrit_damage_rate(r, KI, D)

            #### TEST
            if (cell_id == 1) & (q_point == 1) & !(eltype(σ) <: ForwardDiff.Dual)
                println("σ at cell_id & qp = 1 : ", σ)
                println("σ_prev at cell_id & qp = 1 : ", state.temp_σ)
                println("KI at cell_id & qp = 1 : ", KI)
                println("dDdt at cell_id & qp = 1 : ", dDdt)
            end
            if D <= r.damage.D₀
                println("σ = ",σ)
            end
            ####

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
                δD = shape_value(cvD, q_point, i)

                re[i+nu] += δD * ((D - D_prev)/model.clock.Δt - dDdt ) * dΩ #
            end

            # update qp state
            sol_type = eltype(nodal_vars_el)
            if !(sol_type <: ForwardDiff.Dual) # no state update when AD
                set_temp_state!(r,model.clock,state,σ,ϵ)
            end
        end

        # if (D >= 1) | (D_prev + dDdt*model.clock.Δt >= 1) # loss of cohesion, drucker-Prager rheology
        #     rp = Rheology(nothing,nothing,r.elasticity,r.plasticity)
        #
        #     sol_type = eltype(nodal_vars_el)
        #     if sol_type <: ForwardDiff.Dual
        #         state_dual = PlasticMaterialState{sol_type}() # initializes with zero values
        #         add_state!(state_dual,state) # add current state
        #         σ, _ = compute_stress_tangent(ϵ, rp, state_dual, model.clock, noplast = false)
        #     else
        #         σ, _ = compute_stress_tangent(ϵ, rp, state, model.clock, noplast = false)
        #     end
        #
        #     dDdt = 0.0
        #     for i in 1:nu
        #         δϵ2D = shape_symmetric_gradient(cvu, q_point, i)
        #         δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
        #         δu = shape_value(cvu, q_point, i)
        #         # increment residual
        #         #re[BlockIndex(u▄, i)] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
        #         re[i] += ((δϵ ⊡ σ) - (δu ⋅ bodyforces)) * dΩ
        #     end
        #
        #     for i in 1:nD
        #         #re[BlockIndex(D▄, i)] += ((D - D_prev)/model.clock.Δt - dDdt) * dΩ
        #         #(D <= D_prev) && (regularization_term = abs(D-D_prev))
        #         δD = shape_value(cvD, q_point, i)
        #         re[i+nu] += δD * ((D - D_prev)/model.clock.Δt - dDdt ) * dΩ #
        #     end
        # end
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


    ### TEST ###
    #Ke = compute_jacobian(f!, re2 ; eps = 1e-8)
    # if (cell.current_cellid.x == 1)
    #     println("passé ici !!!")
    #     J = compute_jacobian(f!, nodal_vars_el ; eps = 1e-9)
    #     println("finitediff jacobian for cell_id = 1, iter 2:")
    #     display(J[1:2,:])
    #     println("difference with AD for cell_id = 1, iter 2 :")
    #     display(abs.(Array(Ke)[1:2,:].-J[1:2,:]))
    # end
    ############
    return Ke
end

###################################################
#### DAMAGE with displacement only formulation ####
###################################################

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
