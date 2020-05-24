
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


## ELASTO-VISCO-PLASTIC ##
function assemble_cell!(Ke, re, model::Model{dim,1,Nothing,V,E,P}, cell, cv, n_basefuncs, ue, noplast = false) where {dim,V,E,P}
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
        @timeit "tangent computation" σ, D = compute_stress_tangent(ϵ, r, states[q_point], model.clock, noplast = noplast)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs #TODO change loops order
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))
            δu = shape_value(cv,q_point,i) # keep it 2D to fit body_forces dim
            re[i] += ((δϵ ⊡ σ) - (δu⋅bodyforces)) * dΩ # add internal force and external body forces to residual. ⊡ : double contraction, ⋅ : single contraction
            for j in 1:i
                Δϵ2D = shape_symmetric_gradient(cv, q_point, j)
                Δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,Δϵ2D))
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke) #when j in 1:i

    # Add traction as a negative contribution to the element residual `re`:
    apply_Neumann_bc!(re, model, cell, n_basefuncs)
end
