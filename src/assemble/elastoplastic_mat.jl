###

function assemble_cell!(Ke, re, model::Model{dim,1,Nothing,Nothing,E,P}, cell, cv, n_basefuncs, ue, elastic_only = false) where {dim,E,P}
    # unpack
    cell_id = cell.current_cellid.x
    r, states = model.material_properties[cell_id], model.material_state[cell_id]

    # reinit cellvalues
    reinit!(cv, cell)

    @inbounds for q_point in 1:getnquadpoints(cv)
        # For each integration point, compute stress and material stiffness
        ϵ2D = function_symmetric_gradient(cv, q_point, ue)
        # Total strain recomputed each time because of the newton correction
        ϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,ϵ2D))
        @timeit "tangent computation" σ, D = compute_stress_tangent(ϵ, r, states[q_point], model.clock, elastic_only = elastic_only)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δϵ2D = shape_symmetric_gradient(cv, q_point, i)
            δϵ = SymmetricTensor{2,3}((i,j)->get_3D_func(i,j,δϵ2D))

            re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
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
