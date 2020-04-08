function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, facevalues_u, grid, mp, rheology::Rheology{T,Nothing,Nothing,E,DruckerPrager{T}}, bcn) where {dim,T,E<:Elasticity{T}}

    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:getnbasefunctions(cellvalues_u)]
    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2

    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues_u, cell)
        reinit!(cellvalues_p, cell)

        cellid = cell.current_cellid.x
        Gmod,Kmod = get_G(mp[cellid].elasticity), get_K(mp[cellid].elasticity)
        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
            dΩ = getdetJdV(cellvalues_u, q_point)

            for i in 1:n_basefuncs_u
                ɛdev[i] = dev(shape_symmetric_gradient(cellvalues_u, q_point, i))
            end

            for i in 1:n_basefuncs_u
                #divδu = shape_divergence(cellvalues_u, q_point, i)
                #δu = shape_value(cellvalues_u, q_point, i)
                for j in 1:i
                    Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ
                end
            end

            for i in 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, q_point, i)
                for j in 1:n_basefuncs_u
                    divδu = shape_divergence(cellvalues_u, q_point, j)
                    Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
                end
                for j in 1:i
                    p = shape_value(cellvalues_p, q_point, j)
                    Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/Kmod * δp * p  * dΩ
                end

            end
        end


        symmetrize_lower!(Ke)

        # We integrate the Neumann boundary using the facevalues.
        # We loop over all the faces in the cell, then check if the face
        # is in our `"traction"` faceset.
        apply_Neumann_bc!(fe, cell, grid, bcn, facevalues_u, n_basefuncs_u)

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
end
