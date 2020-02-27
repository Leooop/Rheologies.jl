import Statistics.mean

function doassemble(cellvalues_u::CellVectorValues{dim}, cellvalues_p::CellScalarValues{dim},
                    facevalues_u::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid,
                    dh::DofHandler, bcn::Maybe(Dict), rheology::Rheology) where {dim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
    Ke = PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    # assemble_up dispatch on rheology type parameters
    assemble_up!(Ke, fe, assembler, dh::DofHandler,
                cellvalues_u, cellvalues_p, facevalues_u, grid,rheology, bcn)

    return K, f
end;

# TODO : implement elasticity assembly in displacement formulation only :
function doassemble(cellvalues_u::CellVectorValues{dim}, facevalues_u::FaceVectorValues{dim},
                    K::SparseMatrixCSC, grid::Grid,
                    dh::DofHandler, bcn::Maybe(Dict), rheology::Rheology) where {dim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    n = ndofs_per_cell(dh)

    fe = zeros(n) # local force vector
    Ke = zeros(n,n) # local stiffness matrix

    # assemble_up dispatch on rheology type parameters and the number of cellvalues_
    assemble_up!(Ke, fe, assembler, dh::DofHandler,
                cellvalues_u, facevalues_u, grid,rheology, bcn)

    return K, f
end;


function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, cellvalues_p, facevalues_u, grid, rheology::Rheology{T,Nothing,E,Nothing}, bcn) where {dim,T,E}

    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:getnbasefunctions(cellvalues_u)]
    Gmod,Kmod = get_G(rheology.elasticity), get_K(rheology.elasticity)
    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2

    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues_u, cell)
        reinit!(cellvalues_p, cell)

        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
            for i in 1:n_basefuncs_u
                ɛdev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
            end
            dΩ = getdetJdV(cellvalues_u, q_point)
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
                    Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/Kmod * δp * p * dΩ
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

# function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, facevalues_u, grid, rheology::Rheology{T,Nothing,E,Nothing}, bcn) where {dim,T,E}
#
#     # cache ɛdev outside the element routine to avoid some unnecessary allocations
#     ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:getnbasefunctions(cellvalues_u)]
#     Gmod,Kmod = get_G(rheology.elasticity), get_K(rheology.elasticity)
#     n_basefuncs_u = getnbasefunctions(cellvalues_u)
#
#     for cell in CellIterator(dh)
#         fill!(Ke, 0)
#         fill!(fe, 0)
#
#         reinit!(cellvalues_u, cell)
#
#         # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
#         @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
#             dΩ = getdetJdV(cellvalues_u, q_point)
#             for i in 1:n_basefuncs_u
#                 #divδu = shape_divergence(cellvalues_u, q_point, i)
#                 #δu = shape_value(cellvalues_u, q_point, i)
#                 for j in 1:i
#                     Ke[i, j] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ
#                 end
#             end
#
#             for i in 1:n_basefuncs_p
#                 δp = shape_value(cellvalues_p, q_point, i)
#                 for j in 1:n_basefuncs_u
#                     divδu = shape_divergence(cellvalues_u, q_point, j)
#                     Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
#                 end
#                 for j in 1:i
#                     p = shape_value(cellvalues_p, q_point, j)
#                     Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/Kmod * δp * p * dΩ
#                 end
#
#             end
#         end
#
#
#         symmetrize_lower!(Ke)
#
#         # We integrate the Neumann boundary using the facevalues.
#         # We loop over all the faces in the cell, then check if the face
#         # is in our `"traction"` faceset.
#         apply_Neumann_bc!(fe, cell, grid, bcn, facevalues_u, n_basefuncs_u)
#
#         assemble!(assembler, celldofs(cell), fe, Ke)
#     end
# end

# function assemble_up!(Ke, fe, cell, cellvalues_u, cellvalues_p, facevalues_u, grid, rheology::Rheology{nothing,E,nothing}, bc_dicts, ɛdev)
#     Gmod,Kmod = get_G(rheology.elasticity), get_K(rheology.elasticity)
#     n_basefuncs_u = getnbasefunctions(cellvalues_u)
#     n_basefuncs_p = getnbasefunctions(cellvalues_p)
#     u▄, p▄ = 1, 2
#     reinit!(cellvalues_u, cell)
#     reinit!(cellvalues_p, cell)
#
#     # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
#     @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
#         for i in 1:n_basefuncs_u
#             ɛdev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
#         end
#         dΩ = getdetJdV(cellvalues_u, q_point)
#         for i in 1:n_basefuncs_u
#             #divδu = shape_divergence(cellvalues_u, q_point, i)
#             #δu = shape_value(cellvalues_u, q_point, i)
#             for j in 1:i
#                 Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ
#             end
#         end
#
#         for i in 1:n_basefuncs_p
#             δp = shape_value(cellvalues_p, q_point, i)
#             for j in 1:n_basefuncs_u
#                 divδu = shape_divergence(cellvalues_u, q_point, j)
#                 Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
#             end
#             for j in 1:i
#                 p = shape_value(cellvalues_p, q_point, j)
#                 Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/Kmod * δp * p * dΩ
#             end
#
#         end
#     end
#
#
#     symmetrize_lower!(Ke)
#
#     # We integrate the Neumann boundary using the facevalues.
#     # We loop over all the faces in the cell, then check if the face
#     # is in our `"traction"` faceset.
#     apply_Neumann_bc!(fe, cell, grid, bc_dicts, facevalues_u, n_basefuncs_u)
#
# end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;


"""
    apply_Neumann_bc!(fe, cell::CellIterator, grid::Grid, bc_dicts::Dict, facevalues, n_basefuncs)

Apply traction boundary conditions on the element force vector `fe`. This is done by looping over the element faces and checking if it belong to a Neumann set in 'bc_dicts.neumann'. The traction applied on each relevant face is then the average of the traction function of `(x,t)` evaluated at the face nodes.
Time dependency of the traction function is not allowed for now.

"""
function apply_Neumann_bc!(fe, cell::CellIterator, grid::Grid, bcn::Maybe(Dict), facevalues, n_basefuncs)
    t = 0.0 # TODO : make use of simulation time when it will be implemented
    bcn != nothing && @inbounds for face in 1:nfaces(cell)
        if onboundary(cell, face)
            for (neumann_set, traction_func) in pairs(bcn)
                if (cellid(cell), face) ∈ getfaceset(grid, neumann_set)
                    reinit!(facevalues, cell, face)
                    face_nodes_coords = get_face_coordinates(cell::CellIterator, face::Int)
                    tractions = traction_func.(face_nodes_coords,Ref(t))
                    face_traction = mean(hcat(tractions...),dims=2)
                    for q_point in 1:getnquadpoints(facevalues)
                        dΓ = getdetJdV(facevalues, q_point)
                        for i in 1:n_basefuncs
                            δu = shape_value(facevalues, q_point, i)
                            fe[i] += (δu ⋅ face_traction) * dΓ
                        end
                    end
                end
            end
        end
    end
end
