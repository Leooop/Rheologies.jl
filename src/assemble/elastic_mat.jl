# Assembly function for purely elastic materials without time integration
# good function
function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, cellvalues_p, facevalues_u, grid, mp, rheology::Rheology{T,Nothing,Nothing,E,Nothing}, bcn::Maybe(Dict)) where {dim,T,E<:Elasticity{T}}

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

        cell_id = cell.current_cellid.x
        Gmod,Kmod = get_G(mp[cell_id].elasticity), get_K(mp[cell_id].elasticity)
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

# with matricial form :
# function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, cellvalues_p, facevalues_u, grid, mp, rheology::Rheology{T,Nothing,Nothing,E,Nothing}, bcn::Maybe(Dict)) where {dim,T,E<:Elasticity{T}}
#
#     # cache ɛdev outside the element routine to avoid some unnecessary allocations
#     ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:getnbasefunctions(cellvalues_u)]
#     n_basefuncs_u = getnbasefunctions(cellvalues_u)
#     n_basefuncs_p = getnbasefunctions(cellvalues_p)
#     u▄, p▄ = 1, 2
#     m = SVector(1,1,0)
#
#     for cell in CellIterator(dh)
#         fill!(Ke, 0)
#         fill!(fe, 0)
#
#         reinit!(cellvalues_u, cell)
#         reinit!(cellvalues_p, cell)
#
#         cell_id = cell.current_cellid.x
#         D = get_D3(mp[cell_id].elasticity)
#         Gmod,Kmod = get_G(mp[cell_id].elasticity), get_K(mp[cell_id].elasticity)
#         # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
#         @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
#             dΩ = getdetJdV(cellvalues_u, q_point)
#             B = get_B3(cellvalues_u,q_point,n_basefuncs_u)
#             Np = get_N(cellvalues_p,q_point,n_basefuncs_p)
#
#             Ke[Block(u▄, u▄)] .+= B'*D*B .* dΩ
#             Ke[Block(u▄, p▄)] .+= 1/(3*Kmod) .* B'*(D .- 3*Kmod.*I(3))*m*Np .* dΩ
#             Ke[Block(p▄, u▄)] .+= Np'*m'*B .* dΩ
#             Ke[Block(p▄, p▄)] .+= 1/Kmod .*Np'*Np .* dΩ
#         end
#
#         # We integrate the Neumann boundary using the facevalues.
#         # We loop over all the faces in the cell, then check if the face
#         # is in our `"traction"` faceset.
#
#         apply_Neumann_bc!(fe, cell, grid, bcn, facevalues_u, n_basefuncs_u)
#
#         assemble!(assembler, celldofs(cell), fe, Ke)
#     end
# end

function assemble_up_matricial!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues_u, cellvalues_p, facevalues_u, grid, mp, rheology::Rheology{T,Nothing,Nothing,E,Nothing}, bcn::Maybe(Dict)) where {dim,T,E<:Elasticity{T}}

    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:getnbasefunctions(cellvalues_u)]
    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2
    m = SVector(1,1,0)

    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues_u, cell)
        reinit!(cellvalues_p, cell)

        cell_id = cell.current_cellid.x
        Gmod,Kmod = get_G(mp[cell_id].elasticity), get_K(mp[cell_id].elasticity)
        D = [2Gmod  0.0   0.0 ;
              0.0  2Gmod  0.0 ;
              0.0   0.0   Gmod ]
        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
            dΩ = getdetJdV(cellvalues_u, q_point)
            B = get_B3(cellvalues_u,q_point,n_basefuncs_u)
            dev_B = get_dev_B3(B)
            Np = get_N(cellvalues_p,q_point,n_basefuncs_p)

            Ke[Block(u▄, u▄)] .+= dev_B'*D*dev_B .* dΩ
            up_block_increment = - B'*m*Np .* dΩ
            Ke[Block(u▄, p▄)] .+= up_block_increment
            Ke[Block(p▄, u▄)] .+= up_block_increment'
            Ke[Block(p▄, p▄)] .+= -(1/Kmod) .* Np'*Np .* dΩ
        end

        # We integrate the Neumann boundary using the facevalues.
        # We loop over all the faces in the cell, then check if the face
        # is in our `"traction"` faceset.

        apply_Neumann_bc!(fe, cell, grid, bcn, facevalues_u, n_basefuncs_u)

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
end
#
# function get_D3(r)
#     E, ν = r.E, r.ν
#     D = E/((1+ν)*(1-2ν)) .* [1-ν  ν  0.0 ;
#                              ν   1-ν 0.0 ;
#                              0.0 0.0 (1-2ν)/2 ]
# end

function get_B3(cvv::CellVectorValues,q_point,n_basefuncs)
    B = zeros(Float64,(3,n_basefuncs))
    for bfunc in 1:2:n_basefuncs
        dNdx, dNdy = shape_gradient(cvv, q_point, bfunc)[1,:]
        B[1,bfunc] = dNdx
        B[2,bfunc+1] = dNdy
        B[3,bfunc] = dNdy
        B[3,bfunc+1] = dNdx
    end
    return SMatrix{3,n_basefuncs}(B)
end

function get_dev_B3(B)
    #dev(shape_symmetric_gradient(cellvalues_u, q_point, i))
    dev_B = Matrix(B)
    for i in 1:2:size(B,2)
        mean_B = mean((B[1,i],B[2,i+1]))
        dev_B[:,i:i+1] .-= mean_B
    end
    nl,nc = size(dev_B)
    return SMatrix{nl,nc}(dev_B)
end


function get_N(cvv::CellVectorValues,q_point,n_basefuncs)
    N = zeros(Float64,(2,n_basefuncs))
    for bfunc in 1:n_basefuncs
        N[:,bfunc] .= shape_value(cvv, q_point, bfunc)
    end
    return SMatrix{2,n_basefuncs}(N)
end

function get_N(csv::CellScalarValues,q_point,n_basefuncs)
    N = zeros(Float64,n_basefuncs)
    for bfunc in 1:n_basefuncs
        N[bfunc] = shape_value(csv, q_point, bfunc)
    end
    return SMatrix{1,n_basefuncs}(N)
end

function assemble_up!(Ke, fe, assembler, dh::DofHandler{dim}, cellvalues, facevalues, grid, mp, r::Rheology{T,Nothing,Nothing,E,Nothing}, bcn::Maybe(Dict)) where {dim,T,E<:Elasticity{T}}

    n_basefuncs = getnbasefunctions(cellvalues)

    for (i,cell) in enumerate(CellIterator(dh))
        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues, cell)

        cell_id = cell.current_cellid.x
        r = mp[i]
        D = get_elastic_stiffness_duretz(r.elasticity)
        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            # Compute B at current integration point :
            B = get_B(cellvalues,q_point,n_basefuncs)

            Ke += B'*D*B .* dΩ
        end
        Ke = Matrix(Ke)
        #symmetrize_lower!(Ke)

        # We integrate the Neumann boundary using the facevalues.
        # We loop over all the faces in the cell, then check if the face
        # is in our `"traction"` faceset.

        apply_Neumann_bc!(fe, cell, grid, bcn, facevalues, n_basefuncs)

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
end
