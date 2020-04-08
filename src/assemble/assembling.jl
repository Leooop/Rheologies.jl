import Statistics.mean

include("elastic_mat.jl")

function doassemble(cellvalues_u::CellVectorValues{dim},
                    cellvalues_p::CellScalarValues{dim},
                    facevalues_u::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid,
                    dh::DofHandler, bcn::Maybe(Dict), mp::StructArray, rheology::Rheology) where {dim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
    Ke = PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    # assemble_up dispatch on rheology type parameters
    assemble_up!(Ke, fe, assembler, dh::DofHandler,
                cellvalues_u, cellvalues_p, facevalues_u, grid, mp, rheology, bcn)

    return K, f
end;


# TODO : implement elasticity assembly in displacement formulation only :
function doassemble(cellvalues_u::CellVectorValues{dim}, facevalues_u::FaceVectorValues{dim},
                    K::SparseMatrixCSC, grid::Grid{dim},
                    dh::DofHandler, bcn::Maybe(Dict), mp, rheology::Rheology) where {dim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    n = ndofs_per_cell(dh)

    fe = zeros(n) # local force vector
    Ke = zeros(n,n) # local stiffness matrix
    println(typeof(Ke))

    # assemble_up dispatch on rheology type parameters and the number of cellvalues_
    assemble_up!(Ke, fe, assembler, dh::DofHandler,
                cellvalues_u, facevalues_u, grid, mp, rheology, bcn)

    return K, f
end;



function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;


"""
    apply_Neumann_bc!(fe, cell::CellIterator, grid::Grid, bc_dicts::Dict, facevalues, n_basefuncs)

Apply traction boundary conditions on the element force vector `fe`. This is done by looping over the element faces and checking if it belong to a Neumann faceset in 'bc_dicts.neumann'. The traction applied on each relevant face is then computed as the average of the traction function of `(x,t)` evaluated at the face nodes.
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
