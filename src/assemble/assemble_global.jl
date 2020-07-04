import Statistics.mean

include("assemble_cell.jl")
# function doassemble(cellvalues_u::CellVectorValues{dim},
#                     cellvalues_p::CellScalarValues{dim},
#                     facevalues_u::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid,
#                     dh::DofHandler, bcn::Maybe(Dict), mp::StructArray, rheology::Rheology) where {dim}
#
#     f = zeros(ndofs(dh))
#     assembler = start_assemble(K, f)
#     nu = getnbasefunctions(cellvalues_u)
#     np = getnbasefunctions(cellvalues_p)
#
#     fe = PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
#     Ke = PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix
#
#     # assemble_up dispatch on rheology type parameters
#     assemble_up!(Ke, fe, assembler, dh::DofHandler,
#                 cellvalues_u, cellvalues_p, facevalues_u, grid, mp, rheology, bcn)
#
#     return K, f
# end;

function doassemble!(model::Model{dim,1,Nothing,Nothing,E,Nothing},nbasefuncs) where {dim,E}
    assembler = start_assemble(model.K, model.RHS)

    # Only one primitive variable here
    n = nbasefuncs[1]
    cv = model.cellvalues_tuple[1]

    # initialize local balance equation terms
    fe = zeros(n) # local force vector
    Ke = zeros(n,n) # local stiffness matrix

    @inbounds for (i,cell) in enumerate(CellIterator(model.dofhandler))
        fill!(Ke, 0)
        fill!(fe, 0)

        @timeit "assemble cell" assemble_cell!(Ke, fe, model, cell, cv, n)

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
end

function doassemble!(model::Model{dim,2,Nothing,Nothing,E,Nothing},nbasefuncs) where {dim,E}
    assembler = start_assemble(model.K, model.RHS)
    n1 = nbasefuncs[1]
    n2 = nbasefuncs[2]
    fe = PseudoBlockArray(zeros(n1 + n2), [n1, n2]) # local force vector
    Ke = PseudoBlockArray(zeros(n1 + n2, n1 + n2), [n1, n2], [n1, n2]) # local stiffness matrix

    # unpack things from model
    dh = model.dofhandler
    cvu, cvp, fvu  = model.cellvalues_tuple[1], model.cellvalues_tuple[2], model.facevalues
    mp = model.material_properties

    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2, dim}) for i in 1:n1]
    u▄, p▄ = 1, 2

    @inbounds for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)

        # assemble_up dispatch on rheology type parameters
        @timeit "assemble cell" assemble_cell!(Ke, fe, model, cell, ɛdev, cvu, cvp, n1, n2, u▄, p▄)
        # Assemble local terms to global ones
        assemble!(assembler, celldofs(cell), fe, Ke)
    end
end

function doassemble!(model::Model{dim,1,D,V,E,P},nbasefuncs, u ; noplast = false) where {dim,D,V,E,P}
    assembler = start_assemble(model.K, model.RHS)

    # Only one primitive variable here
    n = nbasefuncs[1]
    cv = model.cellvalues_tuple[1]

    # initialize local balance equation terms
    re = zeros(n) # local force vector
    Ke = zeros(n,n) # local stiffness matrix

    @inbounds for (i,cell) in enumerate(CellIterator(model.dofhandler))
        fill!(Ke, 0)
        fill!(re, 0)

        eldofs = celldofs(cell)
        ue = u[eldofs]

        #@timeit "assemble cell"
        assemble_cell!(Ke, re, model, cell, cv, n, ue, noplast)

        assemble!(assembler, eldofs, re, Ke)
    end
end

### individual assembly of res function and derivative K used for linesearch
function doassemble_res!(model::Model{dim,1,D,V,E,P},nbasefuncs, u ; noplast = false) where {dim,D,V,E,P}
    assembler = start_assemble(model.K, model.RHS)

    # Only one primitive variable here
    n = nbasefuncs[1]
    cv = model.cellvalues_tuple[1]

    # initialize local balance equation terms
    re = zeros(n) # local force vector
    Ke = zeros(n,n) # local stiffness matrix

    @inbounds for (i,cell) in enumerate(CellIterator(model.dofhandler))
        fill!(Ke, 0)
        fill!(re, 0)

        eldofs = celldofs(cell)
        ue = u[eldofs]

        #@timeit "assemble cell"
        assemble_cell!(Ke, re, model, cell, cv, n, ue, noplast)

        assemble!(assembler, eldofs, re, Ke)
    end
end

doassemble_K!(model::Model{dim,1,D,V,E,P},nbasefuncs, u ; noplast = false) where {dim,D,V,E,P} = doassemble_K!(model.K, u, model, nbasefuncs ; noplast = false)
function doassemble_K!(K, u::AbstractVector, cellvalues::CellVectorValues{dim},
                    facevalues::FaceVectorValues{dim}, grid::Grid,
                    dh::DofHandler,bcn::Maybe(Dict), mp, states, t; noplast = false) where {dim}

    assembler = start_assemble(K)
    # Only one primitive variable here
    n = nbasefuncs[1]
    cv = getnbasefunctions(cellvalues[1])
    Ke = zeros(n,n) # local stiffness matrix

    @inbounds for (i,cell) in enumerate(CellIterator(model.dofhandler))
        fill!(Ke, 0)

        eldofs = celldofs(cell)
        ue = u[eldofs]

        #@timeit "assemble cell"
        assemble_cell!(Ke, re, model, cell, cv, n, ue, noplast)

        assemble!(assembler, eldofs, re, Ke)
    end
end
function doassemble_K(K::SparseMatrixCSC, u::AbstractVector, cellvalues::CellVectorValues{dim},
                    facevalues::FaceVectorValues{dim}, grid::Grid,
                    dh::DofHandler,bcn::Maybe(Dict), mp, states, t ; elastic_only = false) where {dim}
    assembler = start_assemble(K)
    nu = getnbasefunctions(cellvalues)
    ke = zeros(nu, nu) # element tangent matrix

    @inbounds for (i, cell) in enumerate(CellIterator(dh))
        state, r = states[i], mp[i]
        fill!(ke, 0)
        eldofs = celldofs(cell)
        ue = u[eldofs]
        assemble_cell_K!(ke, cell, cellvalues, facevalues, grid, r, bcn,
                       ue, state, t, nu, elastic_only)
        assemble!(assembler, eldofs, ke)
    end
    return K
end

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
function apply_Neumann_bc!(fe, model, cell::CellIterator, n_basefuncs ; inc_sign = +)
    t = (model.clock isa Clock) ? model.clock.current_time : 0.0
    model.neumann_bc != nothing && @inbounds for face in 1:nfaces(cell)
        if onboundary(cell, face)
            for (neumann_set, traction_func) in pairs(model.neumann_bc)
                if (cellid(cell), face) ∈ getfaceset(model.grid, neumann_set)
                    reinit!(model.facevalues, cell, face)
                    face_nodes_coords = get_face_coordinates(cell::CellIterator, face::Int)
                    tractions = traction_func.(face_nodes_coords,Ref(t))
                    face_traction = mean(hcat(tractions...),dims=2)
                    for q_point in 1:getnquadpoints(model.facevalues)
                        dΓ = getdetJdV(model.facevalues, q_point)
                        for i in 1:n_basefuncs
                            δu = shape_value(model.facevalues, q_point, i)
                            (inc_sign == +) && (fe[i] += (δu ⋅ face_traction) * dΓ)
                            (inc_sign == -) && (fe[i] -= (δu ⋅ face_traction) * dΓ)
                        end
                    end
                end
            end
        end
    end
end
