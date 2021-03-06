

# function get_variables_interpolation(variables::Variables{NV}, el_geometry::Type{Cell{dim,N,M}}) where {L,dim,N,M,var1,var2,var3}
#     vars = variables.names
#     (var1 == :u) && (interp_u = Lagrange{dim,refshape(el_geometry),var_interp_order[1]}())
#     (var2 == :p) && (L == 2) && (interp_p = Lagrange{dim,refshape(el_geometry),var_interp_order[2]}())
#     (L == 2) && return interp_u, interp_p
#     (L == 1) && return interp_u
#     @error "more than 2 variables (:u,:p) is not implemented yet"
#
# end

function get_quadrature_rules(order::Int, type::Symbol, el_geom::Type{Cell{dim,N,M}}) where{dim,N,M}
    quad_cell = QuadratureRule{dim,refshape(el_geom)}(type, order)
    quad_face = QuadratureRule{dim-1,refshape(el_geom)}(type, order)
    return quad_cell, quad_face
end

function create_values(qr, qr_face, interp_geom, interp_u)
        # cell and facevalues for u
        cellvalues_u = CellVectorValues(qr, interp_u, interp_geom)
        facevalues_u = FaceVectorValues(qr_face, interp_u, interp_geom)

    return (cellvalues_u,), facevalues_u
end;

function create_values(qr, qr_face, interp_geom, interp_u, interp_p)
    # cell and facevalues for u
    cellvalues_u = CellVectorValues(qr, interp_u, interp_geom)
    facevalues_u = FaceVectorValues(qr_face, interp_u, interp_geom)
    # cellvalues for p
    cellvalues_p = CellScalarValues(qr, interp_p, interp_geom)

    return (cellvalues_u, cellvalues_p), facevalues_u
end;

get_dim(c::Cell{dim,N,M}) where{dim,N,M} = dim;
get_dim(g::Grid) = get_dim(g.cells[1]);

function create_dofhandler(grid::Grid, variables::PrimitiveVariables)
    dh = DofHandler(grid)
    for (i,varname) in enumerate(variables.names)
        if varname == :u
            push!(dh, varname, get_dim(grid), variables.interpolations[i])
        else
            push!(dh, varname, 1, variables.interpolations[i])
        end
    end
    close!(dh)
    return dh
end;

function getset(grid::Grid,name::String)
    for (iset,set) in enumerate((grid.cellsets, grid.facesets, grid.nodesets))
        for key in keys(set)
            if key == name
                iset == 1 && return getcellset(grid,name)
                iset == 2 && return getfaceset(grid,name)
                iset == 3 && return getnodeset(grid,name)
            end
        end
    end
    @warn "set name \"$name\" not found in grid object"
end

function get_set_type(grid::Grid,name::String)
    for (iset,set) in enumerate((grid.cellsets, grid.facesets, grid.nodesets))
        for key in keys(set)
            if key == name
                iset == 1 && return :cellset
                iset == 2 && return :faceset
                iset == 3 && return :nodeset
            end
        end
    end
    @warn "set name \"$name\" not found in grid object"
end

function create_dirichlet_bc(dh::DofHandler, bc_dirichlet::Dict)
    dbc = ConstraintHandler(dh)
    for (key,value) in pairs(bc_dirichlet)
        add!(dbc, Dirichlet(value[1], getset(dh.grid, key), value[2], value[3]))
    end
    close!(dbc)
    t = 0.0
    update!(dbc, t)
    return dbc
end

function modify_dirichlet_bc(dh::DofHandler, ch::ConstraintHandler)
    dbc0_vec = ch.dbcs
    dbc = ConstraintHandler(dh)
    for dbc0 in dbc0_vec
        add!(dbc, Dirichlet(dbc0.field_name, dbc0.faces, dbc0.f, dbc0.components))
    end
    close!(dbc)
    t = 0.0
    update!(dbc, t)
    return dbc
end

function setup_model(grid::Grid, variables::PrimitiveVariables{NV},
                     quad_order::Int, quad_type::Symbol,
                     bc_dicts::BoundaryConditions, body_forces::BodyForces, initial_state, rheology::Rheology{T,TD,TV,TE,TP}) where {NV,T,TD,TV,TE,TP}

    # load packages
    # if (NV == 2) & (TD <:Damage)
    #     @info("loading NLsolve package")
    #     @eval(Rheologies, import NLsolve)
    #     @info("loading LineSearches package")
    #     @eval(Rheologies, import LineSearches)
    # end

    # get elements geometry
    el_geom = getcelltype(grid)

    # By default geometry interpolation is quadratic only if there are nodes in edges centers.
    interp_geom = JuAFEM.default_interpolation(el_geom) #

    # quadrature rules
    qr, qr_face = get_quadrature_rules(quad_order, quad_type, el_geom)

    # Dof handler and setup dirichlet bc
    dh = create_dofhandler(grid, variables)
    bcd = create_dirichlet_bc(dh, bc_dicts.dirichlet)

    # cellvalues
    cellvalues, facevalues = create_values(qr, qr_face, interp_geom, variables.interpolations...)

    # sparsity pattern
    K = create_sparsity_pattern(dh);
    RHS = zeros(ndofs(dh))

    mp = create_material_properties(grid, rheology)
    ms = create_material_states(mp,variables,grid,cellvalues[1],initial_state)
    bf = create_body_forces_field(grid, body_forces)

    return dh, bcd, cellvalues, facevalues, mp, ms, bf, K, RHS
end

function get_face_coordinates(cell::CellIterator, face::Int)
    face_nodes_id = J.faces(cell.grid.cells[cell.current_cellid.x])[face]
    return [cell.grid.nodes[nodeid].x for nodeid in face_nodes_id]
end

function get_cell_centroid(grid,cellid)
    nodesid = Base.vect(grid.cells[cellid].nodes...)
    nodes = StructArray(grid.nodes[nodesid])
    return mean(nodes.x)
end


create_material_properties(grid::Grid, rheology::Rheology{T}) where {T<:AbstractFloat} = fill(rheology,getncells(grid))

function create_material_properties(grid::Grid{dim}, rheology::Rheology{Function,D,V,E,P}) where {dim,D,V,E,P}
    x_test = rand(dim)
    tparams = gettypeparameters(rheology(x_test))
    mp = Rheology{tparams...}[]
    for cellid in 1:getncells(grid)
        centroid_x = get_cell_centroid(grid,cellid)
        push!(mp,rheology(centroid_x))
    end
    return mp
end

create_body_forces_field(grid::Grid, body_forces::BodyForces{T}) where {T<:Tensor} = fill(body_forces,getncells(grid))

function create_body_forces_field(grid::Grid{dim}, body_forces::BodyForces{T}) where {dim,T<:Function}
    bf = BodyForces{Float64}[]
    for cellid in 1:getncells(grid)
        centroid_x = get_cell_centroid(grid,cellid)
        push!(bf,body_forces.components(centroid_x))
    end
    return bf
end

function convert_to_undamaged_material_properties(grid,mp0)
    r0 = mp0[1]
    TT,_,TV,TE,TP = gettypeparameters(r0)
    mp = Rheology{TT,Nothing,TV,TE,TP}[]
    if TP <: Plasticity
        @inbounds for cellid in 1:getncells(grid)
            r0_cell = mp0[cellid]
            r_cell = Rheology(nothing, r0_cell.viscosity, r0_cell.elasticity, r0_cell.plasticity)
            push!(mp,r_cell)
        end
    elseif TP == Nothing # add a plasticity field to specialize into a method sharing the same state type as the damaged model
        plas = DruckerPrager(μ = 0.6) #
        @inbounds for cellid in 1:getncells(grid)
            r0_cell = mp0[cellid]
            r_cell = Rheology(nothing, r0_cell.viscosity, r0_cell.elasticity, plas)
            push!(mp,r_cell)
        end
    end
    return mp
end

ms_type(r::Rheology{T,Nothing,V,E,Nothing},vars::PrimitiveVariables{N}) where {T,V,E<:Elasticity,N} = BasicMaterialState()
ms_type(r::Rheology{T,Nothing,V,E,P},vars::PrimitiveVariables{N}) where {T,V,E<:Elasticity, P<:Plasticity,N} = PlasticMaterialState()
ms_type(r::Rheology{T,D,V,E,P},vars::PrimitiveVariables{1}) where {T,D<:Damage,V,E<:Elasticity,P} = DamagedPlasticMaterialState(r)
ms_type(r::Rheology{T,D,V,E,P},vars::PrimitiveVariables{2}) where {T,D<:Damage,V,E<:Elasticity,P} = PlasticMaterialState()
# convenience function :
ms_type(mp::Vector,vars) = ms_type(mp[1],vars)

function create_material_states(mp,vars,grid::Grid,cv,::Nothing)
    nqp = getnquadpoints(cv)
    return [[ms_type(mp,vars) for _ in 1:nqp] for _ in 1:getncells(grid)]
end

function create_material_states(mp,vars::PrimitiveVariables{N},grid::Grid,cv,initial_state::Dict) where {N}
    nqp = getnquadpoints(cv)
    ms = Vector{Vector{typeof(ms_type(mp,vars))}}()
    for cellid in 1:getncells(grid)
        r = mp[cellid]
        cell_ms = create_cell_material_states(r,vars,grid,initial_state,cellid,nqp)
        push!(ms,cell_ms)
    end
    return ms
end

function create_cell_material_states(r,vars::PrimitiveVariables{N},grid,initial_state,cellid,nqp) where {N}
    cell_ms = Vector{typeof(ms_type(r,vars))}(undef,nqp)
    for qp in 1:nqp
        cell_ms[qp] = create_qp_material_state(r,N,grid,initial_state,cellid,qp)
    end
    return cell_ms
end

function create_qp_material_state(r::Rheology{T,Nothing,V,E,Nothing},n_prim_vars,grid,initial_state,cellid,qp) where {T,V<:Viscosity,E<:Elasticity}
    for (key, value) in pairs(initial_state)
        if key in (:ϵ,:σ)
            if typeof(value) <: AbstractArray
                (key == :ϵ) && (ϵ = value[cellid][qp])
                (key == :σ) && (σ = value[cellid][qp])
            elseif typeof(value) <: Function
                centroid_x = get_cell_centroid(grid,cellid)
                (key == :ϵ) && (ϵ = value(centroid_x))
                (key == :σ) && (σ = value(centroid_x))
            end
        else
            @warn "initial state key $key is skipped. Usable keys for $(rheology_summary(r)) rheology are :σ and :ϵ"
        end
    end
    !@isdefined(ϵ) && (ϵ = zero(SymmetricTensor{2, 3}))
    !@isdefined(σ) && (σ = zero(SymmetricTensor{2, 3}))
    return BasicMaterialState(ϵ,σ,ϵ,σ)
end

function create_qp_material_state(r::Rheology{T,Nothing,V,E,P},n_prim_vars,grid,initial_state,cellid,qp) where {T,V,E<:Elasticity, P<:Plasticity}
    for (key, value) in pairs(initial_state)
        if key in (:ϵᵖ,:ϵ̅ᵖ,:ϵ,:σ)
            if typeof(value) <: AbstractArray
                (key == :ϵᵖ) && (ϵᵖ = value[cellid][qp])
                (key == :ϵ̅ᵖ) && (ϵ̅ᵖ = value[cellid][qp])
                (key == :ϵ) && (ϵ = value[cellid][qp])
                (key == :σ) && (σ = value[cellid][qp])
            elseif typeof(value) <: Function
                centroid_x = get_cell_centroid(grid,cellid)
                (key == :ϵᵖ) && (ϵᵖ = value(centroid_x))
                (key == :ϵ̅ᵖ) && (ϵ̅ᵖ = value(centroid_x))
                (key == :ϵ) && (ϵ = value(centroid_x))
                (key == :σ) && (σ = value(centroid_x))
            end
        else
            @warn "initial state key $key is skipped. Usable keys for $(rheology_summary(r)) rheology are :σ and :ϵ"
        end
    end
    !@isdefined(ϵᵖ) && (ϵᵖ = zero(SymmetricTensor{2, 3}))
    !@isdefined(ϵ̅ᵖ) && (ϵ̅ᵖ = 0.0)
    !@isdefined(ϵ) && (ϵ = zero(SymmetricTensor{2, 3}))
    !@isdefined(σ) && (σ = zero(SymmetricTensor{2, 3}))
    return PlasticMaterialState(ϵᵖ,ϵ̅ᵖ,ϵ,σ,ϵᵖ,ϵ̅ᵖ,ϵ,σ)
end

function create_qp_material_state(r::Rheology{T,TD,V,E,P},n_prim_vars,grid,initial_state,cellid,qp) where {T,TD<:Damage,V,E<:Elasticity,P}
    for (key, value) in pairs(initial_state)
        if key in (:D,:ϵᵖ,:ϵ̅ᵖ,:ϵ,:σ)
            if typeof(value) isa AbstractArray
                (key == :D) && (D = value[cellid][qp])
                (key == :ϵᵖ) && (ϵᵖ = value[cellid][qp])
                (key == :ϵ̅ᵖ) && (ϵ̅ᵖ = value[cellid][qp])
                (key == :ϵ) && (ϵ = value[cellid][qp])
                (key == :σ) && (σ = value[cellid][qp])
            elseif typeof(value) isa Function
                centroid_x = get_cell_centroid(grid,cellid)
                (key == :D) && (D = value(centroid_x))
                (key == :ϵᵖ) && (ϵᵖ = value(centroid_x))
                (key == :ϵ̅ᵖ) && (ϵ̅ᵖ = value(centroid_x))
                (key == :ϵ) && (ϵ = value(centroid_x))
                (key == :σ) && (σ = value(centroid_x))
            end
        else
            @warn "initial state key $key is skipped. Usable keys for $(rheology_summary(r)) rheology are :σ and :ϵ"
        end
    end

    !@isdefined(D) && (D = r.damage.D₀ + 1e-9)
    if D <= r.damage.D₀
        @warn "initial damage $D is lower or equal to D₀. Set to (D₀ + 1e-9)"
        D = r.damage.D₀ + 1e-9
    end

    !@isdefined(ϵᵖ) && (ϵᵖ = zero(SymmetricTensor{2, 3}))
    !@isdefined(ϵ̅ᵖ) && (ϵ̅ᵖ = 0.0)
    !@isdefined(ϵ) && (ϵ = zero(SymmetricTensor{2, 3}))
    !@isdefined(σ) && (σ = zero(SymmetricTensor{2, 3}))
    (n_prim_vars == 1) && (return DamagedPlasticMaterialState(D,ϵᵖ,ϵ̅ᵖ,ϵ,σ,D,ϵᵖ,ϵ̅ᵖ,ϵ,σ))
    (n_prim_vars == 2) && (return PlasticMaterialState(ϵᵖ,ϵ̅ᵖ,ϵ,σ,ϵᵖ,ϵ̅ᵖ,ϵ,σ))
end


gettypeparameters(::Rheology{T,D,V,E,P}) where {T,D,V,E,P} = (T,D,V,E,P)

function set_initial_solution_vector!(u,dh,initial_values::Vector)
    u .= initial_values
end

function set_initial_solution_vector!(u,dh,initial_values::Dict)
    vars = keys(initial_values)
    field_names = dh.field_names

    # check that initiated variables are defined
    for var in vars
        if var ∉ field_names
            @error "variable $var defined in initial_values is not a primitive variable of the model.\n
            Please, remove it."
        end
    end

    # check that vars functions of space can be evaluated (enough geometry nodes)
    n_geom_nodes_per_cell = length(dh.grid.cells[1].nodes)
    nbasefuncs = getnbasefunctions.(dh.field_interpolations)
    prescribed_fields = [pf for pf in field_names if pf ∈ vars]
    nbfs = [nbf for (i,nbf) in enumerate(nbasefuncs) if field_names[i] ∈ vars]
    max_basefuncs, id = findmax(nbfs)
    # abort if geometry contains less nodes than one of the initial variables in the dict
    if max_basefuncs > n_geom_nodes_per_cell
        @error "field $(prescribed_fields[id]) interpolation requires more nodes than the geometry definition.\n
        Some fields nodes don't have coordinates to evaluate the function provided. \n
        Please decrease the interpolation order of this field ($(dh.field_interpolations[dh.field_interpolations .== prescribed_fields[id]])), \n
        or increase the geometry's nodes per element"
    end

    # find fields dofs offsets
    fields_dims = dh.field_dims
    fields_el_ndofs = nbasefuncs.*fields_dims
    fields_offsets = [0] #first field is not offset
    for i in 2:length(fields_el_ndofs) # following fields offsets are a cumulative sum of offsets
        push!(fields_offsets,fields_offsets[i-1] + fields_el_ndofs[i-1])
    end
    # find valued variables dofs offsets and elements dofs number:
    vars_offsets = zeros(Int,length(vars))
    vars_el_ndofs = zeros(Int,length(vars))
    vars_dims = zeros(Int,length(vars))
    for (i,var) in enumerate(vars)
        for (j,field) in enumerate(field_names)
            if var == field
                vars_offsets[i] = fields_offsets[j]
                vars_el_ndofs[i] = fields_el_ndofs[j]
                vars_dims[i] = fields_dims[j]
            end
        end
    end
    #loop over cells
    for cell in CellIterator(dh)
        for (i,var) in enumerate(vars)
            var_el_ndofs = vars_el_ndofs[i]
            var_dim = vars_dims[i]

            # account for field dimension
            var_1dim_coords = cell.coords[1:Int(var_el_ndofs/var_dim)]
            var_dofs_coords = repeat(var_1dim_coords, inner=var_dim)

            # compute field value with the coordinates
            var_dofs_values = initial_values[var].(var_dofs_coords)

            var_dofs = cell.celldofs[1+vars_offsets[i]:vars_el_ndofs[i]+vars_offsets[i]]
            u[var_dofs] .= var_dofs_values
        end
    end
    return nothing
end


# "Quadrilateral"
# function get_dofs_coordinates(cell,ndofs,offset,cell_type::Cell{2,4,4})
#     if ndofs == 4 # linear interpolation of the field
#         return cell_coords
#     else
#         @error "no congruence between field and geometry"
#     end
# end
#
# "QuadraticQuadrilateral"
# function get_dofs_coordinates(cell,ndofs,offset,cell_type::Cell{2,9,4})
#     if ndofs == 9 # linear interpolation of the field
#         return cell_coords
#     elseif ndofs == 4
#         return cell_coords
#     else
#         @error "no congruence between field and geometry"
#     end
# end

"""
    get_field_dofs(field::Symbol, model::Model)

Returns a vector of the dofs associated with the primitive variable `field` associated with `model`.
"""
function get_field_dofs(field::Symbol, model)
    dh = model.dofhandler
    field_id = findfirst(J.getfieldnames(dh).==field) # field priority during assembly
    fields_ncelldofs = getnbasefunctions.(model.cellvalues_tuple)
    cell_ndofs = sum(fields_ncelldofs) # total number of dofs per cell
    field_ncelldofs = fields_ncelldofs[field_id] # number of dofs per cell of the field
    fields_offset_incell = Int[1] # cell dofs index at which the fields begins
    for i in 2:J.nfields(dh)
        push!(fields_offset_incell, fields_offset_incell[i-1] + fields_ncelldofs[i-1])
    end
    push!(fields_offset_incell,cell_ndofs+1) # add a length(dofs) + 1 value to allow looping over all dofs

    ncells = getncells(dh.grid)
    all_cell_dofs = dh.cell_dofs
    field_dofs = Int[]
    for cell_id in 1:ncells
        cell_dofs = all_cell_dofs[dh.cell_dofs_offset[cell_id] : dh.cell_dofs_offset[cell_id+1] - 1]
        append!(field_dofs, cell_dofs[fields_offset_incell[field_id]:fields_offset_incell[field_id+1] - 1])
    end
    return unique(field_dofs)
end

function get_boundary_matrix(facesets)
    boundary_set = Set{Tuple{Int,Int}}()
    for (k,v) in pairs(facesets)
        push!(boundary_set, v...)
    end
    boundary = J.boundaries_to_sparse(boundary_set)
    return boundary
end


# import and extend JuAFEM generate_grid function to work with gmsh files
#import JuAFEM: generate_grid
function generate_grid(meshfile::String)

    el_type_str, higher_order_el_mat, meshnodes, facesets, nodesets = GmshParser.extract_mesh(meshfile)
    el_type = eval(Base.Meta.parse(el_type_str))
    #println("typeof(el_type) : ", el_type)
    boundary_matrix = get_boundary_matrix(facesets)

    nodes = [Node((x[4], x[5])) for x in meshnodes]
    elements = [el_type((convert.(Ref(Int),higher_order_el_mat[i,2:end])...,)) for i in 1:size(higher_order_el_mat,1)]

    return Grid(elements, nodes, Dict{String,Set{Int}}(), nodesets, facesets, boundary_matrix)
end
