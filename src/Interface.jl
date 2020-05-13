

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

function setup_model(grid::Grid, variables::PrimitiveVariables,
                     quad_order::Int, quad_type::Symbol,
                     bc_dicts::BoundaryConditions, rheology::Rheology)
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
    ms = create_material_states(rheology,grid,cellvalues[1])

    return dh, bcd, cellvalues, facevalues, mp, ms, K, RHS
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

function create_material_states(r::Rheology{T,Nothing,V,E,Nothing},grid::Grid,cv) where {T,V,E}
    nqp = getnquadpoints(cv)
    return [[BasicMaterialState() for _ in 1:nqp] for _ in 1:getncells(grid)]
end
function create_material_states(r::Rheology{T,Nothing,V,E,P},grid::Grid,cv) where {T,V,E,P<:Plasticity}
    nqp = getnquadpoints(cv)
    return [[PlasticMaterialState() for _ in 1:nqp] for _ in 1:getncells(grid)]
end
function create_material_states(r::Rheology{T,D,V,E,P},grid::Grid,cv) where {T,D<:Damage,V,E,P<:Plasticity}
    nqp = getnquadpoints(cv)
    return [[DamagedPlasticMaterialState(r) for _ in 1:nqp] for _ in 1:getncells(grid)]
end

gettypeparameters(::Rheology{T,D,V,E,P}) where {T,D,V,E,P} = (T,D,V,E,P)
