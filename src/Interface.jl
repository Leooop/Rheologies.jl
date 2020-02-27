

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

    return cellvalues_u, facevalues_u
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

function create_dofhandler(grid::Grid, variables::Variables)
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

function setup_model(grid::Grid, variables::Variables,
                     quad_order::Int, quad_type::Symbol,
                     bc_dicts::BoundaryConditions)
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

    return dh, bcd, cellvalues, facevalues, K
end

function get_face_coordinates(cell::CellIterator, face::Int)
    face_nodes_id = J.faces(cell.grid.cells[cell.current_cellid.x])[face]
    return [cell.grid.nodes[nodeid].x for nodeid in face_nodes_id]
end
