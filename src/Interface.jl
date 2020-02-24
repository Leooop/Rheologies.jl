

refshape(::Type{Cell{1,N,2}}) where{N} = RefCube
refshape(::Type{Cell{2,N,4}}) where{N} = RefCube
refshape(::Type{Cell{3,N,6}}) where{N} = RefCube
refshape(::Type{Cell{2,N,3}}) where{N} = RefTetrahedron
refshape(::Type{Cell{3,N,4}}) where{N} = RefTetrahedron

function get_variables_interpolation(variables::NTuple{L,Symbol}, var_interp_order::NTuple{L,Int}, el_geometry::Type{Cell{dim,N,M}}) where {L,dim,N,M}
    interp_u = 0.0
    (L == 2) && (interp_p = 0.0)
    for (i,var) in enumerate(variables)
        (var == :u) && (interp_u = Lagrange{dim,refshape(el_geometry),var_interp_order[i]}())
        (var == :p) && (interp_p = Lagrange{dim,refshape(el_geometry),var_interp_order[i]}())
    end
    if L == 2
        return interp_u, interp_p
    elseif L == 1
        return interp_u
    else
        @error "more than 2 variables (:u,:p) is not implemented yet"
    end
end

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

    return cellvalues_u, cellvalues_p, facevalues_u
end;

get_dim(c::Cell{dim,N,M}) where{dim,N,M} = dim;
get_dim(g::Grid) = get_dim(g.cells[1]);

function create_dofhandler(grid, ipu::Lagrange)
    dh = DofHandler(grid)
    push!(dh, :u, get_dim(grid), ipu) # displacement
    close!(dh)
    return dh
end;

function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    push!(dh, :u, get_dim(grid), ipu) # displacement
    push!(dh, :p, 1, ipp) # pressure
    close!(dh)
    return dh
end;

create_dofhandler(grid, ip_vars::Tuple) = create_dofhandler(grid, ip_vars...);

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

function setup_model(grid::Grid, variables::Tuple, var_interp_order::Tuple,
                     quad_order::Int, quad_type::Symbol,
                     bc_dicts::BoundaryConditions)
    el_geom = getcelltype(grid)
    interp_geom = default_interpolation(el_geom) # By default the interpolation is quadratic only if there are nodes in edges centers.
    interp_vars = get_variables_interpolation(variables, var_interp_order, el_geom)
    qr, qr_face = get_quadrature_rules(quad_order, quad_type, el_geom)

    dh = create_dofhandler(grid, interp_vars)
    bcd = create_dirichlet_bc(dh, bc_dirichlet)

    # cellvalues
    cellvalues_u, cellvalues_p, facevalues_u = create_values(qr, qr_face, interp_geom, interp_vars...)
    return dh, bcd, cellvalues_u, cellvalues_p, facevalues_u
end
