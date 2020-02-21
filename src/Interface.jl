

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

function get_quadrature_rule(order::Int, type::Symbol, el_geom::Type{Cell{dim,N,M}}) where{dim,N,M}
    QuadratureRule{dim,refshape(el_geom)}(type, order)
end
