
using Base: @kwdef
import Base: show, getproperty

# aliases
Maybe(T)=Union{T,Nothing}
const F64orFunc = Union{Float64,Function}

### PRIMITIVE VARIABLES ###
@kwdef struct PrimitiveVariables{Nvars}
    names::Union{Symbol,NTuple{Nvars,Symbol}}
    interpolations::Union{Interpolation,NTuple{Nvars,Interpolation}}
    function PrimitiveVariables{Nvars}(names,ips) where {Nvars}
        (typeof(names) == Symbol) && (@assert Nvars == 1)
        (typeof(names) == Tuple) && (@assert Nvars == length(names) == length(ips))
        return new(names,ips)
    end
end

refshape(::Type{Cell{1,N,2}}) where{N} = RefCube
refshape(::Type{Cell{2,N,4}}) where{N} = RefCube
refshape(::Type{Cell{3,N,6}}) where{N} = RefCube
refshape(::Type{Cell{2,N,3}}) where{N} = RefTetrahedron
refshape(::Type{Cell{3,N,4}}) where{N} = RefTetrahedron

function PrimitiveVariables{NV}(names,ip_order,el_geom::Type{Cell{dim,N,M}}) where {NV,dim,N,M}
    NV == 1 && return PrimitiveVariables{NV}(names, (Lagrange{dim,refshape(el_geom),ip_order[1]}(),))
    NV == 2 && return PrimitiveVariables{NV}(names, (Lagrange{dim,refshape(el_geom),ip_order[1]}(),
                                          Lagrange{dim,refshape(el_geom),ip_order[2]}()) )
    @error "More than two primitive variables (:u and :p) not implemented yet"
end

# function Variables(names::Symbol,ips::Interpolation)
#     Variables{1}(names,ips)
# end
# function Variables(names::NTuple{N,Symbol},ips::NTuple{N,Interpolation}) where {N}
#     Variables{N}(names,ips)
# end

### BOUNDARY CONDITIONS ###

@kwdef struct BoundaryConditions
    dirichlet::Dict
    neumann::Maybe(Dict)
end
BoundaryConditions(::Nothing,::Dict) = @error "A minimal set of Dirichlet boundary conditions must be provided in order to prevent rigid motion of the material"

### BODY FORCES ###
struct BodyForces{T}
    components::T
    function BodyForces(c::T) where {T}
        if T <: AbstractVector
            dim = length(c)
            if T <: Vector
                c = Vec{dim}(c)
                new{typeof(c)}(c)
            elseif T <: Tensor
                return new{T}(c)
            end
        elseif T <: Function
            return new{T}(c)
        else
            @error "wrong field type. Must be an AbstractVector"
        end
    end
end

# other types :

include("material_behaviors_types.jl")
include("rheology_type.jl")
include("clock_type.jl")
include("solvers_types.jl")
include("material_state_types.jl")
include("output_writer_types.jl")
include("model_type.jl")



## TODO : Check whether a material parameters type is needed
### MATERIAL PARAMETERS ###

# struct MaterialProperties{T}
#     mp :: Vector{Rheology}
#     function MaterialProperties(mp)
#         T = eltype(mp)
#         new{T}(mp)
#     end
# end
# MaterialProperties{T}(mp) where T = MaterialProperties(mp)
# function Base.getproperty(MP::MaterialProperties,sym::Symbol)
#     if sym == :mp
#         Core.getfield(MP::MaterialProperties,sym::Symbol)
#     else
#         Core.getproperty(Core.getfield(MP,:mp),sym)
#     end
# end
#
# Base.setproperty!(MP::MaterialProperties,prop::Symbol,val) = Core.setproperty!(Core.getfield(MP,:mp),prop,val)
#
# function Base.propertynames(MP::MaterialProperties)
#     properties = Set{Symbol}()
#     for fields in fieldMP.mp
#         push!(properties, propertynames())
#     end
#
# end
