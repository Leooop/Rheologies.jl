
using Base: @kwdef

### BEHAVIORS ###
# Higher level supertype :
abstract type AbstractBehavior{T} end

abstract type Viscosity{T} <: AbstractBehavior{T} end
abstract type Elasticity{T} <: AbstractBehavior{T} end
abstract type Plasticity{T} <: AbstractBehavior{T} end

@kwdef struct IncompressibleElasticity{T} <: Elasticity{T}
    E::T
    ν::T
    function IncompressibleElasticity(E,ν)
        if ν <= 0.49
            ν = 0.4999999
            @warn "Poisson ratio can not be lower than 0.49 to ensure that the material is
            close to incompressible. Set to 0.4999999 by default"
        end
        new{typeof(E)}(E,convert(typeof(E),ν))
    end
end

@kwdef struct CompressibleElasticity{T} <: Elasticity{T}
    E::T
    ν::T
end

function Base.show(io::IO, ::MIME"text/plain", elast::E) where{E<:Elasticity}
    println(io, "$(typeof(elast))")
    println(io, "         E : $(elast.E)")
    println(io, "         ν : $(elast.ν)")
end

get_G(e::Elasticity) = e.E/(2*(1+e.ν))
get_K(e::Elasticity) = e.E/(3*(1-2*e.ν))
## RHEOLOGY ###
Maybe(T)=Union{T,Nothing}

struct Rheology{ T<:AbstractFloat,
        V<:Maybe(Viscosity{T}),
        E<:Maybe(Elasticity{T}),
        P<:Maybe(Plasticity{T}) }
    viscosity::V
    elasticity::E
    plasticity::P
end

Rheology(T ;viscosity, elasticity, plasticity) = Rheology(T,viscosity, elasticity, plasticity)

function Base.show(io::IO, ::MIME"text/plain", rheology::Rheology)
    println(io, "⌗  Rheology instance")
    print(io, "    -> viscosity  : ")
    show(io,rheology.viscosity); print(io,"\n")
    print(io, "    -> elasticity : ")
    show(io,MIME"text/plain"(),rheology.elasticity);
    print(io, "    -> plasticity : ")
    show(io,rheology.viscosity);
end

rheology_summary(r::Rheology{T,Nothing,E,Nothing}) where{T,E<:Elasticity} = "elastic"
rheology_summary(r::Rheology{T,V,E,Nothing}) where{T,V<:Viscosity,E<:Elasticity} = "visco-elastic"
rheology_summary(r::Rheology{T,V,Nothing,Nothing}) where{T,V<:Viscosity} = "viscous"
rheology_summary(r::Rheology{T,Nothing,E,P}) where{T,E<:Elasticity,P<:Plasticity} = "elasto-plastic"
rheology_summary(r::Rheology{T,V,E,P}) where{T,V<:Viscosity,E<:Elasticity,P<:Plasticity} = "visco-elasto-plastic"
### BOUNDARY CONDITIONS ###

@kwdef struct BoundaryConditions
    dirichlet::Dict
    neumann::Maybe(Dict)
end
BoundaryConditions(::Nothing,::Dict) = @error "A minimal set of Dirichlet boundary conditions must be provided in order to prevent rigid motion of the spatial domain"

### SOLVED VARIABLES ###

@kwdef struct Variables{Nvars}
    names::Union{Symbol,NTuple{Nvars,Symbol}}
    interpolations::Union{Interpolation,NTuple{Nvars,Interpolation}}
    function Variables{Nvars}(names,ips) where {Nvars}
        (typeof(names) == Symbol) && (@assert Nvars == 1)
        (typeof(names) == Tuple) && (@assert Nvars == length(names) == length(ips))
        return new(names,ips)
    end
end
# function Variables(names::Symbol,ips::Interpolation)
#     Variables{1}(names,ips)
# end
# function Variables(names::NTuple{N,Symbol},ips::NTuple{N,Interpolation}) where {N}
#     Variables{N}(names,ips)
# end
refshape(::Type{Cell{1,N,2}}) where{N} = RefCube
refshape(::Type{Cell{2,N,4}}) where{N} = RefCube
refshape(::Type{Cell{3,N,6}}) where{N} = RefCube
refshape(::Type{Cell{2,N,3}}) where{N} = RefTetrahedron
refshape(::Type{Cell{3,N,4}}) where{N} = RefTetrahedron

function Variables{NV}(names,ip_order,el_geom::Type{Cell{dim,N,M}}) where {NV,dim,N,M}
    NV == 1 && return Variables{NV}(names, Lagrange{dim,refshape(el_geom),ip_order}())
    NV == 2 && return Variables{NV}(names, (Lagrange{dim,refshape(el_geom),ip_order[1]}(),
                                          Lagrange{dim,refshape(el_geom),ip_order[2]}()) )
    @error "More than 2 variables (:u,:p) is not implemented yet"
end



### Time Handler ###

@kwdef mutable struct Clock{T}
    t0::T = 0.0
    current_time::T = 0.0
    cfl::T
    Δt::T = 0.0
    #Clock(t0,current_time,cfl,Δt) = new{typeof(t0)}(t0,current_time,cfl,Δt)
end

Clock{Nothing}() = Clock{Nothing}(nothing,nothing,nothing,nothing)
#Clock(cfl) = Clock(zero(typeof(clf)),zero(typeof(clf)),cfl,zero(typeof(clf)))

function Base.show(io::IO, ::MIME"text/plain", c::Clock{T}) where {T}
    if T <: AbstractFloat
        println(io, "⌗  Clock instance")
        println(io, "    -> t0 : $(c.t0)")
        println(io, "    -> current_time : $(c.current_time)")
        println(io, "    -> clf : $(c.cfl)")
        println(io, "    -> Δt : $(c.Δt)")
    elseif T == Nothing
        println(io, "⌗  Empty Clock instance : No time dependency of the problem")
    end
end
