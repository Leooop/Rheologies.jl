
using Base: @kwdef
import Base: show, getproperty, setproperty!
### BEHAVIORS ###
# Higher level supertype :
abstract type AbstractBehavior{T} end

abstract type Damage{T} <: AbstractBehavior{T} end
abstract type Viscosity{T} <: AbstractBehavior{T} end
abstract type Elasticity{T} <: AbstractBehavior{T} end
abstract type Plasticity{T} <: AbstractBehavior{T} end

abstract type IdealPlasticity{T} <: Plasticity{T} end
abstract type NonIdealPlasticity{T} <: Plasticity{T} end

abstract type ViscoPlasticity{T} <: NonIdealPlasticity{T} end

Maybe(T)=Union{T,Nothing}

## Elastic concrete types :
@kwdef mutable struct IncompressibleElasticity{T} <: Elasticity{T}
    E::Union{Float64,Function}
    ν::Float64 = 0.4999999
    function IncompressibleElasticity(E,ν::Float64)
        if  (ν >= 0.50) || (ν < 0.49)
            ν = 0.4999999
            @warn "Poisson ratio can not be lower than 0.49 and >=0.5 to ensure that the material is close to incompressible. Set to 0.4999999 by default"
        end
        T2 = isa(E,Function) ? Function : typeof(E)
        new{T2}(E,ν)
    end
end
IncompressibleElasticity{T}(E,ν) where {T} = IncompressibleElasticity(E,ν)

@kwdef mutable struct CompressibleElasticity{T} <: Elasticity{T}
    E::Union{Float64,Function}
    ν::Union{Float64,Function}
    function CompressibleElasticity(E,ν)
        T2 = any(isa.((E,ν),Function)) ? Function : typeof(E)
        new{T2}(E,ν)
    end
end
CompressibleElasticity{T}(E,ν) where {T} = CompressibleElasticity(E,ν)

function Base.show(io::IO, ::MIME"text/plain", elast::E) where{E<:Elasticity}
    print(io, "$(typeof(elast))\n",
    "├── E : $(elast.E)\n",
    "└── ν : $(elast.ν)\n")
end

get_G(e::Elasticity{T}) where {T<:AbstractFloat} = e.E/(2*(1+e.ν))
get_K(e::Elasticity{T}) where {T<:AbstractFloat} = e.E/(3*(1-2*e.ν))



## Plasticity concrete types :
@kwdef mutable struct AssociatedDruckerPrager{T} <: IdealPlasticity{T}
    μ::Maybe(Union{Float64,Function}) = nothing
    ϕ::Maybe(Union{Float64,Function}) = nothing
    ψ::Maybe(Union{Float64,Function}) = nothing
    C::Maybe(Union{Float64,Function}) = 0.0
    function AssociatedDruckerPrager(μ,ϕ,ψ,C)
        T2 = any(isa.((μ,ϕ,ψ,C),Function)) ? Function : Float64
        if ϕ == nothing
            if μ != nothing
                if μ isa Function
                    return new{T2}(μ, x->atand(μ(x)), x->atand(μ(x)), C)
                elseif μ isa Real
                    μ = Float64(μ)
                    return new{T2}(μ, atand(μ), atand(μ), C)
                else
                    @error "Unexpected type of field μ"
                end
            else
                @error "At least ϕ or μ must be provided"
            end
        elseif ϕ isa Real
            if μ isa Function
                @error "Ambiguity between μ and ϕ, please provide only one of them"
            elseif μ isa Real
                (ψ == nothing) && (@error "Ambiguity between μ and ϕ, please provide only one of them")
                (ψ isa Real) && (return new{T2}(Float64(μ),Float64(ϕ),Float64(ψ),Float64(C)))
            elseif μ == nothing
                ϕ = Float64(ϕ)
                return new{T2}(tand(ϕ), ϕ, ϕ, C)
            else
                @error "Unexpected type of field μ"
            end
        elseif ϕ isa Function
            if μ isa Real
                @error "Ambiguity between μ and ϕ, please provide only one of them"
            elseif μ == nothing
                return new{T2}(x->tand(ϕ(x)), ϕ, ϕ, C)
            elseif μ isa Function
                @error "Ambiguity between μ and ϕ, please provide only one of them"
            else
                @error "Unexpected type of field μ"
            end
        end
    end
end
AssociatedDruckerPrager{T}(μ,ϕ,ψ,C) where {T} = AssociatedDruckerPrager(μ,ϕ,ψ,C)



@kwdef mutable struct NonAssociatedDruckerPrager{T} <: IdealPlasticity{T}
    μ::Union{Float64,Function} = 0.0
    ϕ::Union{Float64,Function} = 0.0
    ψ::Union{Float64,Function} = 0.0
    C::Union{Float64,Function} = 0.0
    function NonAssociatedDruckerPrager(μ,ϕ,ψ,C)
        T2 = any(isa.((μ,ϕ,ψ,C),Function)) ? Function : typeof(μ)
        if ϕ == 0
            if μ != 0
                if μ isa Function
                    return new{T2}(μ, x->atand(μ(x)), ψ, C)
                elseif μ isa Float64
                    return new{T2}(μ, atand(μ), ψ, C)
                else
                    @error "Unexpected type of field μ"
                end
            else
                @error "At least ϕ or μ must be provided"
            end
        else
            #μ != 0 && @error "Ambiguity between μ and ϕ, please provide only one of them"
            if ϕ isa Function
                return new{T2}(x->tand(ϕ(x)), ϕ, ψ, C)
            elseif ϕ isa Float64
                return new{T2}(tand(ϕ), ϕ, ψ, C)
            else
                @error "Unexpected type of field ϕ"
            end
        end
    end
end
NonAssociatedDruckerPrager{T}(μ,ϕ,ψ,C) where {T} = NonAssociatedDruckerPrager(μ,ϕ,ψ,C)

function Base.show(io::IO, ::MIME"text/plain",
                   plast::Union{AssociatedDruckerPrager,NonAssociatedDruckerPrager})
    print(io, "$(typeof(plast))\n",
    "├── μ : $(plast.μ)\n",
    "├── ϕ : $(plast.ϕ)\n",
    "├── ψ : $(plast.ψ)\n",
    "└── C : $(plast.C)\n")
end



@kwdef mutable struct ViscoLinearHardeningDP{T} <: ViscoPlasticity{T}
    μ::Union{Float64,Function} = 0.0
    ϕ::Union{Float64,Function} = 0.0
    ψ::Union{Float64,Function} = 0.0
    C::Union{Float64,Function} = 0.0
    H::Union{Float64,Function} = 0.0
    ηᵛᵖ::Union{Float64,Function} = 0.0
    function ViscoLinearHardeningDP(μ,ϕ,ψ,C,H,ηᵛᵖ)
        T2 = any(isa.((μ,ϕ,ψ,C,H,ηᵛᵖ),Function)) ? Function : typeof(μ)
        if ϕ == 0
            if μ != 0
                if μ isa Function
                    return new{T2}(μ, x->atand(μ(x)), ψ, C, H, ηᵛᵖ)
                elseif μ isa Float64
                    return new{T2}(μ, atand(μ), ψ, C, H, ηᵛᵖ)
                else
                    @error "Unexpected type of field μ"
                end
            else
                @error "At least ϕ or μ must be provided"
            end
        else
            #μ != 0 && @error "Ambiguity between μ and ϕ, please provide only one of them"
            if ϕ isa Function
                return new{T2}(x->tand(ϕ(x)), ϕ, ψ, C, H, ηᵛᵖ)
            elseif ϕ isa Float64
                return new{T2}(tand(ϕ), ϕ, ψ, C, H, ηᵛᵖ)
            else
                @error "Unexpected type of field ϕ"
            end
        end
    end
end
ViscoLinearHardeningDP{T}(μ,ϕ,ψ,C,H,ηᵛᵖ) where {T} = ViscoLinearHardeningDP(μ,ϕ,ψ,C,H,ηᵛᵖ)


@kwdef mutable struct VonMises{T} <: IdealPlasticity{T}
    C::Union{Float64,Function}
    function VonMises(C)
        T2 = isa(C,Function) ? Function : typeof(C)
        new{T2}(C)
    end
end
VonMises{T}(C) where {T} = VonMises(C)

function Base.show(io::IO, ::MIME"text/plain",
                   plast::VonMises)
    print(io, "$(typeof(plast))\n",
    "└── C : $(plast.C)")
end

## Damage concrete types :

@kwdef mutable struct BRSDamage{T} <: Damage{T}
    μ::Union{Float64,Function} = 0.6 # Friction coef
    β::Union{Float64,Function} = 0.1 # Correction factor
    K₁c::Union{Float64,Function} = 1.74e6 # Critical stress intensity factor (Pa.m^(1/2))
    a::Union{Float64,Function} # Initial flaw size (m)
    ψ::Union{Float64,Function} # crack angle to the principal stress (radians)
    D₀::Union{Float64,Function} # Initial flaw density
    n::Union{Float64,Function} = 34.0 # Stress corrosion index
    l̇₀::Union{Float64,Function} = 0.24 # Ref. crack growth rate (m/s)
    H::Union{Float64,Function} = 50e3 # Activation enthalpy (J/mol)
    A::Union{Float64,Function} = 5.71 # Preexponential factor (m/s)
    function BRSDamage(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)
        T2 = any(isa.((μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A),Function)) ? Function : typeof(β)
        new{T2}(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)
    end
end
BRSDamage{T}(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A) where {T} = BRSDamage(μ,β,K₁c,a,ψ,D₀,n,l̇₀,H,A)


## RHEOLOGY ###

struct Rheology{T, D<:Maybe(Damage),
                   V<:Maybe(Viscosity),
                   E<:Maybe(Elasticity),
                   P<:Maybe(Plasticity) }
    damage::D
    viscosity::V
    elasticity::E
    plasticity::P
    function Rheology(damage, viscosity, elasticity, plasticity)
        fields = (damage, viscosity, elasticity, plasticity)
        cond1 = .!(<:).(typeof.(fields),AbstractBehavior{Function})
        cond2 = (<:).(typeof.(fields),AbstractBehavior)
        idxs = findall(cond1)
        if length(idxs) == 4
            if any(cond1 .& cond2)
                no_func_field = fields[findfirst(cond1 .& cond2)]
                T = typeof(getproperty(no_func_field,propertynames(no_func_field)[1]))
            else
                @error "Rheology instance can not be empty"
            end
        else
            T = Function
        end
        new{T,typeof(damage),typeof(viscosity),typeof(elasticity),typeof(plasticity)}(damage, viscosity, elasticity, plasticity)
    end
end
Rheology{T,D,V,E,P}(damage, viscosity, elasticity, plasticity) where {T,D,V,E,P} =
        Rheology(damage, viscosity, elasticity, plasticity)
Rheology{T}(damage, viscosity, elasticity, plasticity) where {T<:AbstractFloat} =
        Rheology(damage, viscosity, elasticity, plasticity)
Rheology(; damage, viscosity, elasticity, plasticity) = Rheology(damage, viscosity, elasticity, plasticity)

" Evaluate the functionaly defined rheology at coordinates x"
function (r::Rheology{Function})(x)
    args = []
    for prop in propertynames(r)
        behavior = deepcopy(getproperty(r,prop))
        if typeof(getproperty(r,prop)) <: AbstractBehavior{Function}
            for b_prop in propertynames(behavior)
                bp = getproperty(behavior,b_prop)
                if isa(bp, Function)
                    try
                    setproperty!(behavior,b_prop,bp(x))
                    catch err
                        @error "An error occured during the evaluation of a field function. Make sure that the length of the coordinate argument of your defined functions matchs the dimensions of your grid"
                        throw(err)
                    end
                end
            end
            push!(args,updatetypeparameter(behavior))
        else
            push!(args,behavior)
        end
    end
    return Rheology(args...)
end

(r::Rheology{T})(x) where {T<:Real} = @info "Rheology instance is already filled with Real values"
#(r::Rheology{T})(x) where {T<:Real} = @info "Rheology instance is already filled with Real values"


function Base.show(io::IO, ::MIME"text/plain", rheology::Rheology{T,D,V,E,P}) where {T,D,V,E,P}
    if T <: Real
        println(io, "⌗  $(typeof(rheology)) instance ")
    elseif T == Function
        println(io, "⌗  $(typeof(rheology)) instance with functional parameters ")
    end
    print(io, "    -> damage  : ")
    show(io,rheology.damage); print(io,"\n")
    print(io, "    -> viscosity  : ")
    show(io,rheology.viscosity); print(io,"\n")
    print(io, "    -> elasticity : ")
    show(io,MIME"text/plain"(),rheology.elasticity);
    print(io, "    -> plasticity : ")
    show(io,MIME"text/plain"(),rheology.plasticity);
end

rheology_summary(r::Rheology{Nothing,Nothing,E,Nothing}) where{E<:Elasticity} = "elastic"
rheology_summary(r::Rheology{Nothing,V,E,Nothing}) where{V<:Viscosity,E<:Elasticity} = "visco-elastic"
rheology_summary(r::Rheology{Nothing,V,Nothing,Nothing}) where{V<:Viscosity} = "viscous"
rheology_summary(r::Rheology{Nothing,Nothing,E,P}) where{E<:Elasticity,P<:Plasticity} = "elasto-plastic"
rheology_summary(r::Rheology{Nothing,V,E,P}) where{V<:Viscosity,E<:Elasticity,P<:Plasticity} = "visco-elasto-plastic"
rheology_summary(r::Rheology{D,Nothing,E,Nothing}) where{D<:Damage,E<:Elasticity} = "damaged-elastic"
rheology_summary(r::Rheology{D,Nothing,E,P}) where{D<:Damage,E<:Elasticity,P<:Plasticity} = "damaged-elasto-plastic"

## TODO : Check whether a material parameters type is needed
### MATERIAL PARAMETERS ###

struct MaterialProperties{T}
    mp :: StructArray
    function MaterialProperties(mp)
        T = eltype(mp)
        new{T}(mp)
    end
end
MaterialProperties{T}(mp) where T = MaterialProperties(mp)
function Base.getproperty(MP::MaterialProperties,sym::Symbol)
    if sym == :mp
        Core.getfield(MP::MaterialProperties,sym::Symbol)
    else
        Core.getproperty(Core.getfield(MP,:mp),sym)
    end
end

Base.setproperty!(MP::MaterialProperties,prop::Symbol,val) = Core.setproperty!(Core.getfield(MP,:mp),prop,val)

function Base.propertynames(MP::MaterialProperties)
    properties = Set{Symbol}()
    for fields in fieldMP.mp
        push!(properties, propertynames())
    end

end




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
    NV == 1 && return Variables{NV}(names, (Lagrange{dim,refshape(el_geom),ip_order[1]}(),))
    NV == 2 && return Variables{NV}(names, (Lagrange{dim,refshape(el_geom),ip_order[1]}(),
                                          Lagrange{dim,refshape(el_geom),ip_order[2]}()) )
    @error "More than 2 variables (:u and :p) is not implemented yet"
end



### Time Handler ###

@kwdef mutable struct Clock{T}
    tspan::Tuple{T,T}
    current_time::T = 0.0
    iter::Int = 0
    Δt_max::T = 1.0
    Δt::T = 1.0
    time_vec::Vector{T} = [0.0]
    Clock(tspan::Tuple{T,T},current_time::T,iter::Int,Δt_max::T,Δt::T,time_vec::Vector{T}) where {T} = new{T}(tspan,current_time,iter,Δt_max,Δt,time_vec)
end
Clock{T}(tspan,current_time,iter,Δt_max,Δt,time_vec) where {T} = Clock(tspan,current_time,iter,Δt_max,Δt,time_vec)
Clock(tspan,current_time,iter,Δt_max,Δt) = Clock(tspan,current_time,iter,Δt_max,Δt,[tspan[1]])
Clock(tspan::Tuple,Δt_max) = Clock(tspan,tspan[1],0,Δt_max,Δt_max)
Clock(tspan::Tuple) = Clock(tspan,tspan[1],0,1.0,1.0)
Clock{Nothing}() = Clock{Nothing}(nothing,nothing,nothing,nothing)

#Clock(t_max) = Clock(zero(typeof(clf)),zero(typeof(clf)),Δt_max,zero(typeof(clf)))

function Base.show(io::IO, ::MIME"text/plain", c::Clock{T}) where {T}
    if T <: AbstractFloat
        print(io, "⌗  Clock instance\n",
        "├── tspan : $(c.tspan)\n",
        "├── current_time : $(c.current_time)\n",
        "├── iter : $(c.iter)\n",
        "├── Δt_max : $(c.Δt_max)\n",
        "├── Δt : $(c.Δt)\n",
        "└── time_vec : $(c.time_vec)\n")
    elseif T == Nothing
        println(io, "⌗  Empty Clock instance : No time dependency of the problem")
    end
end



## UPDATE PARAMETERS :

function updatetypeparameter(b::CompressibleElasticity)
    properties = getproperty.(Ref(b),propertynames(b))
    return CompressibleElasticity(properties...)
end
function updatetypeparameter(b::IncompressibleElasticity)
    properties = getproperty.(Ref(b),propertynames(b))
    return IncompressibleElasticity(properties...)
end
function updatetypeparameter(b::AssociatedDruckerPrager)
    properties = getproperty.(Ref(b),propertynames(b))
    return AssociatedDruckerPrager(properties...)
end
function updatetypeparameter(b::NonAssociatedDruckerPrager)
    properties = getproperty.(Ref(b),propertynames(b))
    return NonAssociatedDruckerPrager(properties...)
end
function updatetypeparameter(b::ViscoLinearHardeningDP)
    properties = getproperty.(Ref(b),propertynames(b))
    return ViscoLinearHardeningDP(properties...)
end
function updatetypeparameter(b::VonMises)
    properties = getproperty.(Ref(b),propertynames(b))
    return VonMises(properties...)
end
function updatetypeparameter(b::BRSDamage)
    properties = getproperty.(Ref(b),propertynames(b))
    return BRSDamage(properties...)
end

# would be gould to generate these functions for each Behavior <: AbstractBehavior using code generation, but the following raises an UndefVarError "CompressibleElasticity not defined"

# for behavior in subtypes(AbstractBehavior)
#     for concrete_behavior in subtypes(behavior)
#         typ = Symbol(concrete_behavior)
#         println(typ)
#         eval(
#             :(function updatetypeparameter(b::$(typ))
#                 properties = getproperty.(Ref(b),propertynames(b))
#                 return $(typ)(properties...)
#             end)
#         )
#     end
# end
