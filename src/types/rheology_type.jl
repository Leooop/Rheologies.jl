## RHEOLOGY ###

@kwdef struct Rheology{T, D, V, E, P}
    damage::D = nothing
    viscosity::V = nothing
    elasticity::E = nothing
    plasticity::P = nothing
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
Rheology{D,V,E,P}(damage, viscosity, elasticity, plasticity) where {D,V,E,P} =
        Rheology{Float64,D,V,E,P}(damage, viscosity, elasticity, plasticity)
Rheology{T}(damage, viscosity, elasticity, plasticity) where {T} =
        Rheology(damage, viscosity, elasticity, plasticity)
#Rheology(; damage, viscosity, elasticity, plasticity) = Rheology(damage, viscosity, elasticity, plasticity)

" Evaluate the functionaly defined rheology at coordinates x"
function (r::Rheology{Function})(x)
    behaviors = []
    for prop in propertynames(r)
        behavior = getproperty(r,prop)
        properties = []
        if typeof(behavior) <: AbstractBehavior{Function}
            for b_prop in propertynames(behavior)
                bp = getproperty(behavior,b_prop)
                if isa(bp, Function)
                    try
                    push!(properties,bp(x))
                    catch err
                        @error "An error occured during the evaluation of a field function. Make sure that the length of the coordinate argument of your defined functions matchs the dimensions of your grid"
                        throw(err)
                    end
                elseif b_prop == :Dᵉ
                    push!(properties,get_elastic_stiffness_tensor(properties[3],properties[5]))
                else
                    push!(properties,bp)
                end
            end
            push!(behaviors,updatetypeparameter(behavior,properties))
        else
            push!(behaviors,behavior)
        end
    end
    return Rheology(behaviors...)
end

(r::Rheology{T})(x) where {T<:Real} = @info "Rheology instance is already filled with Real values"


function Base.show(io::IO, ::MIME"text/plain", rheology::Rheology{T,D,V,E,P}) where {T,D,V,E,P}
    if T <: Real
        println(io, "⌗  $(typeof(rheology)) instance ")
    elseif T == Function
        println(io, "⌗  $(typeof(rheology)) instance with functional parameters ")
    end
    print(io, "    -> damage  : ")
    show(io,MIME"text/plain"(),rheology.damage); print(io,"\n")
    print(io, "    -> viscosity  : ")
    show(io,MIME"text/plain"(),rheology.viscosity); print(io,"\n")
    print(io, "    -> elasticity : ")
    show(io,MIME"text/plain"(),rheology.elasticity);
    print(io, "    -> plasticity : ")
    show(io,MIME"text/plain"(),rheology.plasticity);
end

rheology_summary(r::Rheology{T,Nothing,Nothing,E,Nothing}) where{T,E<:Elasticity} = "elastic"
rheology_summary(r::Rheology{T,Nothing,V,E,Nothing}) where{T,V<:Viscosity,E<:Elasticity} = "visco-elastic"
rheology_summary(r::Rheology{T,Nothing,V,Nothing,Nothing}) where{T,V<:Viscosity} = "viscous"
rheology_summary(r::Rheology{T,Nothing,Nothing,E,P}) where{T,E<:Elasticity,P<:Plasticity} = "elasto-plastic"
rheology_summary(r::Rheology{T,Nothing,V,E,P}) where{T,V<:Viscosity,E<:Elasticity,P<:Plasticity} = "visco-elasto-plastic"
rheology_summary(r::Rheology{T,D,Nothing,E,Nothing}) where{T,D<:Damage,E<:Elasticity} = "damaged-elastic"
rheology_summary(r::Rheology{T,D,Nothing,E,P}) where{T,D<:Damage,E<:Elasticity,P<:Plasticity} = "damaged-elasto-plastic"
