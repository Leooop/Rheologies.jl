module Rheologies

    using Reexport

    #@reexport using Tensors
    @reexport using JuAFEM
    const J = JuAFEM
    import JuAFEM: QuadratureRule

    using InteractiveUtils: subtypes
    using SparseArrays
    using BlockArrays
    using StructArrays
    using StaticArrays

    include("types.jl")
    include("interface.jl")
    include("model.jl")
    include("iterations.jl")
    include("assemble/assembling.jl")

    # exports #
    # types
    export MaterialProperties, Rheology, BoundaryConditions, Variables, Clock
    export Damage
    export Viscosity
    export Elasticity, IncompressibleElasticity, CompressibleElasticity
    export Plasticity, IdealPlasticity, AssociatedDruckerPrager,
           NonAssociatedDruckerPrager, VonMises, ViscoLinearHardeningDP
    export Damage, BRSDamage

    # functions
    export run_simulation, setup_model, create_material_properties

end # module
