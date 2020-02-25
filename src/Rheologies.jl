module Rheologies

    using JuAFEM ; const J = JuAFEM
    using Tensors
    using SparseArrays
    using BlockArrays
    import Base: show

    include("types.jl")
    include("interface.jl")
    include("assembling.jl")
    include("model.jl")

    # exports #
    # types
    export Rheology, BoundaryConditions
    export Viscosity
    export Elasticity, IncompressibleElasticity, CompressibleElasticity
    export Plasticity

    # functions
    export solve

end # module
