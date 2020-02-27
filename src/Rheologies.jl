module Rheologies

    using JuAFEM ; const J = JuAFEM
    import JuAFEM: QuadratureRule
    using Tensors
    using SparseArrays
    using BlockArrays
    import Base: show

    include("types.jl")
    include("interface.jl")
    include("model.jl")
    include("iterations.jl")
    include("assembling.jl")

    # exports #
    # types
    export Rheology, BoundaryConditions, Variables, Clock
    export Viscosity
    export Elasticity, IncompressibleElasticity, CompressibleElasticity
    export Plasticity

    # functions
    export run_simulation, setup_model

end # module
