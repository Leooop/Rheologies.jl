module Rheologies

    using Reexport

    #@reexport using Tensors
    @reexport using JuAFEM
    const J = JuAFEM
    import JuAFEM: QuadratureRule

    #using Parameters
    using InteractiveUtils: subtypes
    using TimerOutputs
    using Printf
    using SparseArrays
    using BlockArrays
    using StaticArrays
    using StructArrays
    using LinearAlgebra
    #import MUMPSjInv
    #import Preconditioners

    include("convenience_functions.jl")
    include("types/types.jl")
    include("physical_functions.jl")
    include("tangent_operators.jl")
    include("interface.jl")
    include("main.jl")
    include("linear_solve.jl")
    include("nonlinear_iterations.jl")
    include("assemble/assemble_global.jl")
    include("GmshParser.jl")

    # exports #
    # types
    export Model
    export MaterialState, BasicMaterialState, PlasticMaterialState, DamagedPlasticMaterialState
    export AbstractSolver, AbstractLinearSolver, AbstractNonLinearSolver, BackslashSolver, MUMPS, NewtonRaphson, ConjugateGradients,  MINRES, ILU, GMRESIterativeSolvers
    export Rheology, BoundaryConditions, BodyForces, Variables, Clock, PrimitiveVariables
    export Damage
    export Viscosity
    export Elasticity
    export Plasticity, DruckerPrager, VonMises
    export Damage, BRSDamage
    export VTKOutputWriter

    # functions
    export solve, setup_model, create_material_properties

end # module
