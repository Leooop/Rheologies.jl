### MODEL ###
# initialize all JuAFEM objects
#const CVT = Union{Tuple{CellVectorValues,CellScalarValues},Tuple{CellVectorValues}}
"""
    Model{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,BF,SP,C}

A `Model` instance contains all the informations needed to run a simulation.
Type parameters associated to it are used to dispatch methods with respect to model specificities, especially {dim,N,D,V,E,P,S} (see description below).

# Fields
- `grid::Grid{dim}` : JuAFEM.Grid instance
- `dofhandler::DH` : JuAFEM.DofHandler instance containing information about primitive variables discretization over the grid
- `dirichlet_bc::CH` : JuAFEM.ConstraintHandler instance containing dirichlet boundary conditions
- `neumann_bc::NBC` : dictionnary with `set_names`, `traction(x::Vector)` pairs
- `body_forces::Vector{BC}` : body forces for every cell
- `cellvalues_tuple::CVT` : A tuple with as many elements as primitive variables in the model, with each one containing in cell integrations related informations
- `facevalues::FV` : Same as cellvalues but for tractions integration on faces
- `material_properties::Vector{Rheology{Float64,D,V,E,P}}` : A Vector of `Rheology`s instance associated to each cell
- `material_state::Vector{Vector{MS}}` : `MaterialState`s for each integration point
- `K::SP` : global stiffness matrix
- `RHS::Vector{Float64}` : global right hand side vector (external force vector or residual vector in incremental formulation)
- `clock::C` : `Clock` instance containing all time related informations
- `solver::S` : `Solver` subtype instance containing informations about the solving procedure
- `multithreading::Bool` : `True` to use multithreaded assembly

# Type parameters
- `dim::Int` : spatial dimensions of the model
- `N::Int` : number of primitives variables. Commonly 1 (displacement) or 2 (displacement and pressure)
- `D<:Damage` : micromechanical properties of the rheology
- `V<:Viscosity` : viscous properties of the rheology
- `E<:Elasticity` : elastic properties of the rheology
- `P<:Plasticity` : plastic properties of the rheology
- `S<:AbstractSolver` : linear or non linear solver
"""
struct Model{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,BF,SP,C}
    grid::Grid{dim}
    dofhandler::DH
    dirichlet_bc::CH
    neumann_bc::NBC # TODO : use a struct similar to ConstraintHandler may be relevant
    body_forces::Vector{BF}
    cellvalues_tuple::CVT
    facevalues::FV
    material_properties::Vector{Rheology{Float64,D,V,E,P}}
    material_state::Vector{Vector{MS}}
    K::SP
    RHS::Vector{Float64}
    clock::C
    solver::S
    multithreading::Bool
    function Model(grid::Grid{dim},
                   dh::DH,
                   dbc::CH,
                   nbc::NBC,
                   bf::Vector{BF},
                   cv_tuple::CVT,
                   fv::FV,
                   mp::Vector{Rheology{Float64,D,V,E,P}},
                   ms::Vector{Vector{MS}},
                   K::SP,
                   RHS::Vector{Float64},
                   clock::C,
                   solver::S,
                   multithreading::Bool) where {dim,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,BF,SP,C}

        N = length(cv_tuple)
        return new{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,BF,SP,C}(grid,dh,dbc,nbc,bf,cv_tuple,fv,mp,ms,K,RHS,clock,solver,multithreading)
    end
end

# Model{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,SP,C}(grid,dh,dbc,nbc,cv_tuple,fv,mp,ms,K,RHS,clock,solver,multithreading)  = Model(grid::Grid{dim},
#                dh::DH,
#                dbc::CH,
#                nbc::NBC,
#                cv_tuple::CVT,
#                fv::FV,
#                mp::Vector{Rheology{Float64,D,V,E,P}},
#                ms::Vector{Vector{MS}},
#                K::SP,
#                RHS::Vector{Float64},
#                clock::C,
#                solver::S,
#                multithreading::Bool) where {dim,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,SP,C}
#

# Constructor using input file user defined variables
"""
    Model(grid::Grid, variables::PrimitiveVariables, quad_order::Int, quad_type::Symbol, bc::BoundaryConditions, body_forces::BodyForces, rheology::Rheology, clock::Clock, solver::Solver, multithreading::Bool = false)

`Model` constructor. Keyword arguments can also be used but all arguments are mandatory.

# Fields
- `grid::Grid{dim}` : JuAFEM.Grid instance
- `variables::PrimitiveVariables` : objet containing informations about primitive variables and their interpolation
- `quad_order::Int` : quadrature order to evaluate the element integrals
- `quad_type::Symbol` : quadrature type. Can be `:legendre` or `:lobatto`
- `bc::BoundaryConditions` : contains dictionaries mapping boundary condition functions to sets of nodes/elements/faces.
- `body_forces::BodyForces` : contains body forces as an homogeneous vectorial quantity or as a function of a space location vector.
- `rheology::Rheology` : contains material properties
- `clock::Clock` : `Clock` instance containing all time related informations
- `solver::AbstractSolver` : `Solver` subtype instance containing informations about the solving procedure
- `multithreading::Bool` : `True` to use multithreaded assembly

"""
function Model(grid::Grid, variables::PrimitiveVariables,
               quad_order::Int, quad_type::Symbol,
               bc::BoundaryConditions, body_forces::BodyForces, rheology::Rheology,
               initial_state, clock::Clock, solver::AbstractSolver, multithreading)

    dh, dbc, cv_tuple, fv, mp, ms, bf, K, RHS = setup_model(grid, variables, quad_order, quad_type, bc, body_forces, initial_state, rheology)
    nbc = bc.neumann
    return Model(grid,dh,dbc,nbc,bf,cv_tuple,fv,mp,ms,K,RHS,clock,solver,multithreading)
end

Model(; grid, variables, quad_order, quad_type, bc, body_forces, rheology, initial_state, clock, solver, multithreading = false) =
         Model(grid, variables, quad_order, quad_type, bc, body_forces, rheology, initial_state, clock, solver, multithreading)

# Ask for a non linear solver if rheology is damaged and/or plastic
# NLRheology = Union{Rheology{T,Nothing,Nothing,E,P} where {T,E,P<:Plasticity},
#                    Rheology{T,D,Nothing,E,P} where {T,D<:Damage,E,P<:Plasticity},
#                    Rheology{T,D,Nothing,E,Nothing} where {T,E,D<:Plasticity}}
# Model(grid, variables, quad_order, quad_type, bc, rheology::NLRheology, clock, solver::S, multithreading) where {S<:LinearSolver} = @error "Non linear rheology, please provide a non linear solver"


Base.show(io::IO, grid::Grid) = Base.show(io::IO, MIME"text/plain"(), grid::Grid)

function Base.show(io::IO, ::MIME"text/plain", model::Model{dim,D,V,E,P,N,S}) where {dim,D,V,E,P,N,S}
    print(io, "⌗  Model to be solved for $(string(model.dofhandler.field_names)) using $(model.dofhandler.field_interpolations) on a grid $(model.grid) \n")
    #print(io,"\n")

end
