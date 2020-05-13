### MODEL ###
# initialize all JuAFEM objects
#const CVT = Union{Tuple{CellVectorValues,CellScalarValues},Tuple{CellVectorValues}}

struct Model{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,SP,C}
    grid::Grid{dim}
    dofhandler::DH
    dirichlet_bc::CH
    neumann_bc::NBC # TODO : use a struct similar to ConstraintHandler may be relevant
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
                   cv_tuple::CVT,
                   fv::FV,
                   mp::Vector{Rheology{Float64,D,V,E,P}},
                   ms::Vector{Vector{MS}},
                   K::SP,
                   RHS::Vector{Float64},
                   clock::C,
                   solver::S,
                   multithreading::Bool) where {dim,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,SP,C}

        N = length(cv_tuple)
        return new{dim,N,D,V,E,P,S,MS,CVT,FV,DH,CH,NBC,SP,C}(grid,dh,dbc,nbc,cv_tuple,fv,mp,ms,K,RHS,clock,solver,multithreading)
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
function Model(grid::Grid, variables::PrimitiveVariables,
               quad_order::Int, quad_type::Symbol,
               bc::BoundaryConditions, rheology::Rheology,
               clock::Clock, solver::Solver, multithreading)

    dh, dbc, cv_tuple, fv, mp, ms, K, RHS = setup_model(grid::Grid, variables::PrimitiveVariables,
                                           quad_order::Int, quad_type::Symbol,
                                           bc::BoundaryConditions, rheology::Rheology)
    nbc = bc.neumann
    return Model(grid,dh,dbc,nbc,cv_tuple,fv,mp,ms,K,RHS,clock,solver,multithreading)
end

Model(; grid, variables, quad_order, quad_type, bc, rheology, clock, solver, multithreading = false) =
         Model(grid, variables, quad_order, quad_type, bc, rheology, clock, solver, multithreading)

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
