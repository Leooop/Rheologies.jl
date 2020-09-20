using Rheologies;
const R = Rheologies;

include(joinpath(@__DIR__, "seeded_crust_gmsh.jl"))

output_path = @__DIR__
## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

Lx = 1.0
Ly = 0.1
H_bdt = 1e3
n_seeds = 0
radius_bounds = (300, 400)
#el_geom = Triangle # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

nx = 100
ny = 10

corner1 = Vec(0.0, 0.0)
corner2 = Vec(Lx, Ly)
el_geom = Quadrilateral # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

# generate the grid and sets:
grid = JuAFEM.generate_grid(el_geom, (nx, ny), corner1, corner2)
# Add facesets / nodesets / cellsets on which bc will be applied (except for the top, bottom, left and right boundary that are implemented by default)



## VARIABLES AND INTERPOLATIONS
variables = PrimitiveVariables{1}((:u,), (2,), el_geom)
#variables = PrimitiveVariables{2}((:u,:p), (2,1), el_geom) # not good !!!!

quad_order = 2
quad_type = :legendre # or :lobatto

## BOUNDARY CONDITIONS ##
# DIRICHLET :
# KEY : set on which to apply bc
# VALUES :
# 1) constrained variable
# 2) constraint function of (x,t)
# 3) constrained components
# NEUMANN :
# KEY : set on which to apply bc
# VALUES : traction as a function of (x,t)
gaussian_shape(x,A,xm,std) = A*exp(-0.5*(x[1]-xm)^2/std^2)
load_rate = -1e6
gaussian_trac(x) = gaussian_shape(x,load_rate,Lx,0.05)
bc = BoundaryConditions(
    dirichlet = Dict("left" => [:u, (x,t)-> Vec(0.0,0.0), [1,2]]),
    neumann = Dict("right" => (x,t)-> Vec(0.0,gaussian_trac(x)*t))
    )
#Dict("top" => (x,t)->Vec(0.0,-1e5),
#"bottom" => (x,t)->Vec(0.0,1e5)) )

bf = BodyForces([0.0, 0.0])#-9.81*2700])

## MATERIAL RHEOLOGY :

# define a rheology distribution function ...
# function multi_seeds(x,x0,radius,val_in,val_out)
#     for i in eachindex(radius)
#         (sqrt((x[1]-x0[i][1])^2 + (x[2]-x0[i][2])^2) <= radius[i]) && return val_in
#     end
#     return val_out
# end
# multi_seeds(x,val_in,val_out) = multi_seeds(x,x0,radius,val_in,val_out)
# and use it ...
seed(x,x0,r,val_in::T,val_out::T) where {T<:Real} = (sqrt((x[1]-x0[1])^2 + (x[2]-x0[2])^2) <= r ?
                                   val_in : val_out)

elas = Elasticity(E = 70e9, ν = 0.3)

visco = Viscosity(η = x -> (x[2] >= Ly - H_bdt) ? 1e24 : 1e17)

Δσ = 1e3
plas = DruckerPrager(
    ϕ = 30.0,
    C = 1e6,#x -> seed(x,[0.1Lx,0.9Ly],0.1*Ly,1e4,1e6),
    H = -1,
    ηᵛᵖ = 1e5,
)

VEP_rheology = Rheology(elasticity = elas, plasticity = plas)

### CLOCK ###
clock = Clock(tspan = (0.0, 40), Δt = 1, Δt_max = 20)

### SOLVER ###
linsolver = BackslashSolver()
nlsolver = NewtonRaphson(max_iter_number = 10, atol = 1e-3, linear_solver = linsolver)


### MODEL ###
VEP_model = Model(
    grid = grid,
    variables = variables,
    quad_order = quad_order,
    quad_type = quad_type,
    bc = bc,
    body_forces = bf,
    rheology = VEP_rheology,
    initial_state = nothing,
    clock = clock,
    solver = nlsolver,
)

### OUTPUT ###
output_name = "point_load_plate"
VEP_outputs = Dict(
    :σxx => (r, s) -> s.σ[1, 1],
    :σyy => (r, s) -> s.σ[2, 2],
    :σzz => (r, s) -> s.σ[3, 3],
    :acum_ep => (r,s)-> s.ϵ̅ᵖ,
    #:η => (r, s) -> r.viscosity.η,
    :C => (r,s)-> r.plasticity.C,
    #:E => (r,s)-> r.elasticity.E,
    :τ => (r,s) -> R.get_τ(dev(s.σ),r.plasticity),
    :F => (r,s) -> R.get_τ(dev(s.σ),r.plasticity) - (1/3 * tr(s.σ) + r.plasticity.ξ*(r.plasticity.C + r.plasticity.H*s.ϵ̅ᵖ)),
    :p => (r,s) -> -1 / 3 * tr(s.σ),
)#,
#:τoverp  => (r,s)-> R.get_τ(dev(s.σ),r.plasticity)/(-1/3 * tr(s.σ)) )

# C = r.plasticity.C
# H = r.plasticity.H
# η = r.plasticity.η
# ξ = r.plasticity.ξ

VEP_ow = VTKOutputWriter(
    VEP_model,
    joinpath(output_path, output_name),
    VEP_outputs,
    interval = 1,
    force_path = true,
)

### SOLVE ###
@time modelsol, u = solve(VEP_model, output_writer = VEP_ow, log = true);
