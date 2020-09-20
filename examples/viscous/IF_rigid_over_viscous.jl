using Rheologies;
const R = Rheologies;

output_path = @__DIR__
## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

Lx = 20e3
Ly = 10e3
H_bdt = 1e3
n_seeds = 0
radius_bounds = (300, 400)
#el_geom = Triangle # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

nx = 50
ny = 50

corner1 = Vec(0.0, 0.0)
corner2 = Vec(Lx, Ly)
el_geom = Quadrilateral # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

# generate the grid and sets:
grid = JuAFEM.generate_grid(el_geom, (nx, ny), corner1, corner2)
# Add facesets / nodesets / cellsets on which bc will be applied (except for the top, bottom, left and right boundary that are implemented by default)
addnodeset!(
    grid,
    "clamped_lateral",
    x -> ((x[1] in (corner1[1], corner2[1])) & (x[2] == Ly / 2)),
);
addnodeset!(grid, "clamped_up", x -> ((x[1] == Lx/2) & (x[2] == Ly)));


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

cos_shape(x,A,Lx) = -A*cos(2π/Lx * x[1])
year = 3600 * 24 * 365.0
v_in_max = 0.01 / year # 1cm/year
cos_shape(x) = cos_shape(x,v_in_max,Lx)

bc = BoundaryConditions(
    dirichlet = Dict(
        "top" => [:u, (x,t)-> 0.0, 2],
        "bottom" => [:u, (x, t) -> cos_shape(x) * t, 2],
        "clamped_up" => [:u, (x, t) -> 0.0, 1]
        ),
    neumann = nothing,
)

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

elas = Elasticity(E = 70e9, ν = 0.3)

visco = Viscosity(η = x -> (x[2] >= Ly - H_bdt) ? 1e24 : 1e17)

VEP_rheology = Rheology(viscosity = visco, elasticity = elas)

### CLOCK ###
clock = Clock(tspan = (0.0, 50year), Δt = 0.1year, Δt_max = 20year)

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
output_name = "viscous_example"
VEP_outputs = Dict(
    :σxx => (r, s) -> s.σ[1, 1],
    :σyy => (r, s) -> s.σ[2, 2],
    :σzz => (r, s) -> s.σ[3, 3],
    :η => (r, s) -> r.viscosity.η,
    :p => (r, s) -> -1 / 3 * tr(s.σ),
)

VEP_ow = VTKOutputWriter(
    VEP_model,
    joinpath(output_path, output_name),
    VEP_outputs,
    interval = 1,
    force_path = true,
)

### SOLVE ###
@time modelsol, u = solve(VEP_model, output_writer = VEP_ow, log = true);
