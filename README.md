# Rheologies

A [Julia](http://julialang.org) package for finite element based simulations of visco-elasto-plastic deformation. It extends the `JuAFEM` package, used as a toolbox for finite element modeling.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Leooop.github.io/Rheologies.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Leooop.github.io/Rheologies.jl/dev)
<!---
[![Build Status](https://travis-ci.com/Leooop/Rheologies.jl.svg?branch=master)](https://travis-ci.com/Leooop/Rheologies.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Leooop/Rheologies.jl?svg=true)](https://ci.appveyor.com/project/Leooop/Rheologies-jl)
[![Codecov](https://codecov.io/gh/Leooop/Rheologies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Leooop/Rheologies.jl)
[![Coveralls](https://coveralls.io/repos/github/Leooop/Rheologies.jl/badge.svg?branch=master)](https://coveralls.io/github/Leooop/Rheologies.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/Leooop/Rheologies.jl.svg)](https://cirrus-ci.com/github/Leooop/Rheologies.jl)
-->


#### Author
- Léo Petit, École Normale Supérieure de Paris, France.

#### License

`Rheologies` is licensed under the [MIT license](./LICENSE.md).

#### Installation

`Rheologies` is not a registered package and can be installed from the package REPL with
```julia
pkg> add https://github.com/Leooop/Rheologies.jl.git
```
or similarly with
```julia
pkg> using Pkg ; Pkg.add("https://github.com/Leooop/Rheologies.jl.git")
```
Requires Julia v1.3 or higher

#### Examples

In addition to these examples, most exported types and functions are documented. Type
`?function_name` or `?type_name`

Elasto-plastic 2D plane strain deformation :

```julia
using Rheologies

## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

# rectangular material discretization and geometry
nx = 80 # n elements along x
ny = 160 # n elements along y
Lx = 0.05 # length along x
Ly = 0.1 # length along y
el_geom = Quadrilateral # element geometry (Quadrilateral or Triangle)

# generate the grid and sets:
grid = generate_grid(el_geom, (nx, ny), Vec(0.0,0.0), Vec(Lx,Ly)) # can take 2 or 4 corners
# Add facesets / nodesets / cellsets on which Dirichlet bc will be applied (top, bottom, left and right boundary are implemented by default)
addnodeset!(grid, "clamped", x -> ( (x[1] == Lx/2) & (x[2] in (0.0, Ly)) ) ); # middle of the top and bottom boundaries

## PRIMITIVE VARIABLES AND INTERPOLATIONS ##
variables = PrimitiveVariables{1}((:u,), (2,), el_geom) # {1} variable :u with second order interpolation on el_geom

## QUADRATURE ##
quad_order = 2 # quadrature order
quad_type = :legendre # quadrature rule (:legendre or :lobatto)

## BOUNDARY CONDITIONS ##
## DIRICHLET :
    # KEY : set on which to apply bc
    # VALUES :
         # 1) constrained variable
         # 2) constraint function of (x,t)
         # 3) constrained components of the variable
## NEUMANN :
    # KEY : set on which to apply neumann bc
    # VALUES : traction as a function of (x,t)

v_in = 1e-6 # inflow velocity
# Apply inward displacement to top and bottom boundaries and fix "clamped" set to prevent rigid body motion
bc = BoundaryConditions(dirichlet = Dict("top" => [:u, (x,t)-> -v_in*t, 2],
                                         "bottom" => [:u, (x,t)-> v_in*t, 2],
                                         "clamped" => [:u, (x,t)-> 0.0, 1]),
                        neumann   = nothing)

# BODY FORCES (not stable for now)
bf = BodyForces([0.0,0.0]) #-9.81*2700

## MATERIAL RHEOLOGY :
seed(x,x0,r,val_in::T,val_out::T) where {T<:Real} = (sqrt((x[1]-x0[1])^2 + (x[2]-x0[2])^2) <= r ?
                                   val_in : val_out)

elas = Elasticity(E = x -> seed(x,Vec(Lx/2,Ly/2),2Lx/nx,30e9,70e9),
                  ν = 0.3)

Δσ = 2e4
plas  = DruckerPrager(ϕ = 30.0,
                      C = 1e6,
                      H = -5e7,
                      ηᵛᵖ = Δσ*Ly/(2*v_in) )

rheology = Rheology(elasticity = elas,
                    plasticity = plas)


## TIME DOMAIN ##
clock = Clock(tspan = (0.0,20.0), Δt = 2) # in seconds

### SOLVER ###
# Use Newton Raphson non linear iterations
nlsolver = NewtonRaphson(atol = 1e-5, linear_solver = BackslashSolver())

### MODEL ###
model = Model( grid = grid,
               variables = variables,
               quad_order = quad_order,
               quad_type = quad_type,
               bc = bc,
               body_forces = bf,
               rheology = rheology,
               clock = clock,
               solver = nlsolver )

### OUTPUTS ###
path = "path/to/folder"

# what to export :
#key : name of the field
#value : function of (r,s) for rheology instance and state instance, state object contains σ and ϵ (+ ϵᵖ and ϵ̅ᵖ for plastic material)
outputs = Dict( :σxx     => (r,s)-> s.σ[1,1],
                :σyy     => (r,s)-> s.σ[2,2],
                :σzz     => (r,s)-> s.σ[3,3],
                :ϵxx     => (r,s)-> s.ϵ[1,1],
                :ϵyy     => (r,s)-> s.ϵ[2,2],
                :ϵzz     => (r,s)-> s.ϵ[3,3],
                :ϵ_dev_inv2 => (r,s) -> sqrt(dev(s.ϵ)⊡dev(s.ϵ)),
                :acum_ep => (r,s)-> s.ϵ̅ᵖ )

ow = VTKOutputWriter(model, path, outputs, interval = 1) # every `interval` iteration (can also be every `frequency` seconds if `frequency` keyword is used instead)

### SOLVE ###
@time model_sol, u = solve(model ; output_writer = ow, log = true) # log enables a performance evaluation of the simulation
 ```
![second deviatoric strain invariant](https://github.com/Leooop/Rheologies.jl/blob/master/sample_pic.png)


Visco-elasto-plastic 2D plane strain simulation of pre-weakened faulting in the upper crust (high viscosity) overlying lower crust (lower viscosity)

```julia
using Rheologies

## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

nx = 240 # n elements along x
ny = 160 # n elements along y
Lx = 30e3
Ly = 20e3
el_geom = Quadrilateral # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

# generate the grid and sets:
grid = generate_grid(el_geom, (nx, ny), Vec(0.0,-Ly), Vec(Lx,0.0)) # can take 2 or 4 corners
# Add facesets / nodesets / cellsets on which bc will be applied (except for the top, bottom, left and right boundary that are implemented by default)
#addnodeset!(grid, "clamped", x -> (x[1] == corner_u[1] && 24.5<x[2]<=25.5));
addnodeset!(grid, "clamped", x -> ( (x[1] in (0.0, Lx)) & (x[2] == 0.0) ) );


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
year = 3600*24*365.0
v_in = 0.01/year # 1cm/year
bc = BoundaryConditions(dirichlet = Dict("left" => [:u, (x,t)-> -v_in*t, 1],
                                         "right" => [:u, (x,t)-> v_in*t, 1],
                                         "clamped" => [:u, (x,t)-> 0.0, 2]),
                        neumann   = nothing)#Dict("top" => (x,t)->Vec(0.0,-1e5),
                                         #"bottom" => (x,t)->Vec(0.0,1e5)) )

bf = BodyForces([0.0,0.0])#-9.81*2700])

## MATERIAL RHEOLOGY :

# define a rheology distribution function ...
function fault2D(x,x0,angle,y_dist,y_max,val_in,val_out)
    if (abs(x[2] + (x0-x[1])*tand(angle)) <= y_dist) & (x[2] >= y_max)
        return val_in
    else
        return val_out
    end
end
fault2D(x,val_in,val_out) = fault2D(x,5Lx/8,60.0,Ly/ny,-Ly/2,val_in,val_out)

# and use it ...
elas = Elasticity(E = 70e9,
                  ν = 0.3)

visco = Viscosity(η = x -> x[2] >= -Ly/2 ? 1e24 : 1e19)

Δσ = 1e3#2e4
plas  = DruckerPrager(ϕ = 30.0,
                      C = x->fault2D(x,1e5,1.5e5),
                      H = x->fault2D(x,-1e7,0.0),
                      ηᵛᵖ = Δσ*Lx/(2*v_in) )

VEP_rheology = Rheology(viscosity = visco,
                        elasticity = elas,
                        plasticity = plas )

### CLOCK ###
clock = Clock(tspan = (0.0,20year), Δt = 1year, Δt_max = 20year)

### SOLVER ###
linsolver = BackslashSolver()
nlsolver = NewtonRaphson(max_iter_number = 20, atol = 1e-2, linear_solver = linsolver)


### MODEL ###
VEP_model = Model( grid = grid,
                   variables = variables,
                   quad_order = quad_order,
                   quad_type = quad_type,
                   bc = bc,
                   body_forces = bf,
                   rheology = VEP_rheology,
                   clock = clock,
                   solver = nlsolver )

### OUTPUT ###
path = "path/to/folder"
VEP_outputs = Dict( :σxx     => (r,s)-> s.σ[1,1],
                    :σyy     => (r,s)-> s.σ[2,2],
                    :σzz     => (r,s)-> s.σ[3,3],
                    :acum_ep => (r,s)-> s.ϵ̅ᵖ,
                    :η => (r,s)-> r.viscosity.η,
                    :C => (r,s)-> r.plasticity.C,
                    :E => (r,s)-> r.elasticity.E )

VEP_ow = VTKOutputWriter(VEP_model, path, VEP_outputs, interval = 1)

### SOLVE ###
@time modelsol, u = solve(VEP_model, output_writer = VEP_ow ,log = true);

 ```

![accumulated plastic strain and displacement field](https://github.com/Leooop/Rheologies.jl/blob/master/ep_pic.png)
