# Rheologies

A [Julia](http://julialang.org) package for finite element based simulations of visco-elasto-plastic deformation. It extends the `JuAFEM` package, used as a toolbox for finite element modeling.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Leooop.github.io/Rheologies.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Leooop.github.io/Rheologies.jl/dev)
[//]: # [![Build Status](https://travis-ci.com/Leooop/Rheologies.jl.svg?branch=master)](https://travis-ci.com/Leooop/Rheologies.jl)
[//]: # [![Build Status](https://ci.appveyor.com/api/projects/status/github/Leooop/Rheologies.jl?svg=true)](https://ci.appveyor.com/project/Leooop/Rheologies-jl)
[//]: # [![Codecov](https://codecov.io/gh/Leooop/Rheologies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Leooop/Rheologies.jl)
[//]: # [![Coveralls](https://coveralls.io/repos/github/Leooop/Rheologies.jl/badge.svg?branch=master)](https://coveralls.io/github/Leooop/Rheologies.jl?branch=master)
[//]: # [![Build Status](https://api.cirrus-ci.com/github/Leooop/Rheologies.jl.svg)](https://cirrus-ci.com/github/Leooop/Rheologies.jl)

#### Author
- Léo Petit, École Normale Supérieure de Paris, France.

#### License

`Rheologies` is licensed under the [MIT license](./LICENSE.md).

#### Installation

`Rheologies` is not a registered package and can be installed from the package REPL with
```
pkg> add https://github.com/Leooop/Rheologies.jl.git
```
Requires Julia v1.3 or higher

#### Exemple

Elasto-plastic 2D plane strain deformation :

```julia
## TIME ##
clock = Clock(tspan = (0.0,10.0), Δt = 0.5) # in seconds

## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

# rectangular material discretization and geometry
nx = 40 # n elements along x
ny = 80 # n elements along y
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

## MATERIAL RHEOLOGY :
seed(x,x0,r,val_in::T,val_out::T) where {T<:Real} = (sqrt((x[1]-x0[1])^2 + (x[2]-x0[2])^2) <= r ?
                                   val_in : val_out)

elas = Elasticity(E = x -> seed(x,Vec(Lx/2,Ly/2),Lx/nx,30e9,70e9),
                  ν = 0.3)

plas  = DruckerPrager(ϕ = 30.0,
                      C = 1e6)

rheology = Rheology(elasticity = elast,
                    plasticity = plas)

### SOLVER ###
# Use Newton Raphson non linear iterations
nlsolver = NewtonRaphson(atol = 1e-5, linear_solver = BackslashSolver())

### MODEL ###
model = Model( grid = grid,
                       variables = variables,
                       quad_order = quad_order,
                       quad_type = quad_type,
                       bc = bc,
                       rheology = elastoplastic_rheology,
                       clock = clock,
                       solver = nlsolver )

### OUTPUTS ###
path = "path/to/folder"

# what to export :
outputs = Dict( :σxx     => (r,s)-> s.σ[1,1],
                :σyy     => (r,s)-> s.σ[2,2],
                :σzz     => (r,s)-> s.σ[3,3],
                :acum_ep => (r,s)-> norm(s.ϵ̅ᵖ))

ow = VTKOutputWriter(model, path, outputs, interval = 1)

### SOLVE ###
@time model_sol, u = solve(model ; output_writer = ow)
 ```
