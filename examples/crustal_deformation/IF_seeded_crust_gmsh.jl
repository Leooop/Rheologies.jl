using Rheologies ; const R = Rheologies

include(joinpath(@__DIR__,"seeded_crust_gmsh.jl"))

output_path = @__DIR__
## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

Lx = 40e3
Ly = 25e3
H_bdt = 10e3
n_seeds = 2
radius_bounds = (300, 400)
el_geom = Triangle # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

lowres=500
highres=40 # near seeds

radius = rand(n_seeds).*(radius_bounds[2]-radius_bounds[1]) .+ radius_bounds[1] # random radius from 1 to 3 km
x0 = [[rand()*Lx , rand()*(H_bdt-2radius[i])+(Ly-H_bdt)+radius[i]] for i in 1:n_seeds]

push!(radius,200.0)
push!(x0, [Lx/2,Ly-radius[end]-100])
#push!(x0, [Lx/2,20e3])

# generate the grid and sets:
meshfile = joinpath(output_path,"crust.msh")
generate_mesh3(x0, radius ; Lx, Ly, H_bdt, lowres, highres, el_type = el_geom, file = meshfile)
grid = generate_grid(meshfile) # can take 2 or 4 corners
# Add facesets / nodesets / cellsets on which bc will be applied (except for the top, bottom, left and right boundary that are implemented by default)
#addnodeset!(grid, "clamped", x -> (x[1] == corner_u[1] && 24.5<x[2]<=25.5));
addnodeset!(grid, "clamped", x -> ( any(x[1] .≈ (0.0,Lx)) & (x[2] ≈ Ly-H_bdt) ));
#addnodeset!(grid, "clamped_down", x -> ( (x[1] == 0.0) & (x[2] in (0.0,Ly) ) ));


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
v_out = 0.01/year # 1cm/year
v_in = 2*Ly*v_out/Lx
bc = BoundaryConditions(dirichlet = Dict("left" => [:u, (x,t)-> v_out*t, 1],
                                         "right" => [:u, (x,t)-> -v_out*t, 1],
                                         #"top" => [:u, (x,t)-> Vec(0.0,0.0), [1,2]],
                                         "bottom" => [:u, (x,t)-> 0.0, 2]
                                         #"clamped" => [:u, (x,t)-> 0.0, 2]
                                         ),
                        neumann   = nothing)#Dict("top" => (x,t)->Vec(0.0,-1e5),
                                         #"bottom" => (x,t)->Vec(0.0,1e5)) )

bf = BodyForces([0.0,-9.81*2700])#-9.81*2700])

## MATERIAL RHEOLOGY :

# define a rheology distribution function ...
function multi_seeds(x,x0,radius,val_in,val_out)
    for i in eachindex(radius)
        (sqrt((x[1]-x0[i][1])^2 + (x[2]-x0[i][2])^2) <= radius[i]) && return val_in
    end
    return val_out
end
multi_seeds(x,val_in,val_out) = multi_seeds(x,x0,radius,val_in,val_out)
# and use it ...
elas = Elasticity(E = x -> multi_seeds(x,30e9,70e9),
                  ν = x -> multi_seeds(x,0.35,0.3) )

visco = Viscosity(η = x -> (x[2] >= Ly-H_bdt) ? 1e24 : 1e18)

Δσ = 2e4
plas  = DruckerPrager(ϕ = 30.0,
                      C = x -> multi_seeds(x,1e6,1e5),
                      H = -1e6,
                      ηᵛᵖ = Δσ*Lx/(2*v_out) )

VEP_rheology = Rheology(viscosity = visco,
                        elasticity = elas,
                        plasticity = plas )


# function sigma_func(x,bf)
#     δ = R.δ
#     p = bf.components[2] * Ly-x[2]
#     σ(i,j) = p*δ(i,j)
#     return SymmetricTensor{2,3,Float64}(σ)
# end
# sigma_func(x) = sigma_func(x,bf)
#
# function epsilon_func(x, sigma_func, bf, elast)
#     G = elast.G(x)
#     ν = elast.ν(x)
#     λ = R.λ_from_Eν(G,ν)
#     Dᵉ = R.get_elastic_stiffness_tensor(G,λ)
#     σ = sigma_func(x,bf)
#     return Dᵉ ⊡ σ
# end
#
# initial_state = Dict( :σ => x->sigma_func(x,bf),
#                       :ϵ => x->epsilon_func(x,sigma_func,bf,elas) )
### CLOCK ###
clock = Clock(tspan = (0.0,500year), Δt = 1year, Δt_max = 20year)

### SOLVER ###
linsolver = BackslashSolver()
nlsolver = NewtonRaphson(max_iter_number = 15, atol = 10, linear_solver = linsolver)


### MODEL ###
VEP_model = Model( grid = grid,
                   variables = variables,
                   quad_order = quad_order,
                   quad_type = quad_type,
                   bc = bc,
                   body_forces = bf,
                   rheology = VEP_rheology,
                   initial_state = nothing,
                   clock = clock,
                   solver = nlsolver )

### OUTPUT ###
output_name = "VEP_crust_gmsh"
VEP_outputs = Dict( :σxx     => (r,s)-> s.σ[1,1],
                    :σyy     => (r,s)-> s.σ[2,2],
                    :σzz     => (r,s)-> s.σ[3,3],
                    :acum_ep => (r,s)-> s.ϵ̅ᵖ,
                    :gamma   => (r,s)-> sqrt(2* dev(s.ϵ) ⊡ dev(s.ϵ)),
                    #:η => (r,s)-> r.viscosity.η,
                    #:C => (r,s)-> r.plasticity.C,
                    #:E => (r,s)-> r.elasticity.E,
                    :τ       => (r,s)-> R.get_τ(dev(s.σ),r.plasticity),
                    :p       => (r,s)-> -1/3 * tr(s.σ),
                    :τoverp  => (r,s)-> R.get_τ(dev(s.σ),r.plasticity)/(-1/3 * tr(s.σ)) )

VEP_ow = VTKOutputWriter(VEP_model, joinpath(output_path, output_name), VEP_outputs, interval = 1, force_path = true)

### SOLVE ###
@time modelsol, u = solve(VEP_model, output_writer = VEP_ow ,log = true);
