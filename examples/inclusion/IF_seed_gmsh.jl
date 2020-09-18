using Rheologies ; const R = Rheologies

include(joinpath(@__DIR__,"sample_circ_inclusion.jl"))

output_path = @__DIR__

## SPATIAL DOMAIN ##
dim = 2 # spatial dimensions

    Lx = 1 # same as inclusion.msh mesh
    Ly = 2 # same as inclusion.msh mesh
    radius = Lx/40
    el_geom = Triangle # Quadrilateral or Triangle, equivalent to Cell{dim,nnodes,nfaces}

    # generate the grid and sets:
    meshfile = joinpath(output_path,"inclusion.msh")

    generate_mesh_inclusion(; Lx, Ly, radius ,lowres=0.03, highres=0.002, file=meshfile)
    grid = generate_grid(meshfile)
    #grid = generate_grid(el_geom, (nx, ny), corner1, Vec(corner1[1]+Lx, corner1[2]+Ly)) # can take 2 or 4 corners

    # Add facesets / nodesets / cellsets on which bc will be applied (except for the top, bottom, left and right boundary that are implemented by default)
    #addnodeset!(grid, "clamped", x -> (x[1] == corner_u[1] && 24.5<x[2]<=25.5));
    #addnodeset!(grid, "clamped", x -> ( (abs(x[1]-Lx/2)<0.01) & any(abs.(x[2] .- (0.0,Ly)).< 0.01) ) );
    addnodeset!(grid, "clamped", x -> ( (x[1] ≈ Lx/2) & any(x[2] .≈ (0.0,Ly)) ) );


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
    v_in = 1e-5
    bc = BoundaryConditions(dirichlet = Dict("top" => [:u, (x,t)-> -v_in*t, 2],
                                             "bottom" => [:u, (x,t)-> v_in*t, 2],
                                             "clamped" => [:u, (x,t)-> 0.0, 1]),
                            neumann   = nothing)#Dict("left" => (x,t)->Vec(1e3,0.0),"right" => (x,t)->Vec(-1e3,0.0)))

    bf = BodyForces([0.0,0.0])#-9.81*2700000

    ## MATERIAL RHEOLOGY :
    seed(x,x0,r,val_in::T,val_out::T) where {T<:Real} = (sqrt((x[1]-x0[1])^2 + (x[2]-x0[2])^2) <= r ?
                                       val_in : val_out)

    # Define seed location accordingly to the mesh seed
    x0 = Vec(Lx/2,Ly/2)

    elas = Elasticity(E = x -> seed(x,x0,radius,30e9,70e9),
                      ν = x -> seed(x,x0,radius,0.45,0.3) )

    Δσ = 5e3
    plas  = DruckerPrager(ϕ = 30.0,
                          C = x -> seed(x,x0,radius,1e6,2e6),
                          H = -5e7,
                          ηᵛᵖ = Δσ*Ly/(2*v_in))


    rheology = Rheology(elasticity = elas,
                        plasticity = plas )

    ### CLOCK ###
    clock = Clock(tspan = (0.0,40.0), Δt = 1.0, Δt_max = 10.0, Δt_fact_up = 1.2)

    ### SOLVER ###
    linsolver = BackslashSolver()
    #linsolver = MUMPS()
    # linsolver = ConjugateGradient(package = :IterativeSolvers,
    #                              (max_iter = 1000,))
    #linsolver = MINRES(package = :IterativeSolvers)
    nlsolver = NewtonRaphson(max_iter_number = 10, atol = 1e-5, linear_solver = linsolver)


    ### MODEL ###
    model = Model( grid = grid,
                   variables = variables,
                   quad_order = quad_order,
                   quad_type = quad_type,
                   bc = bc,
                   body_forces = bf,
                   rheology = rheology,
                   initial_state = nothing,
                   clock = clock,
                   solver = nlsolver )

    ### OUTPUT ###
    outputs = Dict(:σxx     => (r,s)-> s.σ[1,1],
                   :σyy     => (r,s)-> s.σ[2,2],
                   :σzz     => (r,s)-> s.σ[3,3],
                   :ϵxx     => (r,s)-> s.ϵ[1,1],
                   :ϵyy     => (r,s)-> s.ϵ[2,2],
                   :ϵzz     => (r,s)-> s.ϵ[3,3],
                   :ϵ_vol   => (r,s)-> tr(s.ϵ),
                   :gamma   => (r,s)-> sqrt(2* dev(s.ϵ) ⊡ dev(s.ϵ)),
                   :E       => (r,s)-> r.elasticity.E,
                   :nu      => (r,s)-> r.elasticity.ν,
                   :norm_ep => (r,s)-> norm(s.ϵᵖ),
                   :τ       => (r,s)-> R.get_τ(dev(s.σ),r.plasticity),
                   :p       => (r,s)-> -1/3 * tr(s.σ),
                   :τoverp  => (r,s)-> R.get_τ(dev(s.σ),r.plasticity)/(-1/3 * tr(s.σ))
                   )
    # mat_outputs = Dict(:sigma_xx     => (r,s)-> s.σ[1,1],
    #               :sigma_yy     => (r,s)-> s.σ[2,2],
    #               :sigma_zz     => (r,s)-> s.σ[3,3],
    #               :eps_xx     => (r,s)-> s.ϵ[1,1],
    #               :eps_yy     => (r,s)-> s.ϵ[2,2],
    #               :eps_zz     => (r,s)-> s.ϵ[3,3],
    #               :norm_ep    => (r,s)-> norm(s.ϵᵖ),
    #               :tau        => (r,s)-> R.get_τ(dev(s.σ),r.plasticity),
    #               :p          => (r,s)-> -1/3 * tr(s.σ),
    #               :tauoverp  => (r,s)-> R.get_τ(dev(s.σ),r.plasticity)/-(1/3 * tr(s.σ))
    #               )

    VTK_ow = VTKOutputWriter(model, joinpath(output_path,"inclusion_elastique_gmsh_EP"), outputs, interval = 5, force_path = true);
    #MAT_ow = MATOutputWriter(model, joinpath(path,"sawcut_EP_MAT"), mat_outputs, interval = 1, force_path = true);
    #JLD2_ow = JLD2OutputWriter(model, joinpath(path,"seed_EP_JLD2"), outputs, interval = 1, force_path = true);
    #ow = MixedOutputWriter(VTK_ow,JLD2_ow);
### SOLVE ###
    @time modelsol, u = solve(model, output_writer = VTK_ow ,log = true);
