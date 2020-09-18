import Gmsh: gmsh

function generate_mesh(x0,radius; Lx=80e3, Ly=30e3, H_bdt = 10e3, lowres=0.02, highres=0.001, file)

    @assert length(x0) == length(radius)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")

    xmin, xmax = 0.0, Lx
    ymin, ymax = 0.0, Ly
    L_bdt = Ly-H_bdt

    p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lowres)
    p2 = gmsh.model.geo.addPoint(xmax, ymin,  0, lowres)
    p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lowres)
    p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lowres)
    pbdt1 = gmsh.model.geo.addPoint(xmin, L_bdt, 0, lowres)
    pbdt2 = gmsh.model.geo.addPoint(xmax, L_bdt, 0, lowres)
    l_bot = gmsh.model.geo.addLine(p1, p2)
    lr = gmsh.model.geo.addLine(p2, p3)
    lr_bot = gmsh.model.geo.addLine(p2, pbdt2)
    lr_top = gmsh.model.geo.addLine(pbdt2,p3)
    l_top = gmsh.model.geo.addLine(p3, p4)
    ll = gmsh.model.geo.addLine(p4, p1)
    ll_top = gmsh.model.geo.addLine(p4, pbdt1)
    ll_bot = gmsh.model.geo.addLine(pbdt1,p1)
    l_bdt1 = gmsh.model.geo.addLine(pbdt1, pbdt2)
    l_bdt2 = gmsh.model.geo.addLine(pbdt2, pbdt1)
    cl_lcrust = gmsh.model.geo.addCurveLoop([l_bot, lr_bot, l_bdt2, ll_bot])
    cl_ucrust = gmsh.model.geo.addCurveLoop([l_bdt1, lr_top, l_top, ll_top])
    #cl_domain = gmsh.model.geo.addCurveLoop([l_bot,lr,l_top,ll])
    # Add two points that will be used as fixed in ux :
    # ptop = gmsh.model.geo.addPoint(0.5, 2, 0, lowres)
    # pbot = gmsh.model.geo.addPoint(0.5, 0, 0, lowres)
    # embed these points in the surface to force the mesh to define nodes at their location.
    # pbot = gmsh.model.geo.addPoint(xmax/2, ymin, 0, lowres)
    # ptop = gmsh.model.geo.addPoint(xmax/2, ymax, 0, lowres)

    #pmid = gmsh.model.geo.addPoint(xmax/2, ymax/2, 0, highres)
    #gmsh.finalize()

    #circ = gmsh.model.occ.addCircle(xmax/2, ymax/2, 0.0, xmax/20)
    T = typeof(p1)
    curved_loops = T[]
    curved_loops_phys = T[]
    surfaces = T[]
    surfaces_phys = T[]
    for i in eachindex(radius)
        x,y = x0[i]
        r = radius[i]
        pmid = gmsh.model.geo.addPoint(x, y, 0, 2highres)
        pc1 = gmsh.model.geo.addPoint(x, y-r, 0, 2highres)
        pc2 = gmsh.model.geo.addPoint(x+r, y, 0, 2highres)
        pc3 = gmsh.model.geo.addPoint(x, y+r, 0, 2highres)
        pc4 = gmsh.model.geo.addPoint(x-r, y, 0, 2highres)
        arc1 = gmsh.model.geo.addCircleArc(pc1,pmid,pc2)
        arc2 = gmsh.model.geo.addCircleArc(pc2,pmid,pc3)
        arc3 = gmsh.model.geo.addCircleArc(pc3,pmid,pc4)
        arc4 = gmsh.model.geo.addCircleArc(pc4,pmid,pc1)
        cl_seed = gmsh.model.geo.addCurveLoop([arc1,arc2,arc3,arc4])
        seed = gmsh.model.geo.addPlaneSurface([cl_seed])

        cl_seed_phys = gmsh.model.addPhysicalGroup(1, [cl_seed])
        seed_phys = gmsh.model.addPhysicalGroup(2, [seed])

        push!(curved_loops, cl_seed)
        push!(surfaces, seed)
        push!(curved_loops_phys, cl_seed_phys)
        push!(surfaces_phys, seed_phys)
    end
    # domain excluding seeds
    upper_crust = gmsh.model.geo.addPlaneSurface([cl_ucrust,curved_loops...])
    upper_crust_phys = gmsh.model.addPhysicalGroup(2, [upper_crust])

    lower_crust = gmsh.model.geo.addPlaneSurface([cl_lcrust])
    lower_crust_phys = gmsh.model.addPhysicalGroup(2, [lower_crust])

    # domain = gmsh.model.geo.addPlaneSurface([cl_domain,curved_loops...])
    # domain_phys = gmsh.model.addPhysicalGroup(2, [domain])

    # gmsh.model.geo.synchronize()
    # gmsh.model.mesh.embed(0, [l_bdt], 2, domain)
    # gmsh.model.mesh.embed(0, [ptop], 1, l3)
    #gmsh.model.mesh.embed(2, [s2], 2, s1)

    bottom = gmsh.model.addPhysicalGroup(1, [l_bot])
    right = gmsh.model.addPhysicalGroup(1, [lr_bot,lr_top])
    top = gmsh.model.addPhysicalGroup(1, [l_top])
    left = gmsh.model.addPhysicalGroup(1, [ll_bot,ll_top])
    #l_bdt = gmsh.model.addPhysicalGroup(1, [l_bdt1])



    #gmsh.model.setPhysicalName(0, clamped, "clamped")
    gmsh.model.setPhysicalName(1, bottom, "bottom")
    gmsh.model.setPhysicalName(1, right, "right")
    gmsh.model.setPhysicalName(1, top, "top")
    gmsh.model.setPhysicalName(1, left, "left")
    #gmsh.model.setPhysicalName(1, l_bdt, "bdt")
    gmsh.model.setPhysicalName(2, lower_crust_phys, "lower_crust")
    gmsh.model.setPhysicalName(2, upper_crust_phys, "upper_crust")
    #gmsh.model.setPhysicalName(2, domain_phys, "domain")
    for i in eachindex(surfaces_phys)
        gmsh.model.setPhysicalName(1, curved_loops_phys[i], "circle_seed_$i")
        gmsh.model.setPhysicalName(2, surfaces_phys[i], "seed_$i")
    end
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    gmsh.finalize()

    return nothing
end

function generate_mesh2(x0,radius; Lx=80e3, Ly=30e3, H_bdt = 10e3, lowres=0.02, highres=0.001, file)

    @assert length(x0) == length(radius)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")

    xmin, xmax = 0.0, Lx
    ymin, ymax = 0.0, Ly
    L_bdt = Ly-H_bdt

    p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lowres)
    p2 = gmsh.model.geo.addPoint(xmax, ymin,  0, lowres)
    p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lowres)
    p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lowres)
    pbdt1 = gmsh.model.geo.addPoint(xmin, L_bdt, 0, lowres)
    pbdt2 = gmsh.model.geo.addPoint(xmax, L_bdt, 0, lowres)
    l_bot = gmsh.model.geo.addLine(p2,p1)
    lr = gmsh.model.geo.addLine(p2, p3)
    lr_bot = gmsh.model.geo.addLine(pbdt2,p2)
    lr_top = gmsh.model.geo.addLine(pbdt2,p3)
    l_top = gmsh.model.geo.addLine(p3, p4)
    ll = gmsh.model.geo.addLine(p4, p1)
    ll_top = gmsh.model.geo.addLine(p4, pbdt1)
    ll_bot = gmsh.model.geo.addLine(p1, pbdt1)
    l_bdt1 = gmsh.model.geo.addLine(pbdt1, pbdt2)
    l_bdt2 = gmsh.model.geo.addLine(pbdt2, pbdt1)
    # cl_lcrust = gmsh.model.geo.addCurveLoop([l_bdt1, lr_bot, l_bot, ll_bot])
    # cl_ucrust = gmsh.model.geo.addCurveLoop([l_bdt1, lr_top, l_top, ll_top])
    cl_domain = gmsh.model.geo.addCurveLoop([l_bot,lr,l_top,ll])
    # Add two points that will be used as fixed in ux :
    # ptop = gmsh.model.geo.addPoint(0.5, 2, 0, lowres)
    # pbot = gmsh.model.geo.addPoint(0.5, 0, 0, lowres)
    # embed these points in the surface to force the mesh to define nodes at their location.
    # pbot = gmsh.model.geo.addPoint(xmax/2, ymin, 0, lowres)
    # ptop = gmsh.model.geo.addPoint(xmax/2, ymax, 0, lowres)

    #pmid = gmsh.model.geo.addPoint(xmax/2, ymax/2, 0, highres)
    #gmsh.finalize()

    #circ = gmsh.model.occ.addCircle(xmax/2, ymax/2, 0.0, xmax/20)
    T = typeof(p1)
    curved_loops = T[]
    curved_loops_phys = T[]
    surfaces = T[]
    surfaces_phys = T[]
    for i in eachindex(radius)
        x,y = x0[i]
        r = radius[i]
        pmid = gmsh.model.geo.addPoint(x, y, 0, 2highres)
        pc1 = gmsh.model.geo.addPoint(x, y-r, 0, 2highres)
        pc2 = gmsh.model.geo.addPoint(x+r, y, 0, 2highres)
        pc3 = gmsh.model.geo.addPoint(x, y+r, 0, 2highres)
        pc4 = gmsh.model.geo.addPoint(x-r, y, 0, 2highres)
        arc1 = gmsh.model.geo.addCircleArc(pc1,pmid,pc2)
        arc2 = gmsh.model.geo.addCircleArc(pc2,pmid,pc3)
        arc3 = gmsh.model.geo.addCircleArc(pc3,pmid,pc4)
        arc4 = gmsh.model.geo.addCircleArc(pc4,pmid,pc1)
        cl_seed = gmsh.model.geo.addCurveLoop([arc1,arc2,arc3,arc4])
        seed = gmsh.model.geo.addPlaneSurface([cl_seed])

        cl_seed_phys = gmsh.model.addPhysicalGroup(1, [cl_seed])
        seed_phys = gmsh.model.addPhysicalGroup(2, [seed])

        push!(curved_loops, cl_seed)
        push!(surfaces, seed)
        push!(curved_loops_phys, cl_seed_phys)
        push!(surfaces_phys, seed_phys)
    end
    # domain excluding seeds
    upper_crust = gmsh.model.geo.addPlaneSurface([cl_ucrust,curved_loops...])
    upper_crust_phys = gmsh.model.addPhysicalGroup(2, [upper_crust])

    lower_crust = gmsh.model.geo.addPlaneSurface([cl_lcrust])
    lower_crust_phys = gmsh.model.addPhysicalGroup(2, [lower_crust])

    # domain = gmsh.model.geo.addPlaneSurface([cl_domain,curved_loops...])
    # domain_phys = gmsh.model.addPhysicalGroup(2, [domain])

    # gmsh.model.geo.synchronize()
    # gmsh.model.mesh.embed(0, [l_bdt], 2, domain)
    # gmsh.model.mesh.embed(0, [ptop], 1, l3)
    #gmsh.model.mesh.embed(2, [s2], 2, s1)

    bottom = gmsh.model.addPhysicalGroup(1, [l_bot])
    right = gmsh.model.addPhysicalGroup(1, [lr_bot,lr_top])
    top = gmsh.model.addPhysicalGroup(1, [l_top])
    left = gmsh.model.addPhysicalGroup(1, [ll_bot,ll_top])
    #l_bdt = gmsh.model.addPhysicalGroup(1, [l_bdt1])



    #gmsh.model.setPhysicalName(0, clamped, "clamped")
    gmsh.model.setPhysicalName(1, bottom, "bottom")
    gmsh.model.setPhysicalName(1, right, "right")
    gmsh.model.setPhysicalName(1, top, "top")
    gmsh.model.setPhysicalName(1, left, "left")
    #gmsh.model.setPhysicalName(1, l_bdt, "bdt")
    gmsh.model.setPhysicalName(2, lower_crust_phys, "lower_crust")
    gmsh.model.setPhysicalName(2, upper_crust_phys, "upper_crust")
    #gmsh.model.setPhysicalName(2, domain_phys, "domain")
    for i in eachindex(surfaces_phys)
        gmsh.model.setPhysicalName(1, curved_loops_phys[i], "circle_seed_$i")
        gmsh.model.setPhysicalName(2, surfaces_phys[i], "seed_$i")
    end
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    gmsh.finalize()

    return nothing
end

function generate_mesh3(x0,radius; Lx=80e3, Ly=30e3, H_bdt = 10e3, lowres=0.02, highres=0.001, el_type = Triangle, file)

    @assert length(x0) == length(radius)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")

    xmin, xmax = 0.0, Lx
    ymin, ymax = 0.0, Ly
    L_bdt = Ly-H_bdt

    p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lowres)
    p2 = gmsh.model.geo.addPoint(xmax, ymin,  0, lowres)
    p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lowres)
    p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lowres)
    pbdt1 = gmsh.model.geo.addPoint(xmin-1, L_bdt, 0, lowres)
    pbdt2 = gmsh.model.geo.addPoint(xmax+1, L_bdt, 0, lowres)
    l_bot = gmsh.model.geo.addLine(p1, p2)
    lr = gmsh.model.geo.addLine(p2, p3)
    l_top = gmsh.model.geo.addLine(p3, p4)
    ll = gmsh.model.geo.addLine(p4, p1)
    l_bdt = gmsh.model.geo.addLine(pbdt1, pbdt2)
    # cl_lcrust = gmsh.model.geo.addCurveLoop([l_bot, lr_bot, l_bdt2, ll_bot])
    # cl_ucrust = gmsh.model.geo.addCurveLoop([l_bdt1, lr_top, l_top, ll_top])
    cl_domain = gmsh.model.geo.addCurveLoop([l_bot,lr,l_top,ll])
    domain = gmsh.model.geo.addPlaneSurface([cl_domain])
    # Add two points that will be used as fixed in ux :
    # ptop = gmsh.model.geo.addPoint(0.5, 2, 0, lowres)
    # pbot = gmsh.model.geo.addPoint(0.5, 0, 0, lowres)
    # embed these points in the surface to force the mesh to define nodes at their location.
    # pbot = gmsh.model.geo.addPoint(xmax/2, ymin, 0, lowres)
    # ptop = gmsh.model.geo.addPoint(xmax/2, ymax, 0, lowres)

    #pmid = gmsh.model.geo.addPoint(xmax/2, ymax/2, 0, highres)
    #gmsh.finalize()

    #circ = gmsh.model.occ.addCircle(xmax/2, ymax/2, 0.0, xmax/20)
    T = typeof(p1)
    curved_loops = T[]
    curved_loops_phys = T[]
    surfaces = T[]
    surfaces_phys = T[]
    for i in eachindex(radius)
        x,y = x0[i]
        r = radius[i]
        pmid = gmsh.model.geo.addPoint(x, y, 0, 2highres)
        pc1 = gmsh.model.geo.addPoint(x, y-r, 0, 2highres)
        pc2 = gmsh.model.geo.addPoint(x+r, y, 0, 2highres)
        pc3 = gmsh.model.geo.addPoint(x, y+r, 0, 2highres)
        pc4 = gmsh.model.geo.addPoint(x-r, y, 0, 2highres)
        arc1 = gmsh.model.geo.addCircleArc(pc1,pmid,pc2)
        arc2 = gmsh.model.geo.addCircleArc(pc2,pmid,pc3)
        arc3 = gmsh.model.geo.addCircleArc(pc3,pmid,pc4)
        arc4 = gmsh.model.geo.addCircleArc(pc4,pmid,pc1)
        cl_seed = gmsh.model.geo.addCurveLoop([arc1,arc2,arc3,arc4])
        seed = gmsh.model.geo.addPlaneSurface([cl_seed])

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0, [pc1,pc2,pc3,pc4], 2, domain)
        gmsh.model.mesh.embed(1, [arc1,arc2,arc3,arc4], 2, domain)
        #cl_seed_phys = gmsh.model.addPhysicalGroup(1, [cl_seed])
        #seed_phys = gmsh.model.addPhysicalGroup(2, [seed])

        push!(curved_loops, cl_seed)
        push!(surfaces, seed)
        #push!(curved_loops_phys, cl_seed_phys)
        #push!(surfaces_phys, seed_phys)
    end
    # domain excluding seeds
    # upper_crust = gmsh.model.geo.addPlaneSurface([cl_ucrust,curved_loops...])
    # upper_crust_phys = gmsh.model.addPhysicalGroup(2, [upper_crust])
    #
    # lower_crust = gmsh.model.geo.addPlaneSurface([cl_lcrust])
    # lower_crust_phys = gmsh.model.addPhysicalGroup(2, [lower_crust])


    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, [l_bdt], 2, domain)

    # gmsh.model.geo.synchronize()
    #gmsh.model.mesh.embed(2, [s2], 2, s1)

    bottom = gmsh.model.addPhysicalGroup(1, [l_bot])
    right = gmsh.model.addPhysicalGroup(1, [lr])
    top = gmsh.model.addPhysicalGroup(1, [l_top])
    left = gmsh.model.addPhysicalGroup(1, [ll])
    domain_phys = gmsh.model.addPhysicalGroup(2, [domain])
    #l_bdt = gmsh.model.addPhysicalGroup(1, [l_bdt1])



    #gmsh.model.setPhysicalName(0, clamped, "clamped")
    gmsh.model.setPhysicalName(1, bottom, "bottom")
    gmsh.model.setPhysicalName(1, right, "right")
    gmsh.model.setPhysicalName(1, top, "top")
    gmsh.model.setPhysicalName(1, left, "left")
    #gmsh.model.setPhysicalName(1, l_bdt, "bdt")
    # gmsh.model.setPhysicalName(2, lower_crust_phys, "lower_crust")
    # gmsh.model.setPhysicalName(2, upper_crust_phys, "upper_crust")
    gmsh.model.setPhysicalName(2, domain_phys, "domain")
    for i in eachindex(surfaces_phys)
        gmsh.model.setPhysicalName(1, curved_loops_phys[i], "circle_seed_$i")
        gmsh.model.setPhysicalName(2, surfaces_phys[i], "seed_$i")
    end
    gmsh.model.geo.synchronize()
    #gmsh.option.setNumber("Mesh.Algorithm", 2)
    #gmsh.option.setNumber("Mesh.RecombineAll", 1)
    if el_type == Quadrilateral
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
        gmsh.model.mesh.setRecombine(2, domain)
    end
    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    gmsh.finalize()

    return nothing
end

function t15()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Copied from t1.py...
    lc = 1e-2
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)
    gmsh.model.geo.addPoint(0, .3, 0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(3, 2, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    # We change the mesh size to generate a coarser mesh
    lc = lc * 4
    gmsh.model.geo.mesh.setSize([(0, 1), (0, 2), (0, 3), (0, 4)], lc)

    # We define a new point
    gmsh.model.geo.addPoint(0.02, 0.02, 0., lc, 5)

    # We have to synchronize before embedding entites:
    gmsh.model.geo.synchronize()

    # One can force this point to be included ("embedded") in the 2D mesh, using the
    # `embed()' function:
    gmsh.model.mesh.embed(0, [5], 2, 1)

    # In the same way, one can use `embed()' to force a curve to be embedded in the
    # 2D mesh:
    gmsh.model.geo.addPoint(0.02, 0.12, 0., lc, 6)
    gmsh.model.geo.addPoint(0.04, 0.18, 0., lc, 7)
    gmsh.model.geo.addLine(6, 7, 5)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, [5], 2, 1)

    # Points and curves can also be embedded in volumes
    gmsh.model.geo.extrude([(2, 1)], 0, 0, 0.1)

    p = gmsh.model.geo.addPoint(0.07, 0.15, 0.025, lc)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, [p], 3, 1)

    gmsh.model.geo.addPoint(0.025, 0.15, 0.025, lc, p + 1)
    l = gmsh.model.geo.addLine(7, p + 1)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, [l], 3, 1)

    # Finally, we can also embed a surface in a volume:
    gmsh.model.geo.addPoint(0.02, 0.12, 0.05, lc, p + 2)
    gmsh.model.geo.addPoint(0.04, 0.12, 0.05, lc, p + 3)
    gmsh.model.geo.addPoint(0.04, 0.18, 0.05, lc, p + 4)
    gmsh.model.geo.addPoint(0.02, 0.18, 0.05, lc, p + 5)

    gmsh.model.geo.addLine(p + 2, p + 3, l + 1)
    gmsh.model.geo.addLine(p + 3, p + 4, l + 2)
    gmsh.model.geo.addLine(p + 4, p + 5, l + 3)
    gmsh.model.geo.addLine(p + 5, p + 2, l + 4)

    ll = gmsh.model.geo.addCurveLoop([l + 1, l + 2, l + 3, l + 4])
    s = gmsh.model.geo.addPlaneSurface([ll])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(2, [s], 3, 1)

    # Note that with the OpenCASCADE kernel (see `t16.py'), when the `fragment()'
    # function is applied to entities of different dimensions, the lower dimensional
    # entities will be autmatically embedded in the higher dimensional entities if
    # necessary.

    gmsh.model.mesh.generate(3)

    gmsh.write("t15.msh")

    gmsh.finalize()

end

function test_embed()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Copied from t1.py...
    lc = 1e-2
    xmax, ymax = 1.0, 2.0
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(xmax, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(xmax, ymax, 0, lc, 3)
    gmsh.model.geo.addPoint(0, ymax, 0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(3, 2, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    # In the same way, one can use `embed()' to force a curve to be embedded in the
    # 2D mesh:
    gmsh.model.geo.addPoint(-0.0001, ymax/2, 0., lc, 5)
    gmsh.model.geo.addPoint(xmax+0.0001, ymax/2, 0., lc, 6)
    gmsh.model.geo.addLine(5, 6, 5)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, [5], 2, 1)

    gmsh.model.mesh.generate(2)

    gmsh.write("test_embed.msh")

    gmsh.finalize()

end

function test_embed2()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Copied from t1.py...
    lc = 1e-2
    xmax, ymax = 1.0, 2.0
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(xmax, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(xmax, ymax, 0, lc, 3)
    gmsh.model.geo.addPoint(0, ymax, 0, lc, 4)
    gmsh.model.geo.addPoint(0.0001, ymax/2, 0., lc, 5)
    gmsh.model.geo.addPoint(xmax-0.0001, ymax/2, 0., lc, 6)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 6, 2)
    gmsh.model.geo.addLine(6, 5, 3)
    gmsh.model.geo.addLine(5, 1, 4)

    gmsh.model.geo.addLine(6, 3, 5)
    gmsh.model.geo.addLine(3, 4, 6)
    gmsh.model.geo.addLine(4, 5, 7)

    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addCurveLoop([-3, 5, 6, 7], 2)
    gmsh.model.geo.addPlaneSurface([1,2], 1)

    # In the same way, one can use `embed()' to force a curve to be embedded in the
    # 2D mesh:
    # gmsh.model.geo.addPoint(0.0001, ymax/2, 0., lc, 5)
    # gmsh.model.geo.addPoint(xmax-0.0001, ymax/2, 0., lc, 6)
    # gmsh.model.geo.addLine(5, 6, 5)
    #
    # gmsh.model.geo.synchronize()
    # gmsh.model.mesh.embed(1, [5], 2, 1)

    gmsh.model.mesh.generate(2)

    gmsh.write("test_embed2.msh")

    gmsh.finalize()

end
