import Gmsh: gmsh

function generate_mesh_inclusion(; Lx=1.0, Ly=2.0 ,radius=Lx/40 ,lowres=0.02, highres=0.001, file)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")

    xmin, xmax = 0.0, Lx
    ymin, ymax = 0.0, Ly

    p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lowres)
    p2 = gmsh.model.geo.addPoint(xmax, ymin,  0, lowres)
    p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lowres)
    p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lowres)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    # Add two points that will be used as fixed in ux :
    # ptop = gmsh.model.geo.addPoint(0.5, 2, 0, lowres)
    # pbot = gmsh.model.geo.addPoint(0.5, 0, 0, lowres)
    # embed these points in the surface to force the mesh to define nodes at their location.
    # pbot = gmsh.model.geo.addPoint(xmax/2, ymin, 0, lowres)
    # ptop = gmsh.model.geo.addPoint(xmax/2, ymax, 0, lowres)

    #pmid = gmsh.model.geo.addPoint(xmax/2, ymax/2, 0, highres)

    #circ = gmsh.model.occ.addCircle(xmax/2, ymax/2, 0.0, xmax/20)
    pmid = gmsh.model.geo.addPoint(xmax/2, ymax/2, 0, 5*highres)
    pc1 = gmsh.model.geo.addPoint(xmax/2, ymax/2 - radius, 0, 2highres)
    pc2 = gmsh.model.geo.addPoint(xmax/2 + radius, ymax/2, 0, 2highres)
    pc3 = gmsh.model.geo.addPoint(xmax/2, ymax/2 + radius, 0, 2highres)
    pc4 = gmsh.model.geo.addPoint(xmax/2 - radius, ymax/2, 0, 2highres)
    arc1 = gmsh.model.geo.addCircleArc(pc1,pmid,pc2)
    arc2 = gmsh.model.geo.addCircleArc(pc2,pmid,pc3)
    arc3 = gmsh.model.geo.addCircleArc(pc3,pmid,pc4)
    arc4 = gmsh.model.geo.addCircleArc(pc4,pmid,pc1)
    cl2 = gmsh.model.geo.addCurveLoop([arc1,arc2,arc3,arc4])
    s1 = gmsh.model.geo.addPlaneSurface([cl1,cl2])
    s2 = gmsh.model.geo.addPlaneSurface([cl2])
    #gmsh.model.geo.synchronize()
    # gmsh.model.mesh.embed(0, [pbot], 1, l1)
    # gmsh.model.mesh.embed(0, [ptop], 1, l3)
    #gmsh.model.mesh.embed(2, [s2], 2, s1)

    #clamped = gmsh.model.addPhysicalGroup(0, [ptop, pbot])
    bottom = gmsh.model.addPhysicalGroup(1, [l1])
    right = gmsh.model.addPhysicalGroup(1, [l2])
    top = gmsh.model.addPhysicalGroup(1, [l3])
    left = gmsh.model.addPhysicalGroup(1, [l4])
    circle = gmsh.model.addPhysicalGroup(1, [cl2])
    domain1 = gmsh.model.addPhysicalGroup(2, [s1])
    domain2 = gmsh.model.addPhysicalGroup(2, [s2])

    #gmsh.model.setPhysicalName(0, clamped, "clamped")
    gmsh.model.setPhysicalName(1, bottom, "bottom")
    gmsh.model.setPhysicalName(1, right, "right")
    gmsh.model.setPhysicalName(1, top, "top")
    gmsh.model.setPhysicalName(1, left, "left")
    gmsh.model.setPhysicalName(1, circle, "circle")
    gmsh.model.setPhysicalName(2, domain1, "domain1")
    gmsh.model.setPhysicalName(2, domain2, "domain2")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    gmsh.finalize()

    return nothing
end
