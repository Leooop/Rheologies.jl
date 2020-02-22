

interp_geom = default_interpolation(el_geom) # By default the interpolation is quadratic only if there are nodes in edges centers.
interp_vars = get_variables_interpolation(variables, var_interp_order, el_geom)
qr, qr_face = get_quadrature_rules(quad_order, :legendre, el_geom)

dh = create_dofhandler(grid, interp_vars)
bc = create_dirichlet_bc(dh, bc_dirichlet)

# cellvalues
cellvalues_u, cellvalues_p, facevalues_u = create_values(qr, qr_face, interp_geom, interp_vars...)

## TODO : go on !!!
