### CONVENIENCE FUNCTIONS ###

# Tensorial calculus related functions
δ(i,j) = i == j ? 1.0 : 0.0
Isym_func(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
Isymdev_func(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
get_3D_func(i,j,tensor2D) = (i<=2) && (j<=2) ? tensor2D[i,j] : 0.0
