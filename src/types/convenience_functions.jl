### CONVENIENCE FUNCTIONS ###

# Tensorial calculus related functions
δ(i,j) = i == j ? 1.0 : 0.0
Isym_func(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
Isymdev_func(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
get_3D_func(i,j,tensor2D) = (i<=2) && (j<=2) ? tensor2D[i,j] : 0.0

# elastic moduli relations :
E_from_Gν(G,ν) = 2G*(1 + ν)
E_from_KG(K,G) = 9*K*G / (3K + G)
ν_from_KG(K,G) = (3K - 2G) / (2*(3K + G))
G_from_Eν(E,ν) = E / 2(1 + ν)
G_from_Kν(K,ν) = 3K*(1 - 2ν) / (2*(1+ν))
G_from_λν(λ,ν) = λ*(1 - 2ν) / (2ν)
G_from_EK(E,K) = 3K*E / (9K - E)
K_from_Eν(E,ν) = E / 3(1 - 2ν)
K_from_Gν(G,ν) = 2G*(1 + ν) / (3*(1 - 2ν))
λ_from_Eν(E,ν) = E*ν / ((1+ν)*(1-2ν))
λ_from_KG(K,G) = K - (2/3)*G
λ_from_Gν(G,ν) = 2G*ν / (1 - 2ν)

Dᵉ_func(i,j,k,l,G,λ) = λ*(δ(i,j)*δ(k,l)) + G*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
Cᵉ_func(i,j,k,l,G,λ) = (λ/(2G*(3λ + 2G)))*(δ(i,j)*δ(k,l)) + (1/2G)*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
