### PHYSICAL FUNCTIONS ###

# Second invariant of the deviatoric stress tensor :
get_τ(s::AbstractTensor,plas::DruckerPrager) = sqrt(0.5 * s ⊡ s)
get_τ(s::AbstractTensor,plas::VonMises) = sqrt(3/2 * s ⊡ s)

get_τ(s::AbstractTensor,r::Rheology) = get_τ(s::AbstractTensor,r.plasticity)
