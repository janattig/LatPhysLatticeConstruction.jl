# create NNN for an arbitrary lattice by bond distance
function getBondsNNNByBondDistance(
            unitcell :: U,
            labelfunction :: Function = (b1,b2) -> getDefaultLabel(L,"NNN")
        ) :: Vector{B} where {N,L,S,B<:AbstractBond{L,N},U<:AbstractUnitcell{S,B}}

    # obtain the organized bond lists
    bonds_organized_from = organizedBondsFrom(unitcell)
    bonds_organized_to   = organizedBondsTo(unitcell)

    # create all possible bonds by concanation of nearest neighbors
    bonds_NNN = B[
        concanateBonds(b1, b2, labelfunction)
        for i in 1:numSites(unitcell)
        for b1 in bonds_organized_to[i]
        for b2 in bonds_organized_from[i]
    ]
    # filter out the ones that connect to the same site
    filter!(b -> !(!isPeriodic(b) && from(b)==to(b)), bonds_NNN)

    # return the remaining bonds
    return bonds_NNN
end

# export the function
export getBondsNNNByBondDistance
