# Helper function

# N=0
function concanateBonds(
            b1 :: B,
            b2 :: B,
            labelfunction :: Function = (bo1,bo2) -> getDefaultLabel(L)
        ) :: B where {L,B<:AbstractBond{L,0}}

    # assert that the bonds can indeed by concanated
    @assert to(b1) == from(b2) "Bonds cannot be concanated as "*string(from(b1))*"-"*string(to(b1))*" (-) "*string(from(b2))*"-"*string(to(b2))*" not possible"

    # create a new bond
    return newBond(
        B,
        to(b1) :: Int64,
        to(b2) :: Int64,
        labelfunction(b1,b2) :: L,
        NTuple{0,Int64}()
    )
end
# N=1
function concanateBonds(
            b1 :: B,
            b2 :: B,
            labelfunction :: Function = (bo1,bo2) -> getDefaultLabel(L)
        ) :: B where {L,B<:AbstractBond{L,1}}

    # assert that the bonds can indeed by concanated
    @assert to(b1) == from(b2) "Bonds cannot be concanated as "*string(from(b1))*"-"*string(to(b1))*" (-) "*string(from(b2))*"-"*string(to(b2))*" not possible"

    # create a new bond
    return newBond(
        B,
        to(b1) :: Int64,
        to(b2) :: Int64,
        labelfunction(b1,b2) :: L,
        (
            -wrap(b1)[1] + wrap(b2)[1],
        )
    )
end
# N=2
function concanateBonds(
            b1 :: B,
            b2 :: B,
            labelfunction :: Function = (bo1,bo2) -> getDefaultLabel(L)
        ) :: B where {L,B<:AbstractBond{L,2}}

    # assert that the bonds can indeed by concanated
    @assert to(b1) == from(b2) "Bonds cannot be concanated as "*string(from(b1))*"-"*string(to(b1))*" (-) "*string(from(b2))*"-"*string(to(b2))*" not possible"

    # create a new bond
    return newBond(
        B,
        to(b1) :: Int64,
        to(b2) :: Int64,
        labelfunction(b1,b2) :: L,
        (
            -wrap(b1)[1] + wrap(b2)[1],
            -wrap(b1)[2] + wrap(b2)[2]
        )
    )
end
# N=3
function concanateBonds(
            b1 :: B,
            b2 :: B,
            labelfunction :: Function = (bo1,bo2) -> getDefaultLabel(L)
        ) :: B where {L,B<:AbstractBond{L,3}}

    # assert that the bonds can indeed by concanated
    @assert to(b1) == from(b2) "Bonds cannot be concanated as "*string(from(b1))*"-"*string(to(b1))*" (-) "*string(from(b2))*"-"*string(to(b2))*" not possible"

    # create a new bond
    return newBond(
        B,
        to(b1) :: Int64,
        to(b2) :: Int64,
        labelfunction(b1,b2) :: L,
        (
            -wrap(b1)[1] + wrap(b2)[1],
            -wrap(b1)[2] + wrap(b2)[2],
            -wrap(b1)[3] + wrap(b2)[3]
        )
    )
end


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
