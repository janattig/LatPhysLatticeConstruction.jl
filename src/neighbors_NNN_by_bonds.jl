# Helper function

# N=0
function concanateBonds(
            b1 :: B,
            b2 :: B,
            labelfunction :: Function
        ) where {L,B<:AbstractBond{L,0}}

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
            labelfunction :: Function
        ) where {L,B<:AbstractBond{L,1}}

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
            labelfunction :: Function
        ) where {L,B<:AbstractBond{L,2}}

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
            labelfunction :: Function
        ) where {L,B<:AbstractBond{L,3}}

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

    # obtain the organized bond list
    bonds_organized = organizedBonds(unitcell)

    # return the new bonds
    return B[
        concanateBonds(b1, b2, labelfunction)
        for i in 1:numSites(unitcell)
        for b1 in bonds_organized[i]
        for b2 in bonds_organized[i] if b1!=b2
    ]
end

# export the function
export getBondsNNNByBondDistance
