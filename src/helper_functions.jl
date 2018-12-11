################################################################################
#
#   HELPER FUNCTIONS FOR CONSTRUCTION OF LATTICES / NEIGHBORS
#
#   FILE CONTAINS
#   - bond concanation
#
################################################################################





################################################################################
#
#   CONCANATION OF BONDS
#
################################################################################

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
