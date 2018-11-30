using LatPhysBase
using LatPhysLatticeConstruction
using Test

@testset "LatPhysLatticeConstruction.jl" begin

    #------------------
    # create unitcells
    #------------------

    # in 2d
    unitcell_2d = newUnitcell(
        # Type of the unitcell
        Unitcell{Site{Int64,2}, Bond{Symbol,2}},
        # Bravais lattice vectors
        Vector{Float64}[
            Float64[sqrt(3.0)/2, -0.5],
            Float64[sqrt(3.0)/2, +0.5]
        ],
        # Sites
        Site{Int64,2}[
            newSite(Site{Int64,2}, Float64[0.0, 0.0],         getDefaultLabel(Int64)),
            newSite(Site{Int64,2}, Float64[1/sqrt(3.0), 0.0], getDefaultLabel(Int64))
        ],
        # Bonds
        Bond{Symbol,2}[
            newBond(Bond{Symbol,2}, 1, 2, getDefaultLabel(Symbol), (0, 0)),
            newBond(Bond{Symbol,2}, 1, 2, getDefaultLabel(Symbol), (-1, 0)),
            newBond(Bond{Symbol,2}, 1, 2, getDefaultLabel(Symbol), (0, -1)),
            newBond(Bond{Symbol,2}, 2, 1, getDefaultLabel(Symbol), (0, 0)),
            newBond(Bond{Symbol,2}, 2, 1, getDefaultLabel(Symbol), (1, 0)),
            newBond(Bond{Symbol,2}, 2, 1, getDefaultLabel(Symbol), (0, 1))
        ]
    )

    # in 3d
    unitcell_3d = newUnitcell(
        # Type of the unitcell
        Unitcell{Site{Int64,3}, Bond{Symbol,3}},
        # lattice vectors
        Vector{Float64}[
            Float64[1, 0, 0],
            Float64[0, 1, 0],
            Float64[0, 0, 1]
        ],
        # sites
        Site{Int64,3}[
            newSite(Site{Int64,3}, Float64[0,0,0], getDefaultLabel(Int64))
        ],
        # bonds
        Bond{Symbol,3}[
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (+1,0,0)),
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (-1,0,0)),
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (0,+1,0)),
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (0,-1,0)),
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (0,0,+1)),
            newBond(Bond{Symbol,3}, 1,1, getDefaultLabel(Symbol), (0,0,-1))
        ]
    )




    #------------------
    # periodic lattice
    #------------------
    @testset "Periodic lattices" begin

        # 2d given all extents
        @test getLatticePeriodic(unitcell_2d, (2,2))
        # 2d given only 1 extent
        @test getLatticePeriodic(unitcell_2d, 3)

        # 3d given all extents
        @test getLatticePeriodic(unitcell_3d, (2,2,2))
        # 3d given only 1 extent
        @test getLatticePeriodic(unitcell_3d, 3)

    end



    #------------------
    # by bond distance
    #------------------
    @testset "Lattices by bond distance" begin

        # 2d given only the bond distance
        @test getLatticeByBondDistance(unitcell_2d, 3)
        # 2d given bond distance and origin
        @test getLatticeByBondDistance(unitcell_2d, 3, 2)


        # 3d given only the bond distance
        @test getLatticeByBondDistance(unitcell_3d, 3)
        # 3d given bond distance and origin
        @test getLatticeByBondDistance(unitcell_3d, 3, 1)

    end

end
