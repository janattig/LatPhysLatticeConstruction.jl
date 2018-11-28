################################################################################
#
#   module LatticePhysics_LatticeConstruction
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> NAIVE CONSTRUCTION OF LATTICES
#           - periodic pattern of unitcells
#           - semi-periodic pattern of unitcells
#           - open pattern of unitcells
#           - by bond distance to an origin site
#           - in a shape around an origin site
#
#   --> ADDING NNN TO LATTICE
#
#   --> LATTICE CONSTRUCTION ALGORITHM FROM SPIN SPIRALS PAPER
#
################################################################################

# module start
module LatPhysLatticeConstruction


# uses the Base package
using LatPhysBase







# FUNCTIONS TO DEFINE INDICES IN ARRAYS
include("index_functions.jl")

# CONSTRUCTION OF LATTICES FROM UNITCELLS
# BY PUTTING TOGETHER UNITCELLS
# (currently: only periodic patterns)
include("unitcell_patterns_periodic.jl")

# CONSTRUCTION OF LATTICES FROM UNITCELLS
# BY SPREADING OUTWARD OVER BONDS
include("bond_distance_based.jl")

# module end
end
