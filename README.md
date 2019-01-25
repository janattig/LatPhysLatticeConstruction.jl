# LatPhysLatticeConstruction.jl

Lattice construction module of [`LatticePhysics.jl`](https://github.com/janattig/LatticePhysics.jl).



## Contents

Provides functions to generate (construct) lattices out of unitcells in various ways. Currently implemented:
1.  Periodic patterns of unitcells
2.  Lattices based on bond-distance
3.  Lattices grown in some shape (e.g. sphere or box)


## Installation

You can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to first install the unregistered
dependencies manually with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysBase.jl"
```
to finally install the main package with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeConstruction.jl"
```
