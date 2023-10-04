# ABS-Majoranas

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ABS-Majoranas

This code is used for the calculations and plots in the paper **A minimal quantum-dot-based Kitaev chain with only local superconducting proximity effect** by William Samuelson, Viktor Svensson, and Martin Leijnse.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
s
This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ABS-Majoranas"
```
which auto-activate the project and enable local path handling from DrWatson.

Little effort has been made to make this user friendly. If you have any questions, contact Viktor Svensson.