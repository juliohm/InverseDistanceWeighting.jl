## Copyright (c) 2017, Júlio Hoffimann Mendes <juliohm@stanford.edu>
##
## Permission to use, copy, modify, and/or distribute this software for any
## purpose with or without fee is hereby granted, provided that the above
## copyright notice and this permission notice appear in all copies.
##
## THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
## WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
## ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
## WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
## ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
## OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

__precompile__()

module InverseDistanceWeighting

importall GeoStatsBase
using GeoStatsDevTools

using Reexport
using NearestNeighbors
using StaticArrays
@reexport using Distances

export InvDistWeight

"""
    InvDistWeight(var₁=>param₁, var₂=>param₂, ...)

Inverse distance weighting estimation solver.

## Parameters

* `neighbors` - Number of neighbors (default to all data locations)
* `distance`  - A distance defined in Distances.jl (default to Euclidean()
"""
@estimsolver InvDistWeight begin
  @param neighbors = nothing
  @param distance = Euclidean()
end

function solve(problem::EstimationProblem, solver::InvDistWeight)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # result for each variable
  μs = []; σs = []

  for (var,V) in variables(problem)
    # get user parameters
    if var ∈ keys(solver.params)
      varparams = solver.params[var]
    else
      varparams = InvDistWeightParam()
    end

    # get valid data for variable
    X, z = valid(pdata, var)

    # number of data points for variable
    ndata = length(z)

    # allocate memory
    varμ = Vector{V}(npoints(pdomain))
    varσ = Vector{V}(npoints(pdomain))

    if ndata > 0
      # fit search tree
      kdtree = KDTree(X, varparams.distance)

      # keep track of estimated locations
      estimated = falses(npoints(pdomain))

      # consider data locations as already estimated
      for (loc, datloc) in datamap(problem, var)
        estimated[loc] = true
        varμ[loc] = value(pdata, datloc, var)
        varσ[loc] = zero(V)
      end

      # determine number of nearest neighbors to use
      k = varparams.neighbors == nothing ? ndata : varparams.neighbors

      @assert k ≤ ndata "number of neighbors must be smaller or equal to number of data points"

      # pre-allocate memory for coordinates
      coords = MVector{ndims(pdomain),coordtype(pdomain)}()

      # estimation loop
      for location in SimplePath(pdomain)
        if !estimated[location]
          coordinates!(coords, pdomain, location)

          idxs, dists = knn(kdtree, coords, k)

          weights = one(V) ./ dists
          weights /= sum(weights)

          varμ[location] = weights ⋅ z[idxs]
          varσ[location] = minimum(dists)
        end
      end
    end

    push!(μs, var => varμ)
    push!(σs, var => varσ)
  end

  EstimationSolution(pdomain, Dict(μs), Dict(σs))
end

end
