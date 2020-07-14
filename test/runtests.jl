using GeoStatsBase
using InverseDistanceWeighting
using Plots; gr(size=(1000,400))
using VisualRegressionTests
using Distances
using Test, Pkg

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")
visualtests = !istravis || (istravis && islinux)
if !istravis
  Pkg.add("Gtk")
  using Gtk
end

@testset "Basic problem" begin
  geodata = georef((variable=[1.,0.,1.],), [25. 50. 75.;  25. 75. 50.])
  domain  = RegularGrid(100,100)
  problem = EstimationProblem(geodata, domain, :variable)

  solver = InvDistWeight(:variable => (neighbors=3,))

  solution = solve(problem, solver)

  if visualtests
    @plottest contourf(solution) joinpath(datadir,"solution.png") !istravis
  end
end

@testset "Haversine" begin
  geodata = georef((variable=[4.0,-1.0,3.0],), [50.0 100.0 200.0; -30.0 30.0 10.0])
  domain  = RegularGrid((1.0, -89.0), (359.0, 89.0), dims=(200, 100))
  problem = EstimationProblem(geodata, domain, :variable)

  solver = InvDistWeight(:variable => (distance=Haversine(1.0),))

  solution = solve(problem, solver)

  if visualtests
    @plottest contourf(solution) joinpath(datadir,"haversine.png") !istravis
  end
end
