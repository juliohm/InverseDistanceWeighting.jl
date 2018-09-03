using GeoStatsBase
using GeoStatsDevTools
using InverseDistanceWeighting
using Plots; gr()
using VisualRegressionTests
using Test, Pkg

# list of maintainers
maintainers = ["juliohm"]

# environment settings
istravis = "TRAVIS" ∈ keys(ENV)
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
datadir = joinpath(@__DIR__,"data")

if ismaintainer
  Pkg.add("Gtk")
  using Gtk
end

@testset "Basic problem" begin
  geodata = PointSetData(Dict(:variable => [1.,0.,1.]), [25. 50. 75.;  25. 75. 50.])
  domain = RegularGrid{Float64}(100,100)
  problem = EstimationProblem(geodata, domain, :variable)

  solver = InvDistWeight()

  solution = solve(problem, solver)

  if ismaintainer || istravis
    function plot_solution(fname)
      plot(solution, size=(1000,400))
      png(fname)
    end
    refimg = joinpath(datadir,"solution.png")

    @test test_images(VisualTest(plot_solution, refimg), popup=!istravis, tol=0.1) |> success
  end
end
