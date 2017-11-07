using GeoStats
using InverseDistanceWeighting
using Plots; gr()
using Base.Test
using VisualRegressionTests

# setup GR backend for Travis CI
ENV["GKSwstype"] = "100"

# list of maintainers
maintainers = ["juliohm"]

# environment settings
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
istravislinux = "TRAVIS" ∈ keys(ENV) && ENV["TRAVIS_OS_NAME"] == "linux"
datadir = joinpath(@__DIR__,"data")

@testset "Basic problem" begin
  geodata = GeoDataFrame(DataFrames.DataFrame(x=[25.,50.,75.], y=[25.,75.,50.], variable=[1.,0.,1.]), [:x,:y])
  domain = RegularGrid{Float64}(100,100)
  problem = EstimationProblem(geodata, domain, :variable)

  solver = InvDistWeight()

  solution = solve(problem, solver)

  if ismaintainer || istravislinux
    function plot_solution(fname)
      plot(solution, size=(1000,400))
      png(fname)
    end
    refimg = joinpath(datadir,"solution.png")

    @test test_images(VisualTest(plot_solution, refimg), popup=!istravislinux) |> success
  end
end
