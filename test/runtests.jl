using Clouds
using Test

@testset "Clouds.jl" begin
    @test sum(vec(box2d())) == 16
    # Write your tests here.
end
