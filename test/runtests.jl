isdefined(:FourierOptics) || include("../src/FourierOptics.jl")

module TestFourierOptics

using FourierOptics
using FourierOptics.CoordinateTransforms
using Base.Test

Base.isapprox(x::NTuple{2,Real}, y::NTuple{2,Real}; kwds...) =
    (isapprox(x[1], y[1]; kwds...) && isapprox(x[2], y[2]; kwds...))

@testset "Coordinate Transforms" begin
    stp, x0, y0 = 1.2, -2.7, 3.1
    t = (0.7, -1.8)
    p = (0.9, 11.0)
    α = 0.17
    I = CoordinateTransform()
    A = CoordinateTransform(stp, x0, y0)
    B = inv(A)
    @test origin(A) ≈ (x0, y0)
    @test step(A) ≈ stp
    @test A ≈ (x0,y0) + stp*I
    @test A∘B ≈ I  atol=4*eps(Float64)
    @test B∘A ≈ I  atol=4*eps(Float64)
    @test A/A ≈ I  atol=4*eps(Float64)
    @test A\A ≈ I  atol=4*eps(Float64)
    B = t + A
    @test B(p) ≈ t .+ A(p)
    @test (t - A)(p) ≈ t .- A(p)
    @test (A + t)(p) ≈ A(p .+ t)
    @test (A - t)(p) ≈ A(p .- t)
    @test (α*A)(p) ≈ α.*A(p)
    @test (A*α)(p) ≈ A(α.*p)
end

end
