isdefined(:FourierOptics) || include("../src/FourierOptics.jl")

module TestFourierOptics

using Base.Test

using FourierOptics
using FourierOptics.CoordinateTransforms
using FourierOptics.Regions

Base.isapprox(x::NTuple{2,Real}, y::NTuple{2,Real}; kwds...) =
    (isapprox(x[1], y[1]; kwds...) && isapprox(x[2], y[2]; kwds...))

Base.isapprox(x::NTuple{4,Real}, y::NTuple{4,Real}; kwds...) =
    (isapprox(x[1], y[1]; kwds...) && isapprox(x[2], y[2]; kwds...) &&
     isapprox(x[3], y[3]; kwds...) && isapprox(x[4], y[4]; kwds...))

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

@testset "Regions" begin
    stp, xc, yc = 1.2, -2.7, 3.1
    n1, n2 = 27, 52
    c = 0.7, -1.1
    grdxmin = xc - (n1 - 1)*stp/2
    grdxmax = xc + (n1 - 1)*stp/2
    grdymin = yc - (n2 - 1)*stp/2
    grdymax = yc + (n2 - 1)*stp/2
    boxxmin = xc - n1*stp/2
    boxxmax = xc + n1*stp/2
    boxymin = yc - n2*stp/2
    boxymax = yc + n2*stp/2
    X = linspace(grdxmin, grdxmax, n1)
    Y = linspace(grdymin, grdymax, n2)
    R1 = Region(X, Y)
    @test center(R1) ≈ (xc, yc)
    @test step(R1) ≈ stp
    @test size(R1) ≈ (n1, n2)
    @test extrema(R1) ≈ (grdxmin, grdxmax, grdymin, grdymax)
    @test boundingbox(R1) ≈ (boxxmin, boxxmax, boxymin, boxymax)
    R2 = recenter(R1)
    @test center(R2) ≈ (0,0)
    R2 = recenter(R1, c)
    @test center(R2) ≈ c
    @test center(R2 - c) ≈ (0,0)
end

end
