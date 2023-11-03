using FourierOptics
using Test
using Unitful

@testset "FourierOptics.jl" begin
    @testset "Utilities" begin
        @test one(Float32)*u"m" isa FourierOptics.Meters{Float32}
        @test FourierOptics.standard_length(1m) === 1.0
        @test FourierOptics.standard_length(Float32, 1m) === 1.0f0
        @test FourierOptics.standard_length(200cm) === 2.0
        @test FourierOptics.standard_length(Float32, 3000mm) === 3.0f0

        # Infinity
        let infinity = FourierOptics.infinity
            @test infinity(0.0) === Inf
            @test infinity(Float32) === Inf32
            @test infinity(1.0m) === Inf*m
            @test_throws Exception infinity(42)
            @test_throws Exception infinity(42μm)
        end
    end
    @testset "Fields" begin
        n = 16
        λ = 3μm
        Δx = 0.1mm
        F = @inferred Field(grid_size=n, wavelength=λ, grid_step=Δx)
        @test eltype(F) == Complex{Float64}
        @test ndims(F) == 2
        @test size(F) == (n, n)
        @test axes(F) == (1:n, 1:n)
        @test F[1,1] === 1.0 + 0.0im
        F[end,end] = 0
        @test F[n,n] === 0.0 + 0.0im

        # Test copy constructor.
        for j in eachindex(F) # make sure that the complex amplitude is not just ones
            F[j] = complex(j, -j)
        end
        Fcpy = @inferred copy(F)
        @test typeof(Fcpy) === typeof(F)
        @test FourierOptics.get_amplitude(Fcpy) !== FourierOptics.get_amplitude(F)
        for key in fieldnames(typeof(F))
            @test getfield(F, key) === getfield(Fcpy, key) || getfield(F, key) == getfield(Fcpy, key)
        end
    end
end
