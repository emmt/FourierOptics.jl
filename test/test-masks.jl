using FourierOptics, Test, LinearAlgebra, TypeUtils, Unitful
using FourierOptics:
    GeometricObject, Point, Box, Rectangle, Circle, Polygon,
    center, diameter, radius, vertices

other_type(::Type{Int}) = Int16
other_type(::Type{Float64}) = Float32
other_type(::Type{Float32}) = Float64
other_type(::Type{Quantity{T,D,U}}) where {T,D,U} = Quantity{other_type(T),D,U}

function test(::Type{Point}, x::T, y::T) where T
    # Basic constructors.
    A = @inferred Point(x,y)
    @test A isa Point{T}
    @test A === @inferred Point{T}(x,y)

    # Properties.
    @test A.x === x
    @test A.y === y

    # Base methods.
    @test eltype(A) === T

    # Rebuild from parts.
    @test A === @inferred Point((A.x, A.y))
    @test A === @inferred Point(Tuple(A))
    @test A === @inferred Point(Tuple(A)...)
    @test A === @inferred Point(A...)
    @test Tuple(A) === (A...,)

    # Show.
    str = let io = IOBuffer()
        show(io, A)
        String(take!(io))
    end
    @test startswith(str, "Point(")

    # Conversion to other element type, ≈, and ==.
    @test A === @inferred convert(Point, A)
    @test A === @inferred convert(Point{eltype(A)}, A)
    U = other_type(T)
    B = @inferred Point{U}(A)
    @test B isa Point{U}
    @test eltype(B) === U
    @test B === @inferred Point{U}(A...)
    @test B ≈ A atol=abs(A)/10_000
    if bare_type(eltype(A)) <: Integer
        @test B == A
    end
    B = @inferred Point{U}(A)
    @test B === @inferred convert(Point{U}, A)
    @test A === @inferred convert_eltype(eltype(A), A)
    @test B === @inferred convert_eltype(U, A)

    # Mathematics.
    α = 2f0
    B = @inferred Point{U}(2*y, -5*x)
    @test +A === A
    @test -A === Point(-A.x, -A.y)
    @test A + B === Point(A.x + B.x, A.y + B.y)
    @test A - B === Point(A.x - B.x, A.y - B.y)
    C = @inferred A + B
    @test C.x === A.x + B.x && C.y === A.y + B.y
    C = @inferred α*A
    @test C.x === α*A.x && C.y === α*A.y
    @test C === @inferred A*α
    C = @inferred α\A
    @test C.x === α\A.x && C.y === α\A.y
    C = @inferred A/α
    @test C.x === A.x/α && C.y === A.y/α
    @test A*B === A.x*B.x + A.y*B.y
    @test FourierOptics.inner(A, B) === A.x*B.x + A.y*B.y
    @test FourierOptics.outer(A, B) === A.x*B.y - A.y*B.x

    # Polar coordinates.
    r = @inferred abs(A)
    @test r ≈ sqrt(x^2 + y^2)
    @test abs2(A) ≈ r^2
    @test norm(A) === r
    θ = uconvert(°, @inferred atan(A))
    @test θ ≈ atan(y, x)

    # Keyword-based constructors.
    @test A === @inferred Point(y=y, x=x)
    @test A === @inferred Point(; y, x)
    @test A ≈ @inferred Point(r=r, θ=θ)
    @test A ≈ @inferred Point(; r, θ)
    @test_throws Exception Point(; x)
    @test_throws Exception Point(; y)
    @test_throws Exception Point(; r)
    @test_throws Exception Point(; θ)
    @test_throws Exception Point(; x, y, r)
    @test_throws Exception Point(; r, θ, x)
end

function test(::Type{Box}, xmin::T, ymin::T, xmax::T, ymax::T) where T
    # Basic constructors.
    A = @inferred Box((xmin,ymin),(xmax,ymax))
    @test A isa Box{T}
    @test A === @inferred Box{T}((xmin,ymin),(xmax,ymax))

    # Properties.
    @test A.xmin === xmin
    @test A.ymin === ymin
    @test A.xmax === xmax
    @test A.ymax === ymax

    # Base methods.
    @test eltype(A) === T
    @test first(A) === Point(A.xmin, A.ymin)
    @test last(A) === Point(A.xmax, A.ymax)
    @test !isempty(A)

    # Coordinates are not sorted in a box.
    B = @inferred Box((xmax,ymin),(xmin,ymax))
    @test B.xmin === xmax
    @test B.ymin === ymin
    @test B.xmax === xmin
    @test B.ymax === ymax
    @test isempty(B)
    B = @inferred Box((xmin,ymax),(xmax,ymin))
    @test B.xmin === xmin
    @test B.ymin === ymax
    @test B.xmax === xmax
    @test B.ymax === ymin
    @test isempty(B)
    B = @inferred Box((xmax,ymax),(xmin,ymin))
    @test B.xmin === xmax
    @test B.ymin === ymax
    @test B.xmax === xmin
    @test B.ymax === ymin
    @test isempty(B)

    # Rebuild from parts.
    @test A === @inferred Box(first(A), last(A))
    @test A === @inferred Box(Tuple(A)...)
    @test A === @inferred Box(A...)
    @test Tuple(A) === (A...,)

    # Show.
    str = let io = IOBuffer()
        show(io, A)
        String(take!(io))
    end
    @test startswith(str, "Box(")

    # Conversion to other element type, ≈, and ==.
    @test A === @inferred convert(Box, A)
    @test A === @inferred convert(Box{eltype(A)}, A)
    U = other_type(T)
    B = @inferred Box{U}(A)
    @test B isa Box{U}
    @test eltype(B) === U
    @test B === @inferred Box{U}(A...)
    @test B ≈ A atol=max(map(abs, (xmin,ymin,xmax,ymax))...)/10_000
    if bare_type(eltype(A)) <: Integer
        @test B == A
    end
    B = @inferred Box{U}(A)
    @test B === @inferred convert(Box{U}, A)
    @test convert_eltype(eltype(A), A) === A
    @test B === @inferred convert_eltype(U, A)

    # Mathematics.
    α = 2f0
    @assert α > zero(α) # following tests assume this
    C = @inferred Point(-1.7*oneunit(T), 2.3*oneunit(T))
    @test +A === A
    @test -A === @inferred Box((-A.xmax, -A.ymax), (-A.xmin, -A.ymin))
    @test -A === @inferred Box(-last(A), -first(A))
    @test isempty(-A) == false
    @test (@inferred B + C) === (@inferred Box(first(B) + C, last(B) + C))
    @test (@inferred B + C) === (@inferred Box((B.xmin + C.x, B.ymin + C.y), (B.xmax + C.x, B.ymax + C.y)))
    @test (@inferred C + B) === (@inferred B + C)
    @test (@inferred B - C) === (@inferred Box(first(B) - C, last(B) - C))
    @test (@inferred B - C) === (@inferred Box((B.xmin - C.x, B.ymin - C.y), (B.xmax - C.x, B.ymax - C.y)))
    @test (@inferred C - B) === (@inferred -(B - C))
    @test (@inferred C - B) === (@inferred Box(C - last(B), C - first(B)))
    @test (@inferred C - B) === (@inferred Box((C.x - B.xmax, C.y - B.ymax), (C.x - B.xmin, C.y - B.ymin)))
    @test α*B === @inferred Box((α*B.xmin,α*B.ymin), (α*B.xmax,α*B.ymax))
    @test α*B === @inferred Box(α*first(B), α*last(B))
    @test B*3.1 === 3.1*B
    @test α\B === @inferred Box((α\B.xmin,α\B.ymin), (α\B.xmax,α\B.ymax))
    @test α\B === @inferred Box(α\first(B), α\last(B))
    @test B/α === @inferred Box((B.xmin/α,B.ymin/α), (B.xmax/α,B.ymax/α))
    @test B/α === @inferred Box(first(B)/α, last(B)/α)

    @test -α*B === @inferred Box((-α*B.xmax,-α*B.ymax), (-α*B.xmin,-α*B.ymin))
    @test -α*B === @inferred Box(-α*last(B), -α*first(B))
    @test B*-3.1 === -3.1*B
    @test -α\B === @inferred Box((-α\B.xmax,-α\B.ymax), (-α\B.xmin,-α\B.ymin))
    @test -α\B === @inferred Box(-α\last(B), -α\first(B))
    @test B/-α === @inferred Box((B.xmax/-α,B.ymax/-α), (B.xmin/-α,B.ymin/-α))
    @test B/-α === @inferred Box(last(B)/-α, first(B)/-α)

    # Keyword-based constructors.
    @test A === @inferred Box(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    @test A === @inferred Box(; xmin, ymin, xmax, ymax)
    @test_throws Exception Box(; xmin, ymin, xmax)
    @test_throws Exception Box(; ymax)
end

function test(::Type{Rectangle}, x0::T, y0::T, x1::T, y1::T) where T
    # Basic constructors.
    A = @inferred Rectangle((x0,y0),(x1,y1))
    @test A isa Rectangle{T}
    @test A === @inferred Rectangle{T}((x0,y0),(x1,y1))
    @test A === @inferred Rectangle((x1,y0),(x0,y1)) # order does not matter
    @test A === @inferred Rectangle((x1,y1),(x0,y0)) # order does not matter
    @test A === @inferred Rectangle((x0,y1),(x1,y0)) # order does not matter

    # Properties.
    @test A.x0 === min(x0,x1)
    @test A.y0 === min(y0,y1)
    @test A.x1 === max(x0,x1)
    @test A.y1 === max(y0,y1)

    # Base methods.
    @test eltype(A) === T
    @test first(A) === Point(A.x0, A.y0)
    @test last(A) === Point(A.x1, A.y1)

    # Rebuild from parts.
    @test A === @inferred Rectangle(first(A), last(A))
    @test A === @inferred Rectangle(last(A), first(A)) # order does not matter
    @test A === @inferred Rectangle(Tuple(A)...)
    @test A === @inferred Rectangle(A...)
    @test Tuple(A) === (A...,)

    # Show.
    str = let io = IOBuffer()
        show(io, A)
        String(take!(io))
    end
    @test startswith(str, "Rectangle(")

    # Conversion to other element type, ≈, and ==.
    @test A === @inferred convert(Rectangle, A)
    @test A === @inferred convert(Rectangle{eltype(A)}, A)
    U = other_type(T)
    B = @inferred Rectangle{U}(A)
    @test B isa Rectangle{U}
    @test eltype(B) === U
    @test B === @inferred Rectangle{U}(A...)
    @test B ≈ A atol=max(map(abs, (x0,y0,x1,y1))...)/10_000
    if bare_type(eltype(A)) <: Integer
        @test B == A
    end
    B = @inferred Rectangle{U}(A)
    @test B === @inferred convert(Rectangle{U}, A)
    @test convert_eltype(eltype(A), A) === A
    @test B === @inferred convert_eltype(U, A)

    # Mathematics.
    α = 2f0
    C = @inferred Point(-1.7*oneunit(T), 2.3*oneunit(T))
    @test +A === A
    @test -A === @inferred Rectangle((-A.x0, -A.y0), (-A.x1, -A.y1))
    @test -A === @inferred Rectangle(-first(A), -last(A))
    @test (@inferred B + C) === (@inferred Rectangle(first(B) + C, last(B) + C))
    @test (@inferred B + C) === (@inferred Rectangle((B.x0 + C.x, B.y0 + C.y), (B.x1 + C.x, B.y1 + C.y)))
    @test (@inferred C + B) === (@inferred B + C)
    @test (@inferred B - C) === (@inferred Rectangle(first(B) - C, last(B) - C))
    @test (@inferred B - C) === (@inferred Rectangle((B.x0 - C.x, B.y0 - C.y), (B.x1 - C.x, B.y1 - C.y)))
    @test (@inferred C - B) === (@inferred -(B - C))
    @test (@inferred C - B) === (@inferred Rectangle(C - first(B), C - last(B)))
    @test (@inferred C - B) === (@inferred Rectangle((C.x - B.x0, C.y - B.y0), (C.x - B.x1, C.y - B.y1)))
    @test α*B === @inferred Rectangle((α*B.x0,α*B.y0), (α*B.x1,α*B.y1))
    @test α*B === @inferred Rectangle(α*first(B), α*last(B))
    @test B*3.1 === 3.1*B
    @test α\B === @inferred Rectangle((α\B.x0,α\B.y0), (α\B.x1,α\B.y1))
    @test α\B === @inferred Rectangle(α\first(B), α\last(B))
    @test B/α === @inferred Rectangle((B.x0/α,B.y0/α), (B.x1/α,B.y1/α))
    @test B/α === @inferred Rectangle(first(B)/α, last(B)/α)

    # Keyword-based constructors.
    @test A === @inferred Rectangle(x0=x0, x1=x1, y0=y0, y1=y1)
    @test A === @inferred Rectangle(; x0, x1, y0, y1)
    @test_throws Exception Rectangle(; x0, y0, x1)
    @test_throws Exception Rectangle(; y1)
end

function test(::Type{Circle}, x::T, y::T, r::T) where T
    # Basic constructors.
    A = @inferred Circle((x,y), r)
    @test A isa Circle{T}
    @test_throws ArgumentError Circle((x,y), -r) # r < 0 is invalid
    @test A === @inferred Circle{T}((x,y), r)

    # Properties.
    @test A.x === x
    @test A.y === y
    @test A.r === r

    # Getters.
    @test center(A) === Point(A.x, A.y)
    @test radius(A) === A.r
    @test diameter(A) === 2radius(A)

    # Base methods.
    @test eltype(A) === T

    # Rebuild from parts.
    @test A === @inferred Circle(Tuple(A)...)
    @test A === @inferred Circle(A...)
    @test Tuple(A) === (A...,)

    # Show.
    str = let io = IOBuffer()
        show(io, A)
        String(take!(io))
    end
    @test startswith(str, "Circle(")

    # Conversion to other element type, ≈, and ==.
    @test A === @inferred convert(Circle, A)
    @test A === @inferred convert(Circle{eltype(A)}, A)
    U = other_type(T)
    B = @inferred Circle{U}(A)
    @test B isa Circle{U}
    @test eltype(B) === U
    @test B === @inferred Circle{U}(A...)
    @test B ≈ A atol=max(map(abs, (x,y,r))...)/10_000
    if bare_type(eltype(A)) <: Integer
        @test B == A
    end
    B = @inferred Circle{U}(A)
    @test B === @inferred convert(Circle{U}, A)
    @test convert_eltype(eltype(A), A) === A
    @test B === @inferred convert_eltype(U, A)

    # Mathematics.
    α = 2f0
    C = @inferred Point(-1.7*oneunit(T), 2.3*oneunit(T))
    @test +A === A
    @test -A === @inferred Circle((-A.x, -A.y), r)
    @test (@inferred B + C) === (@inferred Circle(center(B) + C, radius(B)))
    @test (@inferred B + C) === (@inferred Circle((B.x + C.x, B.y + C.y), B.r))
    @test (@inferred C + B) === (@inferred B + C)
    @test (@inferred B - C) === (@inferred Circle(center(B) - C, radius(B)))
    @test (@inferred B - C) === (@inferred Circle((B.x - C.x, B.y - C.y), B.r))
    @test (@inferred C - B) === (@inferred -(B - C))
    @test (@inferred C - B) === (@inferred Circle(C - center(B), radius(B)))
    @test (@inferred C - B) === (@inferred Circle((C.x - B.x, C.y - B.y), B.r))
    @test α*B === @inferred Circle((α*B.x, α*B.y), abs(α)*B.r)
    @test α*B === @inferred Circle(α*center(B), abs(α)*radius(B))
    @test -α*B === @inferred Circle((-α*B.x, -α*B.y), abs(-α)*B.r)
    @test -α*B === @inferred Circle(-α*center(B), abs(-α)*radius(B))
    @test B*3.1 === 3.1*B
    @test α\B === @inferred Circle((α\B.x, α\B.y), abs(α)\B.r)
    @test α\B === @inferred Circle(α\center(B), abs(α)\radius(B))
    @test -α\B === @inferred Circle((-α\B.x, -α\B.y), abs(-α)\B.r)
    @test -α\B === @inferred Circle(-α\center(B), abs(-α)\radius(B))
    @test B/α === @inferred Circle((B.x/α,B.y/α), B.r/abs(α))
    @test B/α === @inferred Circle(center(B)/α, radius(B)/abs(α))
    @test B/-α === @inferred Circle((B.x/-α, B.y/-α), B.r/abs(-α))
    @test B/-α === @inferred Circle(center(B)/-α, radius(B)/abs(-α))

    # Keyword-based constructors.
    @test A === @inferred Circle(; center=(x,y), radius=r)
    if bare_type(T) <: Integer
        @test A == @inferred Circle(; center=(x,y), diameter=2r)
    else
        @test A === @inferred Circle(; center=(x,y), diameter=2r)
    end
    @test_throws Exception Circle(; center=(x,y))
    @test_throws Exception Circle(; radius=r)
    @test_throws Exception Circle(; diameter=2r)
    @test_throws Exception Circle(; center=(x,y), radius=r, diameter=2r)
end

function test(::Type{GeometricObject})
    @testset "Geometric Objects" begin
        @testset "Points (T=$(typeof(first(args))))" for args in ((1,2),
                                                                  (-2.0,3.0),
                                                                  (2.1mm,-4.3mm))
            test(Point, args...)
        end
        @testset "Boxes (T=$(typeof(first(args)))))" for args in ((1,2,3,5),
                                                                  (-2.0,3.0,1.0,4.0),
                                                                  (-1.6mm,1.4mm,2.1mm,3.7mm))
            test(Box, args...)
        end
        @testset "Rectangles (T=$(typeof(first(args)))))" for args in ((1,2,3,5),
                                                                       (-2.0,3.0,1.0,4.0),
                                                                       (-1.6mm,1.4mm,2.1mm,3.7mm))
            test(Rectangle, args...)
        end
        @testset "Circles (T=$(typeof(first(args)))))" for args in ((1,2,3),
                                                                    (-1.0,2.0,3.0),
                                                                    (-1.6mm,1.4mm,2.7mm))
            test(Circle, args...)
        end
    end
end

nothing
