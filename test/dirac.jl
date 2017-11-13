@testset "Basic arrays generators" begin
  @testset "ket" begin
    #standard tests
    @test ket(1, 2) == [1;0]
    #error tests
    @test_throws AssertionError ket(4, 2)
    @test_throws ArgumentError ket(-4, -2)
  end

  @testset "bra" begin
    #standard tests
    @test bra(1, 2) == [1 0 ]
    #error tests
    @test_throws AssertionError bra(4, 2)
    @test_throws ArgumentError bra(-4, -2)
  end

  @testset "ketbra" begin
    #standard tests
    @test ketbra(1, 2, 3) ==  [0 1 0; 0 0 0; 0 0 0]
    #error tests
    @test_throws AssertionError ketbra(3, 2, 2)
    @test_throws AssertionError ketbra(2, 3, 2)
    @test_throws ArgumentError ketbra(-4, -2, -1)
  end

  @testset "proj" begin
    #standard tests#
    result = [0.0+0.0im 0.0+0.0im 0.0+0.0im;
              0.0+0.0im 1.0+0.0im 0.0+0.0im;
              0.0+0.0im 0.0+0.0im 0.0+0.0im]
    @test proj(2, 3) ==  result
    @test proj(1/sqrt(2) * (ket(1, 3)+ket(3, 3))) ≈
                         [0.5+0.0im 0.0+0.0im 0.5+0.0im;
                          0.0+0.0im 0.0+0.0im 0.0+0.0im;
                          0.5+0.0im 0.0+0.0im 0.5+0.0im]
    #error tests
    @test_throws ArgumentError proj(2, -1)
    @test_throws AssertionError proj(3, 2)
  end
end

@testset "Array reshuffles" begin
  @testset "res and unres" begin
    #standard tests
    M = Matrix{Float64}(reshape(1:9, (3, 3))')
    v = Vector{Float64}(collect(1:9))
    A = Complex{Float64}[0.354177+0.0im 0.0891553-0.0251879im 0.0702961+0.0516828im 0.0708664+0.0767941im; 0.0891553+0.0251879im 0.336055+0.0im 0.0420202-0.0109173im 0.0683605-0.00692846im; 0.0702961-0.0516828im 0.0420202+0.0109173im 0.212401+0.0im 0.0939615+0.0553555im; 0.0708664-0.0767941im 0.0683605+0.00692846im 0.0939615-0.0553555im 0.0973671+0.0im]
    @test res(M) == v
    @test unres(v) == M
    @test unres(res(M)) == M
    @test res(unres(v)) == v
    @test unres(res(A)) == A
    #error tests
    @test_throws ArgumentError unres(collect(1:8)*1.)
  end
end
