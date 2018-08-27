@testset "Evolution generator construction" begin
  H = [1. 1. + im 3.; 1. - im 1. im; 3. -im 1.]
  L1 = sparse([1. + 0im 2. 3.; 4. 5. 6.; 6. 7. -6.])
  L2 = sparse([0. + 0im 0. 1.; 0. 0. 0.; 0. 0. 0.])
  resultnoomega = sparse([-52.0+0.0im  -29.0+1.0im    7.5+3.0im  -29.0-1.0im    4.0+0.0im     6.0+0.0im    7.5-3.0im     6.0+0.0im   10.0+0.0im;
                          -29.0+1.0im  -60.5+0.0im   10.0+0.0im    8.0+0.0im  -21.0-1.0im    12.0+0.0im   12.0+0.0im    19.5-3.0im   18.0+0.0im;
                           10.5+3.0im    9.0+0.0im  -73.5+0.0im   12.0+0.0im   14.0+0.0im   -43.0-1.0im   18.0+0.0im    21.0+0.0im  -13.5-3.0im;
                          -29.0-1.0im    8.0+0.0im   12.0+0.0im  -60.5+0.0im  -21.0+1.0im    19.5+3.0im   10.0+0.0im    12.0+0.0im   18.0+0.0im;
                           16.0+0.0im  -13.0-1.0im   24.0+0.0im  -13.0+1.0im  -53.0+0.0im    34.0+0.0im   24.0+0.0im    34.0+0.0im   36.0+0.0im;
                           24.0+0.0im   28.0+0.0im  -57.0-1.0im   34.5+3.0im   37.0+0.0im  -110.0+0.0im   36.0+0.0im    42.0+0.0im  -32.0+0.0im;
                           10.5-3.0im   12.0+0.0im   18.0+0.0im    9.0+0.0im   14.0+0.0im    21.0+0.0im  -73.5+0.0im   -43.0+1.0im  -13.5+3.0im;
                           24.0+0.0im   34.5-3.0im   36.0+0.0im   28.0+0.0im   37.0+0.0im    42.0+0.0im  -57.0+1.0im  -110.0+0.0im  -32.0+0.0im;
                           36.0+0.0im   42.0+0.0im  -31.5-3.0im   42.0+0.0im   49.0+0.0im   -40.0+0.0im  -31.5+3.0im   -40.0+0.0im  -46.0+0.0im])
  @testset "Standard usage" begin
    #no locH case

    @test evolve_generator(H, [L1, L2]) == resultnoomega
    @test evolve_generator(H, [L1, L2], 1/2) == resultnoomega/2
  end
  @testset "Nonmoralized usage" begin
    globalH = nm_glob_ham(H)
    Lnonmoral1 = nm_lind(L1)
    locH1 = nm_loc_ham(Lnonmoral1[2])
    @test evolve_generator(globalH, SparseMatrixCSC{ComplexF64}[], locH1, 1/3) ≈ evolve_generator((1-1/3)*globalH+1/3*locH1,SparseMatrixCSC{ComplexF64}[])
    @test evolve_generator(globalH, [Lnonmoral1[1]], locH1, 1/3) ≈ evolve_generator(globalH+1/2*(locH1),[Lnonmoral1[1]],1/3)

  end

  @testset "Error tests" begin
    @test_throws MethodError evolve_generator(H, [L1, L2], 1im)
    @test_throws ArgumentError evolve_generator(H, [L1, L2], -1)
    @test_throws ArgumentError evolve_generator(H, [L1, L2], 3)
  end
end

@testset "Classical lindblad generators" begin
  H = [1. 1.0+im ; 1.0-im im]
 result = [sparse([1.0+0im 0 ; 0 0 ]),
            sparse([0.0+0im 1.0+im ; 0 0 ]),
            sparse([0.0+0im 0 ; 1.0-im 0 ]),
            sparse([0.0+0im 0 ; 0 im ])]
  resultwithepsilon = [sparse([0.0+0im 1.0+im ; 0 0 ]),
           sparse([0.0+0im 0 ; 1.0-im 0 ])]
  @testset "Standard usage" begin
    @test local_lind(H) == result
    @test local_lind(sparse(H)) == result
    @test local_lind(H, epsilon = 1.1) == resultwithepsilon
    @test local_lind(sparse(H), epsilon = 1.1) == resultwithepsilon
  end

  @testset "Type tests" begin
    @test typeof(local_lind(H)) == Vector{SparseMatrixCSC{ComplexF64}}
    @test typeof(local_lind(sparse(H))) == Vector{SparseMatrixCSC{ComplexF64}}
  end

  @testset "Error tests" begin
    @test_throws ArgumentError local_lind(H, epsilon = -1.)
    @test_throws TypeError local_lind(H, epsilon = -1im)
  end
end
