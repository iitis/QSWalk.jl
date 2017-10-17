facts("Global operator construction") do
  H = [1. 1.+im 3.; 1.-im 1. im; 3. -im 1.]
  L1 = sparse([1.+0im 2. 3.; 4. 5. 6.; 6. 7. -6.])
  L2 = sparse([0.+0im 0. 1.; 0. 0. 0.; 0. 0. 0.])
  resultnoomega = sparse([-52.0+0.0im  -29.0+1.0im    7.5+3.0im  -29.0-1.0im    4.0+0.0im     6.0+0.0im    7.5-3.0im     6.0+0.0im   10.0+0.0im;
                          -29.0+1.0im  -60.5+0.0im   10.0+0.0im    8.0+0.0im  -21.0-1.0im    12.0+0.0im   12.0+0.0im    19.5-3.0im   18.0+0.0im;
                           10.5+3.0im    9.0+0.0im  -73.5+0.0im   12.0+0.0im   14.0+0.0im   -43.0-1.0im   18.0+0.0im    21.0+0.0im  -13.5-3.0im;
                          -29.0-1.0im    8.0+0.0im   12.0+0.0im  -60.5+0.0im  -21.0+1.0im    19.5+3.0im   10.0+0.0im    12.0+0.0im   18.0+0.0im;
                           16.0+0.0im  -13.0-1.0im   24.0+0.0im  -13.0+1.0im  -53.0+0.0im    34.0+0.0im   24.0+0.0im    34.0+0.0im   36.0+0.0im;
                           24.0+0.0im   28.0+0.0im  -57.0-1.0im   34.5+3.0im   37.0+0.0im  -110.0+0.0im   36.0+0.0im    42.0+0.0im  -32.0+0.0im;
                           10.5-3.0im   12.0+0.0im   18.0+0.0im    9.0+0.0im   14.0+0.0im    21.0+0.0im  -73.5+0.0im   -43.0+1.0im  -13.5+3.0im;
                           24.0+0.0im   34.5-3.0im   36.0+0.0im   28.0+0.0im   37.0+0.0im    42.0+0.0im  -57.0+1.0im  -110.0+0.0im  -32.0+0.0im;
                           36.0+0.0im   42.0+0.0im  -31.5-3.0im   42.0+0.0im   49.0+0.0im   -40.0+0.0im  -31.5+3.0im   -40.0+0.0im  -46.0+0.0im])
  context("Standard usage") do
    #no locH case

    @fact global_operator(H, [L1, L2]) --> resultnoomega
    @fact global_operator(H, [L1, L2], 1/2) --> resultnoomega/2
  end

  context("Error tests") do
    @fact_throws MethodError global_operator(H, [L1, L2], 1im)
    @fact_throws ArgumentError global_operator(H, [L1, L2], -1)
    @fact_throws ArgumentError global_operator(H, [L1, L2], 3)
  end
end

facts("Classical lindlbad generators") do
  H = [1. 1.+im ; 1.-im im]
  result = [sparse([1.+0im 0 ; 0 0 ]), 
            sparse([0.+0im 1.+im ; 0 0 ]), 
            sparse([0.+0im 0 ; 1.-im 0 ]), 
            sparse([0.+0im 0 ; 0 im ])]
  resultwithepsilon = [sparse([0.+0im 1.+im ; 0 0 ]), 
            sparse([0.+0im 0 ; 1.-im 0 ])]
  context("Standard usage") do
    @fact classical_lindblad_operators(H) --> result
    @fact classical_lindblad_operators(sparse(H)) --> result
    @fact classical_lindblad_operators(H, epsilon = 1.1) --> resultwithepsilon
    @fact classical_lindblad_operators(sparse(H), epsilon = 1.1) --> resultwithepsilon
  end

  context("Type tests") do
    @fact typeof(classical_lindblad_operators(H)) --> Vector{SparseMatrixCSC{Complex128}}
    @fact typeof(classical_lindblad_operators(sparse(H))) --> Vector{SparseMatrixCSC{Complex128}}
  end

  context("Error tests") do
    @fact_throws ArgumentError classical_lindblad_operators(H, epsilon = -1.)
    @fact_throws TypeError classical_lindblad_operators(H, epsilon = -1im)
  end
end
