facts("basic arrays generators") do
  context("ket") do
    #standard tests
    @fact ket(1, 2) --> [1;0]
    #error tests
    @fact_throws AssertionError ket(4, 2)
    @fact_throws ArgumentError ket(-4, -2)
  end

  context("bra") do
    #standard tests
    @fact bra(1, 2) --> [1 0 ]
    #error tests
    @fact_throws AssertionError bra(4, 2)
    @fact_throws ArgumentError bra(-4, -2)
  end

  context("ketbra") do
    #standard tests
    @fact ketbra(1, 2, 3) -->  [0 1 0; 0 0 0; 0 0 0]
    #error tests
    @fact_throws AssertionError ketbra(3, 2, 2)
    @fact_throws AssertionError ketbra(2, 3, 2)
    @fact_throws ArgumentError ketbra(-4, -2, -1)
  end
end

facts("Array reshuffles") do
  context("proj") do
    #standard tests#
    result = [0.0+0.0im 0.0+0.0im 0.0+0.0im;
              0.0+0.0im 1.0+0.0im 0.0+0.0im;
              0.0+0.0im 0.0+0.0im 0.0+0.0im]
    @fact proj(2, 3) -->  result
    @fact proj(1/sqrt(2) * (ket(1, 3)+ket(3, 3))) -->
                          roughly([0.5+0.0im 0.0+0.0im 0.5+0.0im;
                          0.0+0.0im 0.0+0.0im 0.0+0.0im;
                          0.5+0.0im 0.0+0.0im 0.5+0.0im])
    #error tests
    @fact_throws ArgumentError proj(2, -1)
    @fact_throws AssertionError proj(3, 2)
  end

  context("res and unres") do
    #standard tests
    M = Matrix{Float64}(reshape(1:9, (3, 3))')
    v = Vector{Float64}(collect(1:9))
    A = Complex{Float64}[0.354177+0.0im 0.0891553-0.0251879im 0.0702961+0.0516828im 0.0708664+0.0767941im; 0.0891553+0.0251879im 0.336055+0.0im 0.0420202-0.0109173im 0.0683605-0.00692846im; 0.0702961-0.0516828im 0.0420202+0.0109173im 0.212401+0.0im 0.0939615+0.0553555im; 0.0708664-0.0767941im 0.0683605+0.00692846im 0.0939615-0.0553555im 0.0973671+0.0im]
    @fact res(M) --> v
    @fact unres(v) --> M
    @fact unres(res(M)) --> M
    @fact res(unres(v)) --> v
    @fact unres(res(A)) --> A
    #error tests
    @fact_throws ArgumentError unres(collect(1:8)*1.)
  end
end
