using FactCheck
using QSWalk

facts("Basic linear util functions") do
  context("ket") do
    #standard tests
    @fact ket(1,2) --> [1.0+0.0im;0.0+0.0im]
    #type tests
    @fact typeof(ket(Complex128,1,2)) --> typeof(ket(1,2))
    @fact typeof(ket(Float64,1,2)) --> Vector{Float64}
    #error tests
    @fact_throws AssertionError ket(4,2)
    @fact_throws ArgumentError ket(-4,-2)
  end
  context("bra") do
    #standard tests
    @fact bra(1,2) --> [1.0+0.0im 0.0+0.0im]
    #type tests
    @fact typeof(bra(Complex128,1,2)) --> typeof(bra(1,2))
    @fact typeof(bra(Float64,1,2)) --> Matrix{Float64}
    #error tests
    @fact_throws AssertionError, bra(4,2)
    @fact_throws ArgumentError, bra(-4,-2)
  end
  context("ketbra") do
    #standard tests
    @fact ketbra(1,2,3) -->  [0.0+0.0im 1.0+0.0im 0.0+0.0im;
                              0.0+0.0im 0.0+0.0im 0.0+0.0im;
                              0.0+0.0im 0.0+0.0im 0.0+0.0im]
    #type tests
    @fact typeof(ketbra(Complex128,1,2,3)) --> typeof(ketbra(1,2,3))
    @fact typeof(ketbra(Float64,1,2,3)) --> Matrix{Float64}
    #error tests
    @fact_throws AssertionError, ketbra(3,2,2)
    @fact_throws AssertionError, ketbra(2,3,2)
    @fact_throws ArgumentError, ketbra(-4,-2,-1)
  end
  context("proj") do
    #standard tests#
    result = [0.0+0.0im 0.0+0.0im 0.0+0.0im;
                          0.0+0.0im 1.0+0.0im 0.0+0.0im;
                          0.0+0.0im 0.0+0.0im 0.0+0.0im]
    @fact proj(2,3) -->  result
    @fact proj(1/sqrt(2) * (ket(1,3)+ket(3,3))) -->
                          roughly([0.5+0.0im 0.0+0.0im 0.5+0.0im;
                          0.0+0.0im 0.0+0.0im 0.0+0.0im;
                          0.5+0.0im 0.0+0.0im 0.5+0.0im])
    #type tests
    @fact typeof(proj(Complex128,2,3)) -->  typeof(proj(2,3))
    @fact typeof(proj(Float64,2,3)) --> Matrix{Float64}
    #error tests
    @fact_throws ArgumentError, proj(2,-1)
    @fact_throws AssertionError, proj(3,2)
  end
  context("reshuffle and reshuffle") do
    M = Matrix{Float64}(reshape(1:9, (3,3))')
    v = Vector{Float64}(collect(1:9))
    @fact reshuffle(M) --> v
    @fact unreshuffle(v) --> M
    @fact unreshuffle(reshuffle(M)) --> M
    @fact reshuffle(unreshuffle(v)) --> v
    @fact_throws ArgumentError unreshuffle(collect(1:8)*1.)
  end
end

facts("Global operator preparation") do
  context("globaloperator") do
    #no locH case
    H = [1. 1.+im 3.; 1.-im 1. im; 3. -im 1.]
    L1 = sparse([1.+0im 2. 3.; 4. 5. 6.; 6. 7. -6.])
    L2 = sparse([ 0.+0im 0. 1. ; 0. 0. 0. ; 0. 0. 0.])
    resultnoomega = sparse([-52.0+0.0im  -29.0+1.0im    7.5+3.0im  -29.0-1.0im    4.0+0.0im     6.0+0.0im    7.5-3.0im     6.0+0.0im   10.0+0.0im;
                            -29.0+1.0im  -60.5+0.0im   10.0+0.0im    8.0+0.0im  -21.0-1.0im    12.0+0.0im   12.0+0.0im    19.5-3.0im   18.0+0.0im;
                             10.5+3.0im    9.0+0.0im  -73.5+0.0im   12.0+0.0im   14.0+0.0im   -43.0-1.0im   18.0+0.0im    21.0+0.0im  -13.5-3.0im;
                            -29.0-1.0im    8.0+0.0im   12.0+0.0im  -60.5+0.0im  -21.0+1.0im    19.5+3.0im   10.0+0.0im    12.0+0.0im   18.0+0.0im;
                             16.0+0.0im  -13.0-1.0im   24.0+0.0im  -13.0+1.0im  -53.0+0.0im    34.0+0.0im   24.0+0.0im    34.0+0.0im   36.0+0.0im;
                             24.0+0.0im   28.0+0.0im  -57.0-1.0im   34.5+3.0im   37.0+0.0im  -110.0+0.0im   36.0+0.0im    42.0+0.0im  -32.0+0.0im;
                             10.5-3.0im   12.0+0.0im   18.0+0.0im    9.0+0.0im   14.0+0.0im    21.0+0.0im  -73.5+0.0im   -43.0+1.0im  -13.5+3.0im;
                             24.0+0.0im   34.5-3.0im   36.0+0.0im   28.0+0.0im   37.0+0.0im    42.0+0.0im  -57.0+1.0im  -110.0+0.0im  -32.0+0.0im;
                             36.0+0.0im   42.0+0.0im  -31.5-3.0im   42.0+0.0im   49.0+0.0im   -40.0+0.0im  -31.5+3.0im   -40.0+0.0im  -46.0+0.0im])
    @fact globaloperator(H,[L1,L2]) --> resultnoomega
    @fact globaloperator(H,[L1,L2],1/2) --> resultnoomega/2
    #type test

    #typeerrortest
    @fact_throws ArgumentError, globaloperator(H,[L1,L2],1im)
    @fact_throws ArgumentError, globaloperator(H,[L1,L2],-1)
    @fact_throws ArgumentError, globaloperator(H,[L1,L2],3)
  end

end


facts("User utils") do
  context("classicallindbladoperators") do
    H = [1. 1.+im ; 1.-im im]
    result = [sparse([1.+0im 0 ; 0 0 ]),
              sparse([0.+0im 1.+im ; 0 0 ]),
              sparse([0.+0im 0 ; 1.-im 0 ]),
              sparse([0.+0im 0 ; 0 im ])]
    resultwithepsilon = [sparse([0.+0im 1.+im ; 0 0 ]),
              sparse([0.+0im 0 ; 1.-im 0 ])]
    @fact classicallindbladoperators(H) --> result
    @fact classicallindbladoperators(sparse(H)) --> result
    @fact classicallindbladoperators(H, epsilon=1.1) --> resultwithepsilon
    @fact classicallindbladoperators(sparse(H), epsilon=1.1) --> resultwithepsilon
    #type test
    @fact typeof(classicallindbladoperators(H)) --> Vector{SparseMatrixCSC{Complex128}}
    @fact typeof(classicallindbladoperators(sparse(H))) --> Vector{SparseMatrixCSC{Complex128}}
    #error test
    @fact_throws ArgumentError, classicallindbladoperators(H, -1.)
    @fact_throws ArgumentError, classicallindbladoperators(H, -1im)
  end
end

#=facts("Demoralization user utils") do
  context("partitionsize") do
    partition = [[1,4],[2,3,5],[6],[7,8]]
    @fact (QSWalk.partitionsize(partition)) --> 8
  end

  context("defaultlocalhamiltonian") do
    @fact QSWalk.defaultlocalhamiltonian(1) --> spzeros(Complex128,1,1)
    @fact QSWalk.defaultlocalhamiltonian(3) --> sparse([0. im 0.; -im 0 im; 0 -im 0])
  end

  context("localhamiltonian") do
    #default option
    @fact localhamiltonian([[1],[2,3]]) --> sparse([0.+0im 0 0;0 0 im;0 -im 0])
    @fact localhamiltonian([[1],[2],[3]]) --> spzeros(Complex128,3,3)
    #by size version
    @fact localhamiltonian([[1,3],[2]], x->speye(x)) --> speye(Complex128,3)
    @fact localhamiltonian([[1],[2],[3]], x->speye(x)) --> speye(Complex128,3)
    #by index version
    M1 = [1. 2;3 5]
    M2 = zeros(1,1)+1.
    @fact localhamiltonian([[1,3],[2]], FunctionByIndex(x->[M1,M2][x])) -->
            sparse([ 1.0+0im 0 2; 0 1 0; 3 0 5 ])
  end

  context("Incidences lists") do

  end

  context("makepartition") do
    @fact QSWalk.makepartition([[1,3],[2,3],Int[],[4,6,1]]) --> [[1,2],[3,4],[5],[6,7,8]]
  end

  context("fourier matrix") do
    @fact QSWalk.rectangularfouriermatrix(2,1) --> roughly([1,1])
    @fact QSWalk.rectangularfouriermatrix(1,1) --> roughly([1])
    @fact QSWalk.rectangularfouriermatrix(4,2) --> roughly([1, 1im, -1, -1im])
  end

  context("demoralizedlindbladian") do
    A = sparse([0 1 0; 1 0 1; 0 1 0]+0im)*1.

    function symmetricl(size, column) #działa tylko dla size<=2
      if size == 1
        return [1]
      else
        if column == 1
          return [1,-1]
        else
          return [1, 1]
        end
      end
    end
    #default
    @fact demoralizedlindbladian(A)[1] --> roughly(sparse([0 1 1 0; 1 0 0 1; 1 0 0 -1; 0 1 1 0]))
    @fact demoralizedlindbladian(A)[2] --> QSWalk.makepartition(QSWalk.reversedincidencelist(A))

    #symmetric case
    @fact demoralizedlindbladian(A,eps(1.),(x,y)->symmetricl(x,y))[1] -->
                    roughly(sparse([0 1 1 0; 1 0 0 1; -1 0 0 1; 0 1 1 0]+0im)*1.)
    @fact demoralizedlindbladian(A,eps(1.),(x,y)->symmetricl(x,y))[2] -->
                    QSWalk.makepartition(QSWalk.reversedincidencelist(A))

  end
end

facts("evolution") do
  context("distributionsummation") do
    probability = [0.05,0.1,0.25,0.3,0.01,0.20,0.04,0.05]
    partition = [[1,4],[2,3,5],[6],[7,8]]
    result = [0.35,0.36,0.2,0.09]
    @fact distributionsummation(probability,partition) --> result
  end

  context("Initial states") do
    @fact initialstate([1,3,4],[[1],[2,3,4],[5,6,7],[8,9]]) --> roughly(spdiagm([1./3,0,0,0,1./9,1./9,1./9,1./6,1./6]))
    A1 = ones(1,1)/4
    A2 = [ 1/5 0 1/5; 0 1/10 0 ; 1/5 0 1/5 ]
    A3 = [0.125 -0.125+0im; -0.125 0.125]
    @fact initialstate([A1,A2,A3]) -->
          sparse([0.25+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
           0.0+0.0im 0.2+0.0im 0.0+0.0im 0.2+0.0im 0.0+0.0im 0.0+0.0im;
           0.0+0.0im 0.0+0.0im 0.1+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
           0.0+0.0im 0.2+0.0im 0.0+0.0im 0.2+0.0im 0.0+0.0im 0.0+0.0im;
           0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0.125+0.0im -0.125+0.0im;
           0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im -0.125+0.0im 0.125+0.0im])
  end

  context("Simple evolution functions") do
    #size 4x4
    A = Complex{Float64}[0.354177+0.0im 0.0891553-0.0251879im 0.0702961+0.0516828im 0.0708664+0.0767941im; 0.0891553+0.0251879im 0.336055+0.0im 0.0420202-0.0109173im 0.0683605-0.00692846im; 0.0702961-0.0516828im 0.0420202+0.0109173im 0.212401+0.0im 0.0939615+0.0553555im; 0.0708664-0.0767941im 0.0683605+0.00692846im 0.0939615-0.0553555im 0.0973671+0.0im]
    #expm
    @fact simpleevolve(zeros(16,16), A, 0.) --> roughly(A)
    @fact simpleevolve(zeros(16,16), sparse(A), 0.) --> roughly(A)
    @fact simpleevolve(zeros(16,16), A, [0.,5.,10.]) --> [A,A,A]

    @fact simpleevolve(rand(16,16), A, 0.) --> roughly(A)

    #expmv
    @fact simpleevolve(spzeros(16,16), A, 0.) --> roughly(A)
    @fact simpleevolve(spzeros(16,16), A, [0.,5.,10.])[1] --> roughly(A)
    @fact simpleevolve(spzeros(16,16), A, [0.,5.,10.])[2] --> roughly(A)
    @fact simpleevolve(spzeros(16,16), A, [0.,5.,10.])[3] --> roughly(A)

    @fact simpleevolve(sparse(rand(16,16)), A, 0.) --> roughly(A)
  end
end

facts("User compact functions") do
  context("expmv evolution") do

  end

  context("expm evolution") do

  end

  context("partition") do

  end
end
=#
