facts("Demoralization user utils") do
  context("vertexsetsize") do
    partition_int = [[1, 4], [2, 3, 5], [6], [7, 8]]
    partition_vertex = VertexSet([Vertex(block) for block = partition_int])
    @fact (QSWalk.vertexsetsize(partition_vertex)) --> 8
    # error test
    @fact_throws MethodError QSWalk.vertexsetsize(partition_int)

  end

  context("default_local_hamiltonian") do
    @fact QSWalk.default_local_hamiltonian(1) --> spzeros(Complex128, 1, 1)
    @fact QSWalk.default_local_hamiltonian(3) --> sparse([0. im 0.; -im 0 im; 0 -im 0])
  end

  context("local_hamiltonian") do
    #default option
    @fact local_hamiltonian(VertexSet([[1], [2, 3]])) --> sparse([0.+0im 0 0;0 0 im;0 -im 0])
    @fact local_hamiltonian(VertexSet([[1], [2], [3]])) --> spzeros(Complex128, 3, 3)
    #by size version
    @fact local_hamiltonian(VertexSet([[1, 2], [3, 4]]), Dict(2 =>[0 1; 1 0])) --> sparse([0. 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0])
    @fact local_hamiltonian(VertexSet([[1], [2], [3]]), Dict(1 =>ones(Float64, (1, 1)))) --> speye(Complex128, 3)
    #by index version
    M1 = [1. 2;3 5]
    M2 = zeros(1, 1)+1.
    dict = Dict(Vertex([1, 3]) =>M1, Vertex([2]) =>M2)
    @fact local_hamiltonian(VertexSet([[1, 3], [2]]), dict) -->
            sparse([ 1.0+0im 0 2; 0 1 0; 3 0 5 ])
  end

  context("Incidence and reverse incidence lists") do
    A = [1 2 3; 0 3. 4.; 0 0 5.]
    @fact QSWalk.incidence_list(A) --> [[1], [1, 2], [1, 2, 3]]
    @fact QSWalk.incidence_list(sparse(A)) --> QSWalk.incidence_list(A)
    @fact QSWalk.incidence_list(A; epsilon = 2.5) --> [Int64[], [2], [1, 2, 3]]
    @fact QSWalk.incidence_list(sparse(A)) --> QSWalk.incidence_list(A)
    @fact QSWalk.reversed_incidence_list(A) --> [[1, 2, 3], [2, 3], [3]]
    @fact QSWalk.reversed_incidence_list(sparse(A)) --> QSWalk.reversed_incidence_list(A)
    @fact QSWalk.reversed_incidence_list(A; epsilon = 2.5) --> [[3], [2, 3], [3]]
    @fact QSWalk.reversed_incidence_list(sparse(A); epsilon = 2.5) --> QSWalk.reversed_incidence_list(A; epsilon = 2.5)
    #errors
    @fact_throws ArgumentError QSWalk.reversed_incidence_list(A, epsilon = -1)
    @fact_throws ArgumentError QSWalk.incidence_list(A, epsilon = -1)
    @fact_throws ArgumentError QSWalk.reversed_incidence_list(sparse(A), epsilon = -1)
    @fact_throws ArgumentError QSWalk.incidence_list(sparse(A), epsilon = -1)
  end

  context("vertexset creation") do
    @fact QSWalk.revinc_to_vertexset([[1, 3], [2, 3], Int[], [4, 6, 1]]) --> VertexSet([[1, 2], [3, 4], [5], [6, 7, 8]])
    @fact make_vertex_set([1 2 3; 0 3. 4.; 0 0 5.]) --> VertexSet([[1, 2, 3], [4, 5], [6]])
  end

  context("Fourier matrix") do
    @fact QSWalk.fourier_matrix(4) --> roughly([1 1 1 1; 1 1im -1 -1im;1 -1 1 -1; 1 -1im -1 1im])
    @fact QSWalk.fourier_matrix(2) --> roughly([1 1; 1 -1.])
    @fact QSWalk.fourier_matrix(1) --> roughly(ones(Float64, (1, 1)))
  end

  context("nonmoralizing_lindbladian") do
    A = sparse([0.+0.0im 1 0; 1 0 1; 0 1 0])
    B1, B2 = eye(1), ones(2, 2)
    #default
    #needs to be roughly, since exp computing is inexact
    @fact nonmoralizing_lindbladian(A)[1] --> roughly([0 1 1 0;
                                                     1 0 0 1;
                                                     1 0 0 -1;
                                                     0 1 1 0])
    @fact nonmoralizing_lindbladian(A)[2] --> make_vertex_set(A)

    A = sparse([0.+0.0im 0 0 0 1;
                0 0 1 0 1;
                0 0 0 0 0;
                0 1 1 0 0;
                0 0 0 0 0])
    L, vset = nonmoralizing_lindbladian(A)
    v1, v2, v3, v4, v5 = vlist(vset)
    @fact nonmoralizing_lindbladian(A)[1] --> roughly([0 0 0 0 0 0 1;
                                                     0 0 0 1 0 0 1;
                                                     0 0 0 1 0 0 -1;
                                                     0 0 0 0 0 0 0;
                                                     0 1 1 1 0 0 0;
                                                     0 1 1 -1 0 0 0;
                                                     0 0 0 0 0 0 0])
    @fact nonmoralizing_lindbladian(A)[2] --> make_vertex_set(A)
    @fact nonmoralizing_lindbladian(A, Dict(v1 => B1, v2 => 2*B2, v3 => 3*B1, v4 => 4*B2, v5 => 5*B1))[1] -->
      roughly([0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  4.0+0.0im  4.0+0.0im  4.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  4.0+0.0im  4.0+0.0im  4.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im])



    A = [0 0 0; 0 0 0; 1 1 0]
    @fact nonmoralizing_lindbladian(A)[1] --> roughly([0 0 0 0;
                                                     0 0 0 0;
                                                     1 1 0 0;
                                                     1 -1 0 0])
    @fact nonmoralizing_lindbladian(A)[2] --> make_vertex_set(A)


    @fact nonmoralizing_lindbladian(A, Dict(1  => B1, 2 =>3*B2 ))[1] -->
      roughly([0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im])
    @fact nonmoralizing_lindbladian(A, Dict(1  => B1, 2 =>3*B2 ))[2] -->
      make_vertex_set(A)

    @fact nonmoralizing_lindbladian(A, Dict(1  => B1, 2 =>3*B2 ))[1] -->
      roughly([0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im])

    L, vset = nonmoralizing_lindbladian(A)
    v1, v2, v3 = vlist(vset)
    @fact nonmoralizing_lindbladian(A, Dict(v1 => B1, v2 => 2*B1, v3 => 3*B2))[1] -->
      roughly([0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im])
    @fact nonmoralizing_lindbladian(A, Dict(v1 => B1, v2 => 2*B1, v3 => 3*B2))[2] -->
      make_vertex_set(A)

  end

  context("global_hamiltonian") do
    A = sparse([0 1 0; 1 0 2; 0 2 0])
    #default
    #needs to be roughly, since exp computing is inexact
    @fact global_hamiltonian(A) --> roughly([0 1 1 0;
                                             1 0 0 2;
                                             1 0 0 2;
                                             0 2 2 0])
    @fact global_hamiltonian(A, Dict((1, 2) => (2+1im)*ones(1, 2), (2, 1) =>1im*ones(2, 1))) -->
                                    roughly([0 2+1im 2+1im 0;
                                             2-1im 0 0 2im;
                                             2-1im 0 0 2im;
                                             0 -2im -2im 0])

    v1, v2, v3 = vlist(make_vertex_set(A))
    @fact global_hamiltonian(A, Dict((v1, v2) =>2*ones(1, 2), (v2, v3) =>[1im 2im;]')) -->
                                              [0 2 2 0;
                                               2 0 0 2im;
                                               2 0 0 4im;
                                               0 -2im -4im 0]
  end
end
