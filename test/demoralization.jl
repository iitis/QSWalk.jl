@testset "Demoralization user utils" begin
  @testset "vertexsetsize" begin
    partition_int = [[1, 4], [2, 3, 5], [6], [7, 8]]
    partition_vertex = VertexSet([Vertex(block) for block = partition_int])
    @test (QSWalk.vertexsetsize(partition_vertex)) == 8
    # error test
    @test_throws MethodError QSWalk.vertexsetsize(partition_int)

  end

  @testset "default_nm_loc_ham" begin
    @test QSWalk.default_nm_loc_ham(1) == spzeros(Complex128, 1, 1)
    @test QSWalk.default_nm_loc_ham(3) == sparse([0. im 0.; -im 0 im; 0 -im 0])
  end

  @testset "nm_loc_ham" begin
    #default option
    @test nm_loc_ham(VertexSet([[1], [2, 3]])) == sparse([0.+0im 0 0;0 0 im;0 -im 0])
    @test nm_loc_ham(VertexSet([[1], [2], [3]])) == spzeros(Complex128, 3, 3)
    #by size version
    @test nm_loc_ham(VertexSet([[1, 2], [3, 4]]), Dict(2 =>[0 1; 1 0])) == sparse([0. 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0])
    @test nm_loc_ham(VertexSet([[1], [2], [3]]), Dict(1 =>ones(Float64, (1, 1)))) == speye(Complex128, 3)
    #by index version
    M1 = [1. 2;3 5]
    M2 = zeros(1, 1)+1.
    dict = Dict(Vertex([1, 3]) =>M1, Vertex([2]) =>M2)
    @test nm_loc_ham(VertexSet([[1, 3], [2]]), dict) ==
            sparse([ 1.0+0im 0 2; 0 1 0; 3 0 5 ])
  end

  @testset "Incidence and reverse incidence lists" begin
    A = [1 2 3; 0 3. 4.; 0 0 5.]
    @test QSWalk.incidence_list(A) == [[1], [1, 2], [1, 2, 3]]
    @test QSWalk.incidence_list(sparse(A)) == QSWalk.incidence_list(A)
    @test QSWalk.incidence_list(A; epsilon = 2.5) == [Int64[], [2], [1, 2, 3]]
    @test QSWalk.incidence_list(sparse(A)) == QSWalk.incidence_list(A)
    @test QSWalk.reversed_incidence_list(A) == [[1, 2, 3], [2, 3], [3]]
    @test QSWalk.reversed_incidence_list(sparse(A)) == QSWalk.reversed_incidence_list(A)
    @test QSWalk.reversed_incidence_list(A; epsilon = 2.5) == [[3], [2, 3], [3]]
    @test QSWalk.reversed_incidence_list(sparse(A); epsilon = 2.5) == QSWalk.reversed_incidence_list(A; epsilon = 2.5)
    #errors
    @test_throws ArgumentError QSWalk.reversed_incidence_list(A, epsilon = -1)
    @test_throws ArgumentError QSWalk.incidence_list(A, epsilon = -1)
    @test_throws ArgumentError QSWalk.reversed_incidence_list(sparse(A), epsilon = -1)
    @test_throws ArgumentError QSWalk.incidence_list(sparse(A), epsilon = -1)
  end

  @testset "vertexset creation" begin
    @test QSWalk.revinc_to_vertexset([[1, 3], [2, 3], Int[], [4, 6, 1]]) == VertexSet([[1, 2], [3, 4], [5], [6, 7, 8]])
    @test make_vertex_set([1 2 3; 0 3. 4.; 0 0 5.]) == VertexSet([[1, 2, 3], [4, 5], [6]])
  end

  @testset "Fourier matrix" begin
    @test QSWalk.fourier_matrix(4) ≈ [1 1 1 1; 1 1im -1 -1im;1 -1 1 -1; 1 -1im -1 1im]
    @test QSWalk.fourier_matrix(2) ≈ [1 1; 1 -1.]
    @test QSWalk.fourier_matrix(1) ≈ ones(Float64, (1, 1))
  end

  @testset "nm_lind" begin
    A = sparse([0.+0.0im 1 0; 1 0 1; 0 1 0])
    B1, B2 = eye(1), ones(2, 2)
    #default
    #needs to be roughly, since exp computing is inexact
    @test nm_lind(A)[1] ≈ [0 1 1 0;
                                             1 0 0 1;
                                             1 0 0 -1;
                                             0 1 1 0]
    @test nm_lind(A)[2] == make_vertex_set(A)

    A = sparse([0.+0.0im 0 0 0 1;
                0 0 1 0 1;
                0 0 0 0 0;
                0 1 1 0 0;
                0 0 0 0 0])
    L, vset = nm_lind(A)
    v1, v2, v3, v4, v5 = vlist(vset)
    @test nm_lind(A)[1] ≈ [0 0 0 0 0 0 1;
                                             0 0 0 1 0 0 1;
                                             0 0 0 1 0 0 -1;
                                             0 0 0 0 0 0 0;
                                             0 1 1 1 0 0 0;
                                             0 1 1 -1 0 0 0;
                                             0 0 0 0 0 0 0]
    @test nm_lind(A)[2] == make_vertex_set(A)
    @test nm_lind(A, Dict(v1 => B1, v2 => 2*B2, v3 => 3*B1, v4 => 4*B2, v5 => 5*B1))[1] ≈
             [0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  4.0+0.0im  4.0+0.0im  4.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  4.0+0.0im  4.0+0.0im  4.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]



    A = [0 0 0; 0 0 0; 1 1 0]
    @test nm_lind(A)[1] ≈ [0 0 0 0;
                                             0 0 0 0;
                                             1 1 0 0;
                                             1 -1 0 0]
    @test nm_lind(A)[2] == make_vertex_set(A)


    @test nm_lind(A, Dict(1  => B1, 2 =>3*B2 ))[1] ≈
             [0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im]
    @test nm_lind(A, Dict(1  => B1, 2 =>3*B2 ))[2] ==
      make_vertex_set(A)

    @test nm_lind(A, Dict(1  => B1, 2 =>3*B2 ))[1] ≈
             [0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im]

    L, vset = nm_lind(A)
    v1, v2, v3 = vlist(vset)
    @test nm_lind(A, Dict(v1 => B1, v2 => 2*B1, v3 => 3*B2))[1] ≈
             [0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im;
              3.0+0.0im  3.0+0.0im  0.0+0.0im  0.0+0.0im]
    @test nm_lind(A, Dict(v1 => B1, v2 => 2*B1, v3 => 3*B2))[2] ==
      make_vertex_set(A)

      A = [2 0 0 3;
    im 0 3im 0;
    1 0 0 0;
    im -im 0 0]
    @test nm_lind(A,Dict(1=>ones(1,1), 2=> ones(2,2)))[1] ≈
              [2.0+0.0im 2.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 3.0+0.0im 3.0+0.0im;
               2.0+0.0im 2.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 3.0+0.0im 3.0+0.0im;
               0.0+1.0im 0.0+1.0im 0.0+0.0im 0.0+0.0im 0.0+3.0im 0.0+0.0im 0.0+0.0im;
               0.0+1.0im 0.0+1.0im 0.0+0.0im 0.0+0.0im 0.0+3.0im 0.0+0.0im 0.0+0.0im;
               1.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
               0.0+1.0im 0.0+1.0im 0.0-1.0im 0.0-1.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
               0.0+1.0im 0.0+1.0im 0.0-1.0im 0.0-1.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im]


  end

  @testset "nm_glob_ham" begin
    A = sparse([0 1 0;
               1 0 2;
               0 2 0])
    #default
    #needs to be roughly, since exp computing is inexact
    @test nm_glob_ham(A) ≈ [0 1 1 0;
                                   1 0 0 2;
                                   1 0 0 2;
                                   0 2 2 0]
    @test nm_glob_ham(A, Dict((1, 2) => (2+1im)*ones(1, 2), (2, 1) =>1im*ones(2, 1))) ≈
                                    [0 2+1im 2+1im 0;
                                     2-1im 0 0 2im;
                                     2-1im 0 0 2im;
                                     0 -2im -2im 0]

    v1, v2, v3 = vlist(make_vertex_set(A))
    @test nm_glob_ham(A, Dict((v1, v2) =>2*ones(1, 2), (v2, v3) =>ctranspose([1im 2im;]))) ==
                                              [0 2 2 0;
                                               2 0 0 2im;
                                               2 0 0 4im;
                                               0 -2im -4im 0]
   A = [0  0  1;
        0  0  1;
        2  2  0]
  v1, v2, v3 = vlist(make_vertex_set(A))
  @test nm_glob_ham(A, Dict((v1, v3) =>2*ones(1, 2), (v2, v3) =>[1im 2im;])) ≈
                          [0.0+0.0im  0.0+0.0im  2.0+0.0im  2.0+0.0im;
                           0.0+0.0im  0.0+0.0im  0.0-1.0im  0.0-2.0im;
                           2.0+0.0im  0.0+1.0im  0.0+0.0im  0.0+0.0im;
                           2.0+0.0im  0.0+2.0im  0.0+0.0im  0.0+0.0im]

  @test nm_glob_ham(A, Dict((v1, v3) =>2*ones(1, 2), (v2, v3) =>[1im 2im;]),epsilon=1.5) ≈
                          [0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
                           0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
                           0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im;
                           0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]

  end
end


@testset "Cannonical measurement" begin
  probability = [0.05, 0.1, 0.25, 0.3, 0.01, 0.20, 0.04, 0.05]
  vset = VertexSet([[1, 4], [2, 3, 5], [6], [7, 8]])
  result = [0.35, 0.36, 0.2, 0.09]
  @test nm_measurement(probability, vset) == result

  state = [1/6 1.0 1/6; 1. 2/3 1im; 1/6 -1im 1/6]
  vset = VertexSet([[1, 3], [2]])
  @test nm_measurement(state, vset) ≈ [1./3, 2./3]
  @test_throws AssertionError nm_measurement(eye(2)/2., vset)
  @test_throws AssertionError nm_measurement([1./2, 1./3, 1./6, 0.], vset)

end

@testset "Initial states creation" begin
  vset = VertexSet([[1], [2, 3, 4], [5, 6, 7], [8, 9]])
  @test nm_init(vset[[1, 3, 4]], vset) ≈ spdiagm([1./3, 0, 0, 0, 1./9, 1./9, 1./9, 1./6, 1./6])
  A1 = ones(Complex128, 1, 1)/4
  A2 = [ 1/5+0im 0 1/5; 0 1/10 0 ; 1/5 0 1/5 ]
  A3 = [0.125 -0.125+0im; -0.125 0.125]
  dict = Dict(vset[1] =>A1, vset[3] =>A2, vset[4] =>A3 )
  @test nm_init(dict, vset) ==
        sparse([ 0.25+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.2+0.0im  0.0+0.0im  0.2+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.1+0.0im  0.0+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.2+0.0im  0.0+0.0im  0.2+0.0im     0.0+0.0im     0.0+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im   0.125+0.0im  -0.125+0.0im
                0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  -0.125+0.0im   0.125+0.0im])

  dictwrong = Dict(vset[1] =>A2, vset[3] =>A2, vset[4] =>A3 )
  @test_throws AssertionError nm_init(vlist(VertexSet([[1, 2, 3], [4, 5]])), vset)
  @test_throws AssertionError nm_init(dictwrong, vset)
end
