
@testset "Vertex test" begin
  @testset "Vertex creation" begin
    @test subspace(Vertex([1, 2, 3])) == [1, 2, 3]
    @test subspace(Vertex([1, 2, 3])) == Vertex([1, 2, 3]).subspace
  end

  @testset "Basic Vertex functionality" begin
    v = Vertex([3, 4, 5])
    vcopy = Vertex([3, 4, 5])
    v1 = Vertex([1, 2])
    v1copy = Vertex([1, 2])
    @test subspace(v) == [3, 4, 5]
    @test hash(v) == hash(subspace(v))
    @test hash((v, vcopy)) == hash([subspace(v), subspace(vcopy)])
    @test (v, v1) == (vcopy, v1copy)
    @test v == vcopy
    @test v[2] == 4
    @test length(v) == 3
  end

  @testset "Error check" begin
    @test_throws ArgumentError Vertex([1, 2, 0])
  end
end

@testset "Vertex set" begin
  @testset "VertexSet creation" begin
    @test vlist(VertexSet([Vertex([1, 2, 3]), Vertex([4, 5])])) == [Vertex([1, 2, 3]), Vertex([4, 5])]
    @test vlist(VertexSet([Vertex([1, 2, 3]), Vertex([4, 5])])) == [Vertex([1, 2, 3]), Vertex([4, 5])]
    @test vlist(VertexSet([[1, 2, 3], [4, 5]])) == [Vertex([1, 2, 3]), Vertex([4, 5])]
  end

  @testset "Basic VertexSet functionaliy" begin
    vset = VertexSet([[1, 2, 3], [4, 5]])
    vset2 = VertexSet([[1, 2, 3], [4, 5]])
    @test vlist(vset) ==[Vertex([1, 2, 3]), Vertex([4, 5])]
    @test vset == vset2
    @test vset[2] == Vertex([4, 5])
    @test vset[[1, 2]] == vlist(vset)
    @test length(vset) == 2
  end

  @testset "Error Check" begin
    @test_throws ArgumentError VertexSet([Vertex([1, 2, 0])])
    @test_throws ArgumentError VertexSet([Vertex([1, 2, 3]), Vertex([3, 4])])
  end
end
