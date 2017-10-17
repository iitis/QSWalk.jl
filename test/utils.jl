
facts("Vertex test") do
  context("Vertex creation") do
    @fact subspace(Vertex([1, 2, 3])) --> [1, 2, 3]
    @fact subspace(Vertex([1, 2, 3])) --> Vertex([1, 2, 3]).subspace
  end

  context("Basic Vertex functionality") do
    v = Vertex([3, 4, 5])
    vcopy = Vertex([3, 4, 5])
    v1 = Vertex([1, 2])
    v1copy = Vertex([1, 2])
    @fact subspace(v) --> [3, 4, 5]
    @fact hash(v) --> hash(subspace(v))
    @fact hash((v, vcopy)) --> hash([subspace(v), subspace(vcopy)])
    @fact (v, v1) ==  (vcopy, v1copy) --> true
    @fact v --> vcopy
    @fact v[2] --> 4
    @fact length(v) --> 3
  end

  context("Error check") do
    @fact_throws ArgumentError Vertex([1, 2, 0])
  end
end

facts("Vertex set") do
  context("VertexSet creation") do
    @fact vlist(VertexSet([Vertex([1, 2, 3]), Vertex([4, 5])])) --> [Vertex([1, 2, 3]), Vertex([4, 5])]
    @fact vlist(VertexSet([Vertex([1, 2, 3]), Vertex([4, 5])])) --> [Vertex([1, 2, 3]), Vertex([4, 5])]
    @fact vlist(VertexSet([[1, 2, 3], [4, 5]])) --> [Vertex([1, 2, 3]), Vertex([4, 5])]
  end

  context("Basic VertexSet functionaliy") do
    vset = VertexSet([[1, 2, 3], [4, 5]])
    vset2 = VertexSet([[1, 2, 3], [4, 5]])
    @fact vlist(vset) -->[Vertex([1, 2, 3]), Vertex([4, 5])]
    @fact vset --> vset2
    @fact vset[2] --> Vertex([4, 5])
    @fact vset[[1, 2]] --> vlist(vset)
    @fact length(vset) --> 2
  end

  context("Error Check") do
    @fact_throws ArgumentError VertexSet([Vertex([1, 2, 0])])
    @fact_throws ArgumentError VertexSet([Vertex([1, 2, 3]), Vertex([3, 4])])
  end
end
