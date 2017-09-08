
facts("Vertex test") do
  context("Vertex creation") do
    @fact Vertex([1,2,3]).linspace --> [1,2,3]
    @fact Vertex([1,2,3])() --> Vertex([1,2,3]).linspace
  end

  context("Error check") do
    @fact_throws ArgumentError Vertex([1,2,0])
  end
end

facts("Vertex set") do
  context("VertexSet creation") do
    @fact VertexSet([Vertex([1,2,3]),Vertex([4,5])])() --> [Vertex([1,2,3]),Vertex([4,5])]
    @fact VertexSet([Vertex([1,2,3]),Vertex([4,5])])() --> [Vertex([1,2,3]),Vertex([4,5])]
    @fact VertexSet([[1,2,3],[4,5]])() --> [Vertex([1,2,3]),Vertex([4,5])]
  end

  context("Error Check") do
    @fact_throws ArgumentError VertexSet([Vertex([1,2,0])])
    @fact_throws ArgumentError VertexSet([Vertex([1,2,3]),Vertex([3,4])])
  end
end
