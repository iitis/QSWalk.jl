module QSWalk
using Expokit

include("utils.jl")
include("dirac.jl")
include("operator.jl")
include("demoralization.jl")
include("evolution.jl")
include("compactuserfunctions.jl")

export
  SparseDenseMatrix, 
  SparseDenseVector, 
  Vertex, 
  VertexSet, 
  ket, 
  bra, 
  ketbra, 
  proj, 
  res, 
  unres, 
  classical_lindbladian, 
  evolve_generator, 
  local_hamiltonian, 
  nonmoralizing_lindbladian, 
  make_vertex_set, 
  init_nonmoralized, 
  measurement_nonmoralized, 
  evolve_operator, 
  evolve, 
  global_hamiltonian

end
