var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Welcome-to-QSWalk.jl-1",
    "page": "Home",
    "title": "Welcome to QSWalk.jl",
    "category": "section",
    "text": "QSWalk.jl is a package providing implementation of open continuous-time quantum walk evolution based on GKSL master equation. In particular we have implemented function useful for analysing local interaction, global interaction and nonmoralizing global interaction stochastic quantum walk models.In GitHub we present an examples in .ipynb and .jl extension. The present most of the functionalities of the package. Furthermore, for package explanation we refer the user to our paper describing the package, and for the basic papers describing and introducing quantum stochastic models:Quantum stochastic walks: A generalization of classical random walks and quantum walks by Whitfield, Rodríguez-Rosario, and Aspuru-Guzik,\nSuperdiffusive quantum stochastic walk definable on arbitrary directed graph by Domino, Glos, and Ostaszewski,\nProperties of quantum stochastic walks from the asymptotic scaling exponent by Domino, Glos, Ostaszewski, Pawela, and Sadowski."
},

{
    "location": "gksl.html#",
    "page": "GKSL master equation",
    "title": "GKSL master equation",
    "category": "page",
    "text": "DocTestSetup = quote\n   using QSWalk\nend"
},

{
    "location": "gksl.html#GKSL-master-equation-1",
    "page": "GKSL master equation",
    "title": "GKSL master equation",
    "category": "section",
    "text": "GKSL master equation is general continuous-time open quantum evolution. The master equation base on to form of operators: the Hamiltonian, which describes the evolution of closed system, and the Lindbladian operators which describes the evolution of open system. Basic facts in context of quantum stochastic walks can found here. For local interaction, where each Lindblad operator is a matrix with single nonzero element, we created a local_lind function which splits matrix into mentioned operators.Below we present a documentation of basic function used for simulating GKSL master equation.Order = [:type, :function]\nModules = [QSWalk]\nPages   = [\"gksl.md\"]"
},

{
    "location": "gksl.html#QSWalk.ket",
    "page": "GKSL master equation",
    "title": "QSWalk.ket",
    "category": "function",
    "text": "ket(index, size)\n\nReturn index-th base (column) vector in the size-dimensional vector space. To be consistent with Julia indexing, index = 1, 2, ..., size.\n\nExamples\n\njulia> ket(1, 2)\n2-element SparseArrays.SparseVector{Int64,Int64} with 1 stored entry:\n  [1]  =  1\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.bra",
    "page": "GKSL master equation",
    "title": "QSWalk.bra",
    "category": "function",
    "text": "bra(index, size)\n\nReturn index-th base row vector in the size-dimensional vector space, with index = 1, 2, ..., size.\n\nExamples\n\njulia> bra(1, 2)\n1×2 LinearAlgebra.Adjoint{Int64,SparseArrays.SparseVector{Int64,Int64}}:\n 1  0\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.ketbra",
    "page": "GKSL master equation",
    "title": "QSWalk.ketbra",
    "category": "function",
    "text": "ketbra(irow, icol, size)\n\nReturn matrix acting on size-dimensional vector space, irow, icol = 1, 2, ..., size. The matrix consists of single non-zero element equal to one, located at position (irow, icol).\n\nExamples\n\njulia> ketbra(1, 2, 3)\n3×3 SparseArrays.SparseMatrixCSC{Int64,Int64} with 1 stored entry:\n  [1, 2]  =  1\n\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.proj",
    "page": "GKSL master equation",
    "title": "QSWalk.proj",
    "category": "function",
    "text": "proj(index, size)\n\nReturn projector onto index-th base vector in size-dimensional vector space, with index = 1, 2, ..., size. This is equivalent to ketbra(index, index, size).\n\njulia> proj(1, 2)\n2×2 SparseArrays.SparseMatrixCSC{Int64,Int64} with 1 stored entry:\n  [1, 1]  =  1\n\n\n\n\n\nproj(vector)\n\nReturn projector onto the subspace spanned by vector vector.\n\nExamples\n\njulia> v = 1/sqrt(2) * (ket(1, 3)+ket(3, 3))\n3-element SparseArrays.SparseVector{Float64,Int64} with 2 stored entries:\n  [1]  =  0.707107\n  [3]  =  0.707107\n\njulia> proj(v)\n3×3 SparseArrays.SparseMatrixCSC{Float64,Int64} with 4 stored entries:\n  [1, 1]  =  0.5\n  [3, 1]  =  0.5\n  [1, 3]  =  0.5\n  [3, 3]  =  0.5\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.res",
    "page": "GKSL master equation",
    "title": "QSWalk.res",
    "category": "function",
    "text": "res(matrix)\n\nReturn vectorization of the matrix in the row order. This is equivalent to Base.vec(transpose(matrix)).\n\nExamples\n\njulia> M = reshape(1:9, (3, 3))\'*1im\n3×3 Array{Complex{Int64},2}:\n 0+1im  0+2im  0+3im\n 0+4im  0+5im  0+6im\n 0+7im  0+8im  0+9im\n\njulia> v = res(M)\n9-element reshape(::LinearAlgebra.Transpose{Complex{Int64},Array{Complex{Int64},2}}, 9) with eltype Complex{Int64}:\n 0 + 1im\n 0 + 2im\n 0 + 3im\n 0 + 4im\n 0 + 5im\n 0 + 6im\n 0 + 7im\n 0 + 8im\n 0 + 9im\n\njulia> res(unres(v)) == v\ntrue\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.unres",
    "page": "GKSL master equation",
    "title": "QSWalk.unres",
    "category": "function",
    "text": "unres(vector)\n\nReturn square matrix with elements from vector. The vector is expected to have perfect square number of arguments.\n\nExamples\n\njulia> v = collect(1:9)*im\n9-element Array{Complex{Int64},1}:\n 0 + 1im\n 0 + 2im\n 0 + 3im\n 0 + 4im\n 0 + 5im\n 0 + 6im\n 0 + 7im\n 0 + 8im\n 0 + 9im\n\njulia> unres(v)\n3×3 LinearAlgebra.Transpose{Complex{Int64},Array{Complex{Int64},2}}:\n 0+1im  0+2im  0+3im\n 0+4im  0+5im  0+6im\n 0+7im  0+8im  0+9im\n\njulia> res(unres(v)) == v\ntrue\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.evolve_generator-Tuple{AbstractArray{#s1,2} where #s1<:Number,AbstractArray{#s2,1} where #s2<:(AbstractArray{#s3,2} where #s3<:Number),AbstractArray{#s4,2} where #s4<:Number,Real}",
    "page": "GKSL master equation",
    "title": "QSWalk.evolve_generator",
    "category": "method",
    "text": "evolve_generator(H, L[, localH][, ω])\n\nCreate the generator for the evolution superoperator. Given Hamiltonian H, collection of Lindblad operator L, local Hamiltonian localH and scaling parameter ω, the generator is obtained as a sum\n\n-i(1-ω) (H  1 - 1  H) + ω (-i(localH  1 - 1  localH) + (L  L - 12(L^L  1 + 1  L^T L )))\n\nThe last two arguments are optional.\n\nIf localH is not given, it defaults to sparse zero matrix of the size of H.\n\nIf ω is not given, both parts are taken with the same intensity and the global operator takes the form\n\n-i(H  1 - 1  H) + (-i(localH  1 - 1  localH) + (L  L - 12(L^L  1 + 1  L^T L )))\n\nArguments\n\nH: Hamiltonian, must be hermitian,\nL: collection of Lindblad operators, each must be of the same size as H,\nlocalH: local Hamiltonian, suggested for nonmoralized QS walk, must be hermitian and of the size of H,\nω: scaling parameter, should be in [0, 1].\n\nReturn\n\nThe generator matrix, which can be used in evolve function.\n\nExamples\n\njulia> H, L, localH = [0 1+im; 1-im 0], [0. 1; 0 0], [1. 0.; 0. 1.]\n(Complex{Int64}[0+0im 1+1im; 1-1im 0+0im], [0.0 1.0; 0.0 0.0], [1.0 0.0; 0.0 1.0])\n\njulia> evolve_generator(H, [L], localH, 1/2)\n4×4 Array{Complex{Float64},2}:\n  0.0+0.0im    0.5+0.5im    0.5-0.5im   0.5+0.0im\n -0.5+0.5im  -0.25+0.0im    0.0+0.0im   0.5-0.5im\n -0.5-0.5im    0.0+0.0im  -0.25+0.0im   0.5+0.5im\n  0.0+0.0im   -0.5-0.5im   -0.5+0.5im  -0.5+0.0im\n\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.evolve_operator",
    "page": "GKSL master equation",
    "title": "QSWalk.evolve_operator",
    "category": "function",
    "text": "evolve_operator(evo_gen, time)\n\nReturn an exponent of time×evo_gen. This function is useful in the case of multiple initial states and fixed evo_gen and time, as it is faster to compute the exponent once and use it for evolving on different initial states.\n\nNote: Parameter evo_gen must be of type Matrix. For type SparseMatrixCSC case different numerical approach is used. See function epmv in package Expokit.\n\nExamples\n\njulia> H, L = [0 1; 1 0], [[0 1; 0 0], [0 0; 1 0]]\n([0 1; 1 0], Array{Int64,2}[[0 1; 0 0], [0 0; 1 0]])\n\njulia> evolve_operator(evolve_generator(H, L), 4.0)\n4×4 Array{Complex{Float64},2}:\n 0.499815+0.0im                0.0+0.00127256im  …  0.500185+0.0im\n      0.0+0.00127256im  0.00960957+0.0im                 0.0-0.00127256im\n      0.0-0.00127256im  0.00870607+0.0im                 0.0+0.00127256im\n 0.500185+0.0im                0.0-0.00127256im     0.499815+0.0im\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.evolve-Tuple{AbstractArray{#s4,2} where #s4<:Number,AbstractArray{#s3,2} where #s3<:Number}",
    "page": "GKSL master equation",
    "title": "QSWalk.evolve",
    "category": "method",
    "text": "evolve(evogen, state, time)\nevolve(evogen, state, tpoints)\nevolve(evosuper, state)\n\nSimulate the GKSL master equation according to the equation\n\nresult = exp(time*evogen)state\n\nwhere  denotes the vectorization.\n\nNote: The function returns unvectorized result.\n\nThe evolution can be calculated using three different approaches.\n\nIn the simplest case the function accepts matrix evogen specifying the generator of the evolution, state describing the starting point of the evolution, and time specifying the time of the evolution.\n\nNote: If evogen is of type Matrix, the exponent is calculated using exp function. If evogen is of type SparseMatrixCSC, expmv from Expokit.jl is used.\n\nAlternatively, a list of point of time (tpoints) can be given. Points of time needs to be non-negative (you cannot go back in time). In this case a list of resulting states is returned.\n\nNote: It is up to the user to provide evogen and state fulfilling the appropriate conditions. For the procedure to work correctly evogen should be generated by evolve_generator function and state should be a proper density matrix.\n\nThe third approach can be used if the superoperator is known. In this case argument evosuper can be specified. This argument can be generated by evolve_operator function. This is useful to simulate a fixed model of evolution in the case of multiple initial states and the same time point.\n\nExamples\n\njulia> H, L = [0 1; 1 0], [[0 1; 0 0], [0 0; 1 0]]\n([0 1; 1 0], Array{Int64,2}[[0 1; 0 0], [0 0; 1 0]])\n\njulia> Matrix(evolve(evolve_generator(H, L), proj(1, 2), 4.))\n2×2 Array{Complex{Float64},2}:\n 0.499815+0.0im              0.0+0.00127256im\n      0.0-0.00127256im  0.500185+0.0im\n\njulia> Matrix.(evolve(evolve_generator(H, L), proj(1, 2), [1., 2., 3., 4.]))\n4-element Array{Array{Complex{Float64},2},1}:\n [0.433203+0.0im 0.0+0.107605im; 0.0-0.107605im 0.566797+0.0im]\n [0.485766+0.0im 0.0-0.0171718im; 0.0+0.0171718im 0.514234+0.0im]\n [0.505597+0.0im 0.0-0.00261701im; 0.0+0.00261701im 0.494403+0.0im]\n [0.499815+0.0im 0.0+0.00127256im; 0.0-0.00127256im 0.500185+0.0im]\n\njulia> ev_op = evolve_operator(evolve_generator(H, L), 4.)\n4×4 Array{Complex{Float64},2}:\n 0.499815+0.0im                0.0+0.00127256im  …  0.500185+0.0im\n      0.0+0.00127256im  0.00960957+0.0im                 0.0-0.00127256im\n      0.0-0.00127256im  0.00870607+0.0im                 0.0+0.00127256im\n 0.500185+0.0im                0.0-0.00127256im     0.499815+0.0im\n\njulia> Matrix(evolve(ev_op, proj(1, 2)))\n2×2 Array{Complex{Float64},2}:\n 0.499815+0.0im              0.0+0.00127256im\n      0.0-0.00127256im  0.500185+0.0im\n\n\n\n\n\n"
},

{
    "location": "gksl.html#QSWalk.local_lind",
    "page": "GKSL master equation",
    "title": "QSWalk.local_lind",
    "category": "function",
    "text": "local_lind(A[; epsilon])\n\nSplit the elements of the matrix A into a collection of sparse matrices with exactly one non-zero element. Martices are created if the absolute value of the nonzero element is there are not smaller than epsilon, where epsilon should be nonnegative. The epsilon defaults to eps() if not specified.\n\nExamples\n\njulia> A = [1. 2.; 3. 4.]\n2×2 Array{Float64,2}:\n 1.0  2.0\n 3.0  4.0\n\njulia> local_lind(A)\n4-element Array{SparseArrays.SparseMatrixCSC{Float64,Ti} where Ti<:Integer,1}:\n\n  [1, 1]  =  1.0\n\n  [1, 2]  =  2.0\n\n  [2, 1]  =  3.0\n\n  [2, 2]  =  4.0\n\njulia> local_lind(A, epsilon = 1.5)\n3-element Array{SparseArrays.SparseMatrixCSC{Float64,Ti} where Ti<:Integer,1}:\n\n  [1, 2]  =  2.0\n\n  [2, 1]  =  3.0\n\n  [2, 2]  =  4.0\n\n\n\n\n\n\n"
},

{
    "location": "gksl.html#Full-docs-1",
    "page": "GKSL master equation",
    "title": "Full docs",
    "category": "section",
    "text": "ket\nbra\nketbra\nproj\nres\nunres\nevolve_generator(::AbstractMatrix{<:Number}, ::AbstractVector{<:AbstractMatrix{<:Number}}, ::AbstractMatrix{<:Number}, ::Real)\nevolve_operator\nevolve(::AbstractMatrix{<:Number}, ::AbstractMatrix{<:Number})\nlocal_lind"
},

{
    "location": "demoralization.html#",
    "page": "Demoralization",
    "title": "Demoralization",
    "category": "page",
    "text": "DocTestSetup = quote\n   using QSWalk\nend"
},

{
    "location": "demoralization.html#Nonmoralizing-model-1",
    "page": "Demoralization",
    "title": "Nonmoralizing model",
    "category": "section",
    "text": "Global interaction quantum stochastic walk suffers for creating additional connections. For removing such effect, nonmoralizing quantum stochastic walk was introduced, see here. Such model is constructed in several steps. First, the dimensionality is increased, hence to each vertex multidimensional subspaces is attached. Then, Hamiltonian and Lindblad operator is increased, furthermore additional Hamiltonian called \"local Hamiltonian\". is introduced.Below we present additional functionalities typical for nonmoralizing quantum stochastic walk. By default the the operator is generalized as in the original paper.Order = [:type, :function]\nModules = [QSWalk]\nPages   = [\"demoralization.md\"]"
},

{
    "location": "demoralization.html#QSWalk.Vertex",
    "page": "Demoralization",
    "title": "QSWalk.Vertex",
    "category": "type",
    "text": "type Vertex\n\nType consisting of list of Int, describing the labels of vectors from the canonical basis corresponding to the Vertex. To get the vector label one can use Vertex() function, or Vertex[i] for a unique label.\n\nSee [1] for the more information and usage exmaples.\n\n[1] K. Domino, A. Glos, M. Ostaszewski, Superdiffusive quantum stochastic walk definable on arbitrary directed graph, Quantum Information & Computation, Vol.17 No.11&12, pp. 0973-0986, arXiv:1701.04624.\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.VertexSet",
    "page": "Demoralization",
    "title": "QSWalk.VertexSet",
    "category": "type",
    "text": "type VertexSet\n\nType consisting of a list of Vertex objects. It describes the partition of the linear subspace. Object of this type should be constructed using make_vertex_set or by nm_lind functions. In order to get a list of the vertices from an object of type vertexset, one should use vlist() function, or, for a specific Vertex, an getindex function vertexset[i].\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.nm_measurement",
    "page": "Demoralization",
    "title": "QSWalk.nm_measurement",
    "category": "function",
    "text": "nm_measurement(probability, vertexset)\n\nReturn joint probability of probability, which is real-valued probability vector according to vertexset.\n\nNote: It is up to user to provide proper probability vector.\n\nExamples\n\njulia> probability = [0.05, 0.1, 0.25, 0.3, 0.01, 0.20, 0.04, 0.05]\n8-element Array{Float64,1}:\n 0.05\n 0.1\n 0.25\n 0.3\n 0.01\n 0.2\n 0.04\n 0.05\n\njulia> nm_measurement(probability, VertexSet([[1, 4], [2, 3, 5], [6], [7, 8]]))\n4-element Array{Float64,1}:\n 0.35\n 0.36\n 0.2\n 0.09\n\n\n\n\n\nnm_measurement(state, vertexset)\n\nReturn joint probability of cannonical measurement of density matrix state, according to vertexset.\n\nNote: It is up to user to provide proper density state.\n\nExamples\n\njulia> state = [1/6 0 1/6; 0 2/3 0; 1/6 0 1/6]\n3×3 Array{Float64,2}:\n 0.166667  0.0       0.166667\n 0.0       0.666667  0.0\n 0.166667  0.0       0.166667\n\njulia> nm_measurement(state, VertexSet([[1, 3], [2]]))\n2-element Array{Float64,1}:\n 0.3333333333333333\n 0.6666666666666666\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.nm_loc_ham",
    "page": "Demoralization",
    "title": "QSWalk.nm_loc_ham",
    "category": "function",
    "text": "nm_loc_ham(vertexset[, hamiltoniansByDegree])\n\nReturn Hamiltonian acting locally on each vertex from vertexset linear subspace. hamiltoniansByDegree is a dictionary Dict{Int, SparseDenseMatrix}, which, for a given dimension of vertex linear subspace, yields a hermitian operator. Only matrices for existing dimensions needs to be specified.\n\nNote: Value of vertexset should be generated by make_vertex_set in order to match demoralization procedure. Numerical analysis suggests, that hamiltonians should be complex valued.\n\nExamples\n\njulia> vset = VertexSet([[1, 2], [3, 4]])\nVertexSet(Vertex[Vertex([1, 2]), Vertex([3, 4])])\n\njulia> Matrix(nm_loc_ham(vset))\n4×4 Array{Complex{Float64},2}:\n 0.0+0.0im  0.0+1.0im  0.0+0.0im  0.0+0.0im\n 0.0-1.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+1.0im\n 0.0+0.0im  0.0+0.0im  0.0-1.0im  0.0+0.0im\n\njulia> A = [1 1im; -1im 1]\n2×2 Array{Complex{Int64},2}:\n 1+0im  0+1im\n 0-1im  1+0im\n\njulia> nm_loc_ham(vset, Dict{Int,Matrix{ComplexF64}}(2  => A))\n4×4 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 8 stored entries:\n  [1, 1]  =  1.0+0.0im\n  [2, 1]  =  0.0-1.0im\n  [1, 2]  =  0.0+1.0im\n  [2, 2]  =  1.0+0.0im\n  [3, 3]  =  1.0+0.0im\n  [4, 3]  =  0.0-1.0im\n  [3, 4]  =  0.0+1.0im\n  [4, 4]  =  1.0+0.0im\n\n\n\n\n\nnm_loc_ham(vertexset[, hamiltoniansByVertex])\n\nReturn Hamiltonian acting locally on each vertex from vertexset linear subspace. hamiltoniansByVertex is a dictionary Dict{Vertex, SparseDenseMatrix}, which, for a given vertex, yields a hermitian operator of the size equal to the dimension of the vertex subspace.\n\nNote: Value of vertexset should be generated by make_vertex_set in order to match demoralization procedure. Numerical analysis suggests, that hamiltonians should be complex valued.\n\nExamples\n\njulia> vset = VertexSet([[1, 2], [3, 4]])\nVertexSet(Vertex[Vertex([1, 2]), Vertex([3, 4])])\n\njulia> Matrix(nm_loc_ham(vset))\n4×4 Array{Complex{Float64},2}:\n 0.0+0.0im  0.0+1.0im  0.0+0.0im  0.0+0.0im\n 0.0-1.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im\n 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+1.0im\n 0.0+0.0im  0.0+0.0im  0.0-1.0im  0.0+0.0im\n\njulia> A, B = [1 1im; -1im 1], [0 1; 1 0]\n(Complex{Int64}[1+0im 0+1im; 0-1im 1+0im], [0 1; 1 0])\n\njulia> v1, v2 = vlist(vset)\n2-element Array{Vertex,1}:\n Vertex([1, 2])\n Vertex([3, 4])\n\njulia> nm_loc_ham(vset, Dict{Vertex,Matrix{Number}}(v1  => A, v2  => B))\n4×4 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 6 stored entries:\n  [1, 1]  =  1.0+0.0im\n  [2, 1]  =  0.0-1.0im\n  [1, 2]  =  0.0+1.0im\n  [2, 2]  =  1.0+0.0im\n  [4, 3]  =  1.0+0.0im\n  [3, 4]  =  1.0+0.0im\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.nm_init",
    "page": "Demoralization",
    "title": "QSWalk.nm_init",
    "category": "function",
    "text": "nm_init(init_vertices, vertexset)\n\nCreate initial state in the case of the nonmoralizing evolution based on init_vertices of type Vector{Vertex}. The result is a block diagonal matrix, where each block corresponds to vertex from vertexset. The final state represent an uniform probability over nm_measurement.\n\nNote: The function returns sparse matrix with ComplexF64 field type.\n\nExamples\n\njulia> vset = VertexSet([[1], [2, 3, 4], [5, 6, 7], [8, 9]])\nVertexSet(Vertex[Vertex([1]), Vertex([2, 3, 4]), Vertex([5, 6, 7]), Vertex([8, 9])])\n\njulia> nm_init(vset[[1, 3, 4]], vset)\n9×9 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 6 stored entries:\n  [1, 1]  =  0.333333+0.0im\n  [5, 5]  =  0.111111+0.0im\n  [6, 6]  =  0.111111+0.0im\n  [7, 7]  =  0.111111+0.0im\n  [8, 8]  =  0.166667+0.0im\n  [9, 9]  =  0.166667+0.0im\n\n\n\n\n\nnm_init(init_states, vertexset)\n\nCreate initial state in the case of the nonmoralizing evolution based on init_states of type Dict{Vertex, <:AbstractMatrix{<:Number}}. For each given vertex a block from dictionary is used, otherwise zero matrix is chosen. Each matrix from dictionary should be nonnegative and sum of all traces should equal one. The keys of init_vertices should be a vertices from vertexset. Note that matrix from init_states corresponding to vertex v should be of size length(v)×length(v).\n\nNote: The function returns sparse matrix with ComplexF64 field type.\n\nExamples\n\njulia> vset = VertexSet([[1], [2, 3, 4], [5, 6, 7], [8, 9]])\nVertexSet(Vertex[Vertex([1]), Vertex([2, 3, 4]), Vertex([5, 6, 7]), Vertex([8, 9])])\n\njulia> A1, A2, A3 = ones(ComplexF64, 1, 1)/4, [ 1/5+0im 0 1/5; 0 1/10 0 ; 1/5 0 1/5 ], [0.125 -0.125+0im; -0.125 0.125]\n(Complex{Float64}[0.25+0.0im], Complex{Float64}[0.2+0.0im 0.0+0.0im 0.2+0.0im; 0.0+0.0im 0.1+0.0im 0.0+0.0im; 0.2+0.0im 0.0+0.0im 0.2+0.0im], Complex{Float64}[0.125+0.0im -0.125+0.0im; -0.125+0.0im 0.125+0.0im])\n\njulia> nm_init(Dict(vset[1] =>A1, vset[3] =>A2, vset[4] =>A3), vset)\n9×9 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 10 stored entries:\n  [1, 1]  =  0.25+0.0im\n  [5, 5]  =  0.2+0.0im\n  [7, 5]  =  0.2+0.0im\n  [6, 6]  =  0.1+0.0im\n  [5, 7]  =  0.2+0.0im\n  [7, 7]  =  0.2+0.0im\n  [8, 8]  =  0.125+0.0im\n  [9, 8]  =  -0.125+0.0im\n  [8, 9]  =  -0.125+0.0im\n  [9, 9]  =  0.125+0.0im\n\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.default_nm_loc_ham",
    "page": "Demoralization",
    "title": "QSWalk.default_nm_loc_ham",
    "category": "function",
    "text": "default_nm_loc_ham(size)\n\nReturn default local Hamiltonian of size size×size for the demoralization procedure. The Hamiltonian is sparse with nonzero elements on the first upper diagonal (equal to 1im) and lower diagonal (equal to -1im).\n\nNote: This function is used to provide the default argument for nm_loc_ham function.\n\nExamples\n\njulia> QSWalk.default_nm_loc_ham(4)\n4×4 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 6 stored entries:\n  [2, 1]  =  0.0-1.0im\n  [1, 2]  =  0.0+1.0im\n  [3, 2]  =  0.0-1.0im\n  [2, 3]  =  0.0+1.0im\n  [4, 3]  =  0.0-1.0im\n  [3, 4]  =  0.0+1.0im\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.make_vertex_set",
    "page": "Demoralization",
    "title": "QSWalk.make_vertex_set",
    "category": "function",
    "text": "make_vertex_set(A[, epsilon])\n\nCreates object of type VertexSet which represents how vertices are located in matrix. Should be used only in the nondefault use of evolve_generator function. It is always equal to the second element if output of evolve_generator function.\n\nExamples\n\njulia> A = [1 2 3; 0 3. 4.; 0 0 5.]\n3×3 Array{Float64,2}:\n 1.0  2.0  3.0\n 0.0  3.0  4.0\n 0.0  0.0  5.0\n\njulia> vlist(make_vertex_set(A))\n3-element Array{Vertex,1}:\n Vertex([1, 2, 3])\n Vertex([4, 5])\n Vertex([6])\n\njulia> vlist(make_vertex_set(A, epsilon = 2.5))\n3-element Array{Vertex,1}:\n Vertex([1])\n Vertex([2, 3])\n Vertex([4])\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.vlist",
    "page": "Demoralization",
    "title": "QSWalk.vlist",
    "category": "function",
    "text": "vlist(vset)\n\nReturns the list of vertices for given vset of type VertexSet.\n\njulia> vset = VertexSet([[1], [2, 3, 4], [5, 6, 7], [8, 9]])\nVertexSet(Vertex[Vertex([1]), Vertex([2, 3, 4]), Vertex([5, 6, 7]), Vertex([8, 9])])\n\njulia> vlist(vset)\n4-element Array{Vertex,1}:\n Vertex([1])\n Vertex([2, 3, 4])\n Vertex([5, 6, 7])\n Vertex([8, 9])\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.subspace",
    "page": "Demoralization",
    "title": "QSWalk.subspace",
    "category": "function",
    "text": "subspace(v)\n\nReturns the subspace connected to vertex v.\n\njulia> v = Vertex([1,2,3])\nVertex([1, 2, 3])\n\njulia> subspace(v)\n3-element Array{Int64,1}:\n 1\n 2\n 3\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.fourier_matrix",
    "page": "Demoralization",
    "title": "QSWalk.fourier_matrix",
    "category": "function",
    "text": "fourier_matrix(size)\n\nReturns Fourier matrix of size size×size.\n\nExamples\n\njulia> fourier_matrix(1)\n1×1 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 1 stored entry:\n  [1, 1]  =  1.0+0.0im\n\njulia> fourier_matrix(2)\n2×2 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 4 stored entries:\n  [1, 1]  =  1.0+0.0im\n  [2, 1]  =  1.0+0.0im\n  [1, 2]  =  1.0+0.0im\n  [2, 2]  =  -1.0+1.22465e-16im\n\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#QSWalk.vertexsetsize",
    "page": "Demoralization",
    "title": "QSWalk.vertexsetsize",
    "category": "function",
    "text": "vertexsetsize(vertexset)\n\nReturn the dimension of the linearspace corresponding to given vertexset.\n\nExamples\n\njulia> vertexsetsize(VertexSet([[1, 2, 3], [4, 5]]))\n5\n\n\n\n\n\n"
},

{
    "location": "demoralization.html#Full-docs-1",
    "page": "Demoralization",
    "title": "Full docs",
    "category": "section",
    "text": "Vertex\nVertexSet\nnm_measurement\nnm_loc_ham\nnm_init\ndefault_nm_loc_ham\nmake_vertex_set\nvlist\nsubspace\nfourier_matrix\nvertexsetsize"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Bugs-1",
    "page": "Contributing",
    "title": "Bugs",
    "category": "section",
    "text": "In the case, you noticed some bugs, please start with an issue with a minimal working example of not-working code. If Exception is thrown, please provide an exception message as well.  If no Exception is thrown, but the result is wrong, please provide in the issue message correct answer.In the case you make a pull request, please add a not-working example as a test."
},

{
    "location": "contributing.html#Improvements-1",
    "page": "Contributing",
    "title": "Improvements",
    "category": "section",
    "text": "If you can provide a code, which works faster than already existing, please check its efficiency for various input data. We welcome any ideas concerning the readability and logic of the code as well."
},

{
    "location": "contributing.html#Development-guidelines-1",
    "page": "Contributing",
    "title": "Development guidelines",
    "category": "section",
    "text": "Post an issue.\nWait until the discussion ends.\nCreate assertions on argument types and other requirements.\nInclude necessary references.\nWrite tests."
},

{
    "location": "citing.html#",
    "page": "Citing",
    "title": "Citing",
    "category": "page",
    "text": ""
},

{
    "location": "citing.html#Usage-and-citing-1",
    "page": "Citing",
    "title": "Usage and citing",
    "category": "section",
    "text": "In case of citing, please use the following BibTeX form:@article{glos2018qswalk,\n  title={QSWalk. jl: Julia package for quantum stochastic walks analysis},\n  author={Glos, Adam and Miszczak, Jaros{\\l}aw Adam and Ostaszewski, Mateusz},\n  journal={arXiv preprint arXiv:1801.01294},\n  year={2018}\n}Our package was already used in papers concerning quantum attacksAdam Glos, Jarosław Adam Miszczak, and Mateusz Ostaszewski. \'Limiting properties of stochastic quantum walks on directed graphs.\' Journal of Physics A: Mathematical and Theoretical 51.3 (2017): 035304.\nKrzysztof Domino, Adam Glos, and Mateusz Ostaszewski. \'Superdiffusive quantum stochastic walk definable on arbitrary directed graph\'. Quantum Information & Computation, 17(11-12), 973-986.In case You have used our package for your research, we will be grateful for any information about the paper or the package. With your consent we will provide a link to the paper here."
},

{
    "location": "license.html#",
    "page": "Licence",
    "title": "Licence",
    "category": "page",
    "text": "MIT LicenseCopyright (c) 2017-2018 Adam Glos, Jarosław Adam Miszczak, Mateusz OstaszewskiPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

]}
