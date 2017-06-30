export
  ket,
  bra,
  ketbra,
  proj,
  reshuffle,
  unreshuffle

typealias SparseDenseMatrix{T<:Number} Union{SparseMatrixCSC{T},Matrix{T}}
typealias SparseDenseVector{T<:Number} Union{SparseVector{T},Vector{T}}

type Vertex
  linearspace::Vector{Int}
end

macro argument(ex, msgs...)
    local msg_body = isempty(msgs) ? ex : msgs[1]
    local msg = string(msg_body)
    return :($ex ? nothing : throw(ArgumentError($msg)))
end

"""
    ket([type, ]index, size)

Return `index`-th base (column) vector in the `size`-dimensional vector space.
To be consistent with Julia indexing, `index`=1,2,...,`size`. The `type` defaults
to `Complex128` if not specified

# Examples

```jldoctest
julia> ket(1,2)
2-element Array{Complex{Float64},1}:
 1.0+0.0im
 0.0+0.0im

julia> ket(Float64,1,2)
2-element Array{Float64,1}:
 1.0
 0.0
```
"""
function ket{T<:Number}(::Type{T}, index::Int, size::Int)
  @argument size > 0 "vector size must be positive"
  @assert 1 <= index <= size "index must be greater than 0 and lower than vector size"

  ret = zeros(T, size)
  ret[index] = 1
  ret
end

ket(index::Int, size::Int) = ket(Complex128, index, size)

"""
    bra([type, ]index, size)

Return `index`-th row vector in the `size`-dimensional vector space, with
`index`=1,2,...,`size`. The `type` defaults to `Complex128` if not specified

# Examples

```jldoctest
julia> bra(1,2)
1×2 Array{Complex{Float64},2}:
 1.0-0.0im  0.0-0.0im

 julia> bra(Float64,1,2)
 1×2 Array{Float64,2}:
  1.0  0.0
```
"""
function bra{T<:Number}(::Type{T}, index::Int, size::Int)
  ket(T, index, size)'
end

bra(index::Int, size::Int) = ket(index, size)'

"""
    ketbra([type,] indexrow, indexcolumn, size)

Return matrix acting on `size`-dimensional vector space,
`indexrow`,`indexcolumn`=1,2,...,`size`. The matrix consists of single nonzero
element equal to one at position (`indexrow`,`indexcolumn`). The `type` defaults to `Complex128` if not
specified.

# Examples

```jldoctest
julia> ketbra(1,2,3)
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im

 julia> ketbra(Float64,2,1,3)
 3×3 Array{Float64,2}:
  0.0  0.0  0.0
  1.0  0.0  0.0
  0.0  0.0  0.0
```
"""
function ketbra{T<:Number}(::Type{T},
                              indexrow::Int,
                              indexcolumn::Int,
                              size::Int)
  ket(T, indexrow, size)*bra(T, indexcolumn, size)
end

function ketbra(indexrow::Int, indexcolumn::Int, size::Int)
  ketbra(Complex128, indexrow, indexcolumn, size)
end

"""
    proj([type,] index, size)

Return projector onto `index`-th base vector in `size`-dimensional vector space,
with `index`=1,2,...,`size`. The `type` defaults to `Complex128` if not
specified.


    proj(vector)

Return projector onto the subspace spanned by vector `vector`.

# Examples
```jldoctest
julia> proj(1,2)
2×2 Array{Complex{Float64},2}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> proj(Float64,1,2)
2×2 Array{Float64,2}:
 1.0  0.0
 0.0  0.0

julia> v = 1/sqrt(2) * (ket(1,3)+ket(3,3))
3-element Array{Complex{Float64},1}:
 0.707107+0.0im
      0.0+0.0im
 0.707107+0.0im

julia> QSW.proj(v)
3×3 Array{Complex{Float64},2}:
 0.5+0.0im  0.0+0.0im  0.5+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.5+0.0im  0.0+0.0im  0.5+0.0im
```
"""
function proj{T<:Number}(::Type{T}, index::Int, size::Int)
  ketbra(T, index, index, size)
end

proj(index::Int, size::Int) = proj(Complex128, index, size)

function proj{T<:Number}(vector::SparseDenseVector{T})
  vector*vector'
end

"""

    reshuffle(matrix)

Return vectorization of the `matrix` in the row order. This is equivalent to
`Base.vec(transpose(matrix)`.

# Examples

```jldoctest
julia> M = Matrix{Float64}(reshape(1:9, (3,3))')
3×3 Array{Float64,2}:
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0

julia> reshuffle(M)
9-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 5.0
 6.0
 7.0
 8.0
 9.0

julia> reshuffle(unreshuffle(v)) == v
true
```

"""
function reshuffle{T<:Number}(matrix::SparseDenseMatrix{T})
  Base.vec(transpose(matrix))
end



"""

    unreshuffle(vector)

Return square matrix elements from `vector`. The `vector`
is expected to have perfect square number of arguments to form square matrix.

# Examples
```jldoctest

julia> v = Vector{Float64}(collect(1:9))
9-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 5.0
 6.0
 7.0
 8.0
 9.0

julia> unreshuffle(v)
3×3 Array{Float64,2}:
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0

julia> reshuffle(unreshuffle(v)) == v
true
```

"""

function unreshuffle{T<:Number}(vector::SparseDenseVector{T})
  dim = floor(Int64,sqrt(length(vector)))
  @argument dim*dim == length(vector) "Expected vector with perfect square number of elements."

  reshape(vector, (dim,dim))'
end
