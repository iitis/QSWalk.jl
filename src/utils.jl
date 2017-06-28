export
  ket,
  bra,
  ketbra,
  proj,
  res,
  unres


typealias FieldType Union{Real,Complex}
typealias SparseDenseMatrix{T<:FieldType} Union{SparseMatrixCSC{T},Matrix{T}}
typealias SparseDenseVector{T<:FieldType} Union{SparseVector{T},Vector{T}}


"""
    ket(i,n)

Return i-th base (column) vector in the n-dimensional vector space. To be consistent with
Julia indexing, i=1,2,...,n.

# Examples

```jldoctest
julia> ket(1,2)
2-element Array{Complex{Float64},1}:
 1.0+0.0im
 0.0+0.0im

julia> ket(Complex64,1,2)
2-element Array{Complex{Float32},1}:
 1.0+0.0im
 0.0+0.0im

```
"""
function ket{T<:FieldType}(::Type{T}, i::Int, n::Int)
  ret = zeros(T,n)
  ret[i] = 1
  ret
end

ket(i::Int, n::Int) = ket(Complex128,i,n)

"""
    bra(i,n)

Return i-th row vector in the n-dimensional vector space, with i=1,2,...,n. This vector
is defied as a complex conjugate of the column vector ket(i,n).

# Examples

```jldoctest
julia> bra(2,2)
1×2 Array{Complex{Float64},2}:
 0.0-0.0im  1.0-0.0im

julia> bra(Complex64,2,2)
1×2 Array{Complex{Float32},2}:
 0.0-0.0im  1.0-0.0im
```
"""
function bra{T<:FieldType}(::Type{T}, i::Int, n::Int)
  ket(T,i,n)'
end

bra(i::Int, n::Int) = ket(Complex128,i,n)'

"""
    ketbra(i,j,n)

Return matrix |i><j| acting on n-dimensional vector space, i,j=1,2,...n. By default this
function returns matrix of type Array{Complex{Float64},2}.

If necessary it can be instructed to return matrix of type Array{Complex{Float32},2}.

# Exmaples

```jldoctest
julia> ketbra(1,2,3)
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im

julia> ketbra(Complex64,2,1,3)
3×3 Array{Complex{Float32},2}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im

```
"""
function ketbra{T<:FieldType}(::Type{T}, i::Int, j::Int, n::Int)
  ket(T,i,n)*bra(T,j,n)
end

function ketbra(i::Int,j::Int,n::Int)
  ketbra(Complex128,i,j,n)
end

"""
    proj(i,n)

Return projector onto i-th base vector in n-dimensional vector space.

    proj(v)

Return projector onto the subspace spanned by vector v.

# Examples
```jldoctest
julia> proj(1,2)
2×2 Array{Complex{Float64},2}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> proj(Complex64,1,2)
2×2 Array{Complex{Float32},2}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
```

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

"""
function proj{T<:FieldType}(::Type{T}, i::Int, n::Int)
  ketbra(T,i,i,n)
end

proj(i::Int,n::Int) = proj(Complex128,i,n)

function proj(v::SparseDenseVector)
  v*v'
end



"""

    res(mtx)

Return vectorization of the matrix mtx in the row order. This is equivalent to
Base.vec(transpose(mtx).

# Examples

```jldoctest
julia> v = ketbra(1,2,4)
4×4 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

julia> res(v)
16-element Array{Complex{Float64},1}:
 0.0+0.0im
 1.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
 0.0+0.0im
```

"""
function res(matrix::SparseDenseMatrix)
  Base.vec(transpose(matrix))
end



"""

    unres(vec,[dims])

Return matrix with dimension in dims and elements from vec. If the second argument is
omitted, the vector is expected to have perfect square number of arguments to form
square matrix.

# Examples
```jldoctest

julia> unres(ket(7,9))
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im
 1.0+0.0im  0.0+0.0im  0.0+0.0im

julia> res(unres(ket(7,9))) == ket(7,9)
true

```

"""

function unres(v::SparseDenseVector, dims::Tuple{Int,Int})
  reshape(v,dims)
end

function unres(::Type, v::SparseDenseVector)
  d = floor(Int64,sqrt(length(v)))
  if d*d == length(v)
    transpose(reshape(v,(d,d)))
  else
    throw(ArgumentError("Expected vector with perfect square number of elements."))
  end
end

unres(v::SparseDenseVector) = unres(Complex128,v)
