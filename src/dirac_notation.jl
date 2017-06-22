include("qswalks_types.jl")

export
  ket,
  bra,
  ketbra,
  proj

"""
    ket(i,n)

Return i-th base (column) vector in the n-dimensional vector space. To be consistent with
Julia indexing, i=1,2,...,n.

Elements of the vector can be specified as Complex64 or Complex128.

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
function ket{T<:FieldType}(::Type{T}, i::IndexType, n::IndexType)
  ret = zeros(T,n)
  ret[i] = 1
  ret
end

ket(i::IndexType, n::IndexType) = ket(Complex128,i,n)

"""
    bra(i,n)

Return i-th row vector in the n-dimensional vector space, with i=1,2,...,n. This vector
is defied as a complex conjugate of the column vector ket(i,n).

Elements of the vector can be specifies to as Complex64 or Complex128.

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
function bra{T<:FieldType}(::Type{T}, i::IndexType, n::IndexType)
  ket(T,i,n)'
end

bra(i::IndexType, n::IndexType) = ket(Complex128,i,n)'

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
function ketbra{T<:FieldType}(::Type{T}, i::IndexType, j::IndexType, n::IndexType)
  ket(T,i,n)*bra(T,j,n)
end

function ketbra(i::IndexType,j::IndexType,n::IndexType) 
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

julia> v = 1/sqrt(2)*(ket(1,3)+ket(3,3))
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
function proj{T<:FieldType}(::Type{T}, i::IndexType, n::IndexType)
  ketbra(T,i,i,n)
end

proj(i::IndexType,n::IndexType) = ketbra(Complex128,i,i,n)

function proj{T<:FieldType}(v::Array{T,1})
  v*v'
end
