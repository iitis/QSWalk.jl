include("qswalks_types.jl")

export
  res,
  unres

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
function res{T<:FieldType}(::Type{T}, mtx::Array{T,2} )
  Base.vec(transpose(mtx))
end

function res(mtx::Array{Complex128,2})
  res(Complex128,mtx)
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

function unres{T<:FieldType}(::Type{T}, v::Array{T,1}, dims::Tuple{IndexType,IndexType})
  reshape(v,dims)
end

function unres{T<:FieldType}(::Type, v::Array{T,1})
  d = floor(Int64,sqrt(length(v)))
  if d*d == length(v)
    transpose(reshape(v,(d,d)))
  else
    throw(ArgumentError("Expected vector with perfect square number of elements."))
  end
end

unres(v::Array{Complex128,1}) = unres(Complex128,v)
