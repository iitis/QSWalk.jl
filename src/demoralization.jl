export
  localhamiltonian,
  demoralizedlindbladian

"""

    partitionsize(partition)

# Examples

```jldoctest
julia>
```
"""
function partitionsize(partition::Vector{Vector{Int}})
  sum([length(block) for block=partition])
end

"""

    defaultlocalhamiltonian(size)

# Examples

```jldoctest
julia>
```
"""
function defaultlocalhamiltonian(size::Int)
  if size == 1
    return spzeros(Complex128,1,1)
  else
    spdiagm((im*ones(size-1),-im*ones(size-1)),(1,-1))
  end
end


"""

    localhamiltonian(partition[, hamiltonians, mode])

# Arguments
-
-
-

# Return

# Examples

```jldoctest
julia>
```
"""
function localhamiltonian(partition::Vector{Vector{Int}},
      hamiltonians::Function=x->defaultlocalhamiltonian(x),
      mode::String="size")

  result = spzeros(Complex128,partitionsize(partition),partitionsize(partition))
  if mode == "size"
    for block=partition
      result[block,block] = hamiltonians(length(block))
    end
  elseif mode == "index"
    for (index,block)=enumerate(partition)
      result[block,block] = hamiltonians(index)
    end
  else
    error(ArgumentError("mode should be size or index"))
  end
  result
end

"""

    reversedincidencelist(A[, ϵ])


# Examples

```jldoctest
julia>
```
"""
function reversedincidencelist{T<:FieldType}(A::SparseMatrixCSC{T}, ϵ::Real=eps(1.))
  [(A[:,i].nzind)[find(x -> abs(x)>=ϵ, A[:,i].nzind)] for i=1:size(A,1)]
end

function reversedincidencelist{T<:FieldType}(A::Matrix{T}, ϵ::Real=eps(1.))
  [find(x -> abs(x)>=ϵ, A[:,i]) for i=1:size(A,1)]
end

"""

    incidencelist(A[, ϵ])


# Examples

```jldoctest
julia>
```
"""
function incidencelist(A::SparseMatrixCSC{FieldType}, ϵ::Real=eps(1.))
  [(A[:,i].nzind)[find(x -> abs(x)>=ϵ, A[i,:].nzind)] for i=1:size(A,1)]
end

function incidencelist(A::Matrix{FieldType}, ϵ::Real=eps(1.))
  [find(x -> abs(x)>=ϵ, A[i,:]) for i=1:size(A,1)]
end

"""

    makepartition(revincidencelist)


# Examples

```jldoctest
julia>
```
"""
function makepartition(revincidencelist::Vector{Vector{Int}})
  result = Vector{Int}[]
  start = 1
  for i=revincidencelist
    if length(i)!=0
      push!(result, collect(start:(start+length(i)-1)))
      start+=length(i)
    else
      push!(result, [start])
      start += 1
    end
  end
  result
end

"""

    rectangularfouriermatrix(size, columnnumber)

# Examples

```jldoctest
julia>
```
"""
function rectangularfouriermatrix(size::Int,column::Int)
  result = zeros(Complex128,size)
  result = [ exp(2im*π*(i-1)*(column-1)/size) for i=1:size]
  result
end

"""

    demoralizedlindbladian(A, ϵ)

# Arguments
 -
 -
 -

# Return


# Examples

```jldoctest
julia>
```
"""
function demoralizedlindbladian(A::SparseDenseMatrix, ϵ::Real=eps(1.),
  lindbladians::Function=(x,y)->rectangularfouriermatrix(x,y), mode::String="size")
  revincidencelist = reversedincidencelist(A, ϵ)
  partition = makepartition(revincidencelist)

  L = spzeros(typeof(A[1,1]),partitionsize(partition),partitionsize(partition))
  if mode == "size"
    for i=1:size(A,1), (index,j)=enumerate(revincidencelist[i]), k in partition[j]
        L[partition[i],k] = A[i,j]*lindbladians(length(partition[i]), index)
    end
  elseif mode == "index"
    for i=1:size(A,1), (index,j)=enumerate(revincidencelist[i]), k in partition[j]
        #TODO test, may be wrong
        L[partition[i],k] = A[i,j]*lindbladians(i,j)
    end
  else
    error(ArgumentError("mode should be size or index"))
  end

  L, partition
end
