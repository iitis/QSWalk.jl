export
  localhamiltonian,
  demoralizedlindbladian

type FunctionByIndex
  f::Function
end

function (f::FunctionByIndex)(x::Int)
  f.Function(x)
end

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
-f

# Return

# Examples

```jldoctest
julia>
```
"""
function localhamiltonian(partition::Vector{Vector{Int}},
                          hamiltonians::Function=defaultlocalhamiltonian)

  result = spzeros(Complex128,partitionsize(partition),partitionsize(partition))
  for block=partition
    result[block,block] = hamiltonians(length(block))
  end
  result
end

function localhamiltonian(partition::Vector{Vector{Int}},
                          hamiltonians::FunctionByIndex)
  result = spzeros(Complex128,partitionsize(partition),partitionsize(partition))
  for (index,block)=enumerate(partition)
    result[block,block] = hamiltonians(index)
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
function reversedincidencelist{T<:Number}(A::SparseMatrixCSC{T}, ϵ::Real=eps(1.))
  [(A[:,i].nzind)[find(x -> abs(x)>=ϵ, A[:,i].nzind)] for i=1:size(A,1)]
end

function reversedincidencelist{T<:Number}(A::Matrix{T}, ϵ::Real=eps(1.))
  [find(x -> abs(x)>=ϵ, A[:,i]) for i=1:size(A,1)]
end

"""

    incidencelist(A[, ϵ])


# Examples

```jldoctest
julia>
```
"""
function incidencelist(A::SparseMatrixCSC{Number}, ϵ::Real=eps(1.))
  [(A[:,i].nzind)[find(x -> abs(x)>=ϵ, A[i,:].nzind)] for i=1:size(A,1)]
end

function incidencelist(A::Matrix{Number}, ϵ::Real=eps(1.))
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
  lindbladians::Function=rectangularfouriermatrix)
  revincidencelist = reversedincidencelist(A, ϵ)
  partition = makepartition(revincidencelist)

  L = spzeros(typeof(A[1,1]),partitionsize(partition),partitionsize(partition))
  for i=1:size(A,1), (index,j)=enumerate(revincidencelist[i]), k in partition[j]
      L[partition[i],k] = A[i,j]*lindbladians(length(partition[i]), index)
  end

  L, partition
end

function demoralizedlindbladian(A::SparseDenseMatrix,
  lindbladians::FunctionByIndex, ϵ::Real=eps(1.))
  revincidencelist = reversedincidencelist(A, ϵ)
  partition = makepartition(revincidencelist)

  L = spzeros(typeof(A[1,1]),partitionsize(partition),partitionsize(partition))
  for i=1:size(A,1), (index,j)=enumerate(revincidencelist[i]), k in partition[j]
      L[partition[i],k] = A[i,j]*lindbladians.f(length(partition[i]), index)
  end

  L, partition
end
