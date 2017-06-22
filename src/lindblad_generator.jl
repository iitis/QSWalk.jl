
export
  gener_hamiltonian_from_lindblad,
  gener_class_lind_from_lindblad


function gener_hamiltonian_from_lindblad(L::SparseMatrixCSC{Complex128,Int64})
    H = zeros(L)
    n = size(H,1)
    for i=1:n,j=1:n
        if L[i,j] != 0 || L[j,i] != 0
            H[i,j] = 1
        end
    end
    H
end

#TODO: test
function gener_class_lind_from_lindblad(LL::SparseMatrixCSC{Complex128,Int64})
    L = SparseMatrixCSC{Complex128,Int64}[]
    n = size(LL,1)
    for i=1:n,j=1:n
        if LL[i,j] != 0
            push!(L, sparse(LL[i,j]*ketbra(i,j,n)))
        end
    end
    L
end

