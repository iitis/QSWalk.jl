include("qswalks_types.jl")

export
  integrate_lind,
  evolve

function integrate_lind(H::SparseMatrixCSC{Complex128,Int64}, L::Vector{SparseMatrixCSC{Complex128,Int64}},w::Float64)
    m = size(L,1)
    n = size(H,1)
    F = spzeros(n^2,n^2)
    id = sparse(eye(H))
    for i = 1:m
        F += kron(L[i],conj(L[i]))-0.5*kron(L[i]'*L[i],id)-0.5*kron(id,transpose(L[i])*conj(L[i]))
    end
    F = w*F + (1-w)*im*(kron(id,conj(H))-kron(H,id))
    F
end

function evolve(F::SparseMatrixCSC{Complex128,Int64}, init::Array{Complex{Float64},1}, t::Float64)
    W = expmv( t, F, init)
    A = real(diag(unres(W)))
    round(A,4)
end
