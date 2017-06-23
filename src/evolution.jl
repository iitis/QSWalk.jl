include("utils.jl")

export
  integrate_lind,
  evolve



function evolve(F::SparseMatrixCSC{Complex128,Int64}, init::Array{Complex{Float64},1}, t::Float64)
    W = expmv( t, F, init)
    A = real(diag(unres(W)))
    round(A,4)
end
