facts("evolution") do
    context("Exponent function") do
    @fact evolve_operator(zeros(4, 4), 1.) --> roughly(eye(4))

    # type check
    @fact_throws MethodError evolve_operator(spzeros(4, 4), 1.)
  end

  context("Evolution functions") do
    #size 4x4
    A = Complex{Float64}[0.354177+0.0im 0.0891553-0.0251879im 0.0702961+0.0516828im 0.0708664+0.0767941im; 0.0891553+0.0251879im 0.336055+0.0im 0.0420202-0.0109173im 0.0683605-0.00692846im; 0.0702961-0.0516828im 0.0420202+0.0109173im 0.212401+0.0im 0.0939615+0.0553555im; 0.0708664-0.0767941im 0.0683605+0.00692846im 0.0939615-0.0553555im 0.0973671+0.0im]
    #trivial evolutions
    @fact evolve(zeros(16, 16), A, 0.) --> roughly(A)
    @fact evolve(zeros(16, 16), sparse(A), 0.) --> roughly(A)
    @fact evolve(zeros(16, 16), A, [0., 5., 10.]) --> [A, A, A]

    @fact evolve(rand(16, 16), A, 0.) --> roughly(A)

    @fact evolve(spzeros(16, 16), A, 0.) --> roughly(A)
    @fact evolve(spzeros(16, 16), A, [0., 5., 10.])[1] --> roughly(A)
    @fact evolve(spzeros(16, 16), A, [0., 5., 10.])[2] --> roughly(A)
    @fact evolve(spzeros(16, 16), A, [0., 5., 10.])[3] --> roughly(A)

    @fact evolve(sparse(rand(16, 16)), A, 0.) --> roughly(A)
  end
end
