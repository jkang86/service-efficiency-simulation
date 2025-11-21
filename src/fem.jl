module FEMSolver

export fem_solve

using LinearAlgebra

"""
    fem_solve(c, s, f, a, b, n)

Solves linear BVP:
    -(c(x)u')' + s(x)u = f(x),  u(a)=u(b)=0
Using linear FEM on n subintervals.

Returns (x, u).
"""
function fem_solve(c, s, f, a, b, n)
    h = (b - a)/n
    x = collect(a:h:b)

    # element matrices (from your PDF)
    Ke = [1.0 -1.0; -1.0 1.0]
    Me = (1/6) * [2.0 1.0; 1.0 2.0]
    fe = (1/2) * [1.0; 1.0]

    cbar = (c.(x[1:end-1]) .+ c.(x[2:end])) ./ 2
    sbar = (s.(x[1:end-1]) .+ s.(x[2:end])) ./ 2
    fbar = (f.(x[1:end-1]) .+ f.(x[2:end])) ./ 2

    K = zeros(n-1, n-1)
    M = zeros(n-1, n-1)
    rhs = zeros(n-1)

    # boundary contributions
    K[1,1] = cbar[1]/h
    M[1,1] = sbar[1]*h/3
    rhs[1] = fbar[1]*h/2

    K[end,end] = cbar[end]/h
    M[end,end] = sbar[end]*h/3
    rhs[end] = fbar[end]*h/2

    # interior assembly
    for k in 2:n-1
        K[k-1:k, k-1:k] .+= (cbar[k]/h) * Ke
        M[k-1:k, k-1:k] .+= (sbar[k]*h) * Me
        rhs[k-1:k] .+= (fbar[k]*h) * fe
    end

    u_int = (K + M) \ rhs
    u = vcat(0.0, u_int, 0.0)

    return x, u
end

end # module
