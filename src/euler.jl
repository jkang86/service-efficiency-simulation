module EulerSolver

export euler_method

"""
    euler_method(f, a, b, u0, h)

Euler's Method for IVP:
    u' = f(t, u),  u(a) = u0,  t∈[a,b]

Returns (t, u).
"""
function euler_method(f, a, b, u0, h)
    n = Int((b - a) / h)
    t = range(a, step=h, length=n+1)
    u = zeros(n+1)
    u[1] = u0

    for i in 1:n
        u[i+1] = u[i] + h * f(t[i], u[i])
    end

    return collect(t), u
end

end # module
