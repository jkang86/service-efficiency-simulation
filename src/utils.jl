module Utils

export l2_error

using LinearAlgebra

"""
    l2_error(u_num, u_exact, x)

L2 norm error between numerical and exact solutions.
"""
function l2_error(u_num, u_exact, x)
    return norm(u_exact.(x) - u_num, 2)
end

end # module
