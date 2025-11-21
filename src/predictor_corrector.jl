module PredictorCorrectorSolver

export predictor_corrector, trapezoid_step, adams_bashforth4, adams_moulton3

"""
Trapezoid initialization step (as in your PDF).
"""
function trapezoid_step(f, t0, y0, h)
    return y0 + (h/2) * f(t0, y0)
end

"""
Adams–Bashforth 4-step predictor.
f_vals = [f_{n-3}, f_{n-2}, f_{n-1}, f_n]
"""
function adams_bashforth4(y_n, f_vals, h)
    return y_n + (h/24) * (
        55*f_vals[4] - 59*f_vals[3] + 37*f_vals[2] - 9*f_vals[1]
    )
end

"""
Adams–Moulton 3-step corrector.
Uses predicted f_{n+1}.
"""
function adams_moulton3(y_n, f_vals, h, f_next)
    return y_n + (h/24) * (
        9*f_next + 19*f_vals[4] - 5*f_vals[3] + f_vals[2]
    )
end

"""
    predictor_corrector(f, t0, y0, h, n_steps)

Matches your PDF algorithm exactly:
- First 3 steps via trapezoid initialization
- Then AB4 predictor + AM3 corrector
"""
function predictor_corrector(f, t0, y0, h, n_steps)
    t = collect(t0:h:(t0 + n_steps*h))
    y = zeros(length(t))
    y[1] = y0

    # store last 4 f-values
    f_vals = zeros(4)

    # initialization for y2,y3,y4
    for i in 1:3
        y[i+1] = trapezoid_step(f, t[i], y[i], h)
        f_vals[i] = f(t[i], y[i])
    end
    f_vals[4] = f(t[4], y[4])

    for i in 4:(length(t)-1)
        y_pred = adams_bashforth4(y[i], f_vals, h)
        f_pred = f(t[i+1], y_pred)
        y[i+1] = adams_moulton3(y[i], f_vals, h, f_pred)

        # shift f-values forward
        f_vals = vcat(f_vals[2:4], f(t[i+1], y[i+1]))
    end

    return t, y
end

end # module

