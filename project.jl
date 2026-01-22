using Random
using Distributions
using Statistics
using Plots

# Taguchi loss: L(y) = k (y - m)^2
taguchi_loss(y::Real, m::Real, k::Real) = k * (y - m)^2

struct QueueSimResult
    arrival_times::Vector{Float64}
    service_times::Vector{Float64}
    wait_times::Vector{Float64}
    idle_times::Vector{Float64}
    loss_wait::Vector{Float64}
    loss_idle::Vector{Float64}
end

"""
Simulate a single-server queue.

Note: This keeps the SAME arrival-time construction used in the report:
arrival_times = cumsum(Poisson(λ) draws).
"""
function simulate_queue_report_model(;
    λ::Real = 10,
    μ::Real = 7,
    σ::Real = 2,
    trials::Int = 1000,
    optimal::Real = 2.0,
    k_wait::Real = 10,
    k_idle::Real = 5,
    seed::Int = 42
)::Tuple{QueueSimResult, NamedTuple}

    Random.seed!(seed)

    # Match report approach (integer jumps, then cumulative sum)
    arrival_times = Float64.(cumsum(rand(Poisson(λ), trials)))
    service_times = abs.(rand(Normal(μ, σ), trials))  # ensure nonnegative, as in report

    wait_times = zeros(Float64, trials)
    idle_times = zeros(Float64, trials)
    loss_wait  = zeros(Float64, trials)
    loss_idle  = zeros(Float64, trials)

    server_busy_until = 0.0

    @inbounds for i in 1:trials
        a = arrival_times[i]
        s = service_times[i]

        # If server is busy, customer waits; otherwise server is idle until next arrival
        w = max(server_busy_until - a, 0.0)
        idl = max(a - server_busy_until, 0.0)

        wait_times[i] = w
        idle_times[i] = idl

        # Correct timing update: service starts when both the customer has arrived AND the server is free
        start_service = max(a, server_busy_until)
        server_busy_until = start_service + s

        loss_wait[i] = taguchi_loss(w, optimal, k_wait)
        loss_idle[i] = taguchi_loss(idl, optimal, k_idle)
    end

    result = QueueSimResult(arrival_times, service_times, wait_times, idle_times, loss_wait, loss_idle)

    summary = (
        mean_wait = mean(wait_times),
        sd_wait   = std(wait_times),
        mean_idle = mean(idle_times),
        sd_idle   = std(idle_times),
        mean_loss_wait = mean(loss_wait),
        mean_loss_idle = mean(loss_idle),
    )

    return result, summary
end

# Run it
res, stats = simulate_queue_report_model()
println(stats)

# Optional: histograms like your Figure 1
p1 = histogram(res.wait_times, title="Histogram of Wait Times", xlabel="Wait time (minutes)", ylabel="Frequency")
p2 = histogram(res.idle_times, title="Histogram of Idle Times", xlabel="Idle time (minutes)", ylabel="Frequency")
plot(p1, p2, layout=(1,2))
