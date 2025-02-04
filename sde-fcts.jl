using Parameters, StochasticDiffEq, DataFrames, CSV, UnicodePlots

# deterministic component
function f!(du, u, p, t)

    # remove negative abundances (due to numerical round-off errors)
    ninds = findall(x->x<0,u)
    u[ninds] .= 0.0

    K = floor(Int64,p[1])
    r = p[2:(K+1)]
    α = p[(K+2):(2*K+1)]
    s = p[(2*K+2):(3*K+1)]
    Ψ = reshape(last(p,K*K),K,K)

    # deterministic differentials
    d = (r - α .* u) .* u  +  s .* Ψ * u

    for k in 1:K
        du[k] = d[k]
    end

end

# stochastic component
function g!(du, u, p, t)

    # remove negative abundances (due to numerical round-off errors)
    ninds = findall(x->x<0,u)
    u[ninds] .= 0.0

    K = floor(Int64,p[1])
    β = p[(3*K+2):(4*K+1)]

    d = β .* .√u

    for k in 1:K
        du[k] = d[k]
    end

end
