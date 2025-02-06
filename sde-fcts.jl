using Parameters, StochasticDiffEq, LinearAlgebra, CSV, DataFrames, Plots, PlotThemes

"""
# model parameters

- T = maximum time to run simulation to
- K = number of hosts
- S = number of microbial taxa
- Ψ = social transmission network
- Ξ = environmental transmission network
- r = intrinsic growth rates
- α = competition coefficients
- s = social transmission
- e = environmental transmission
- β = magnitudes of demographic stochasticity
- N₀ = initial abundances

"""
@with_kw mutable struct PARS

    # SYNTAX

    # the `::` assigns a type to the variable
    # the `=` assigns a default value

    # so `T::Float64 = 10.0` defines a 
    # continuous value `T` taking `10.0` as a default

    # maximum time to run simulation to
    T::Float64 = 10.0

    # number of hosts
    K::Int64 = 4

    # number of microbial taxa
    S::Int64 = 2

    # social transmission network
    Ψ::Matrix{Float64} = Symmetric(abs.(randn(K,K))).*(ones(K,K)-I)

    # environmental transmission network
    Ξ::Matrix{Float64} = Symmetric(abs.(randn(K,K))).*(ones(K,K)-I)
    
    # intrinsic growth rates
    # r[i,k] is growth rate of microbe i in host k
    r::Matrix{Float64} = 10*ones(S,K)

    # competition coefficients
    # α[k][i,j] is coeff between microbes i and j in host k
    α::Vector{Matrix{Float64}} = 0.1*fill(abs.(randn(S,S)),K)

    # social transmission
    s::Vector{Float64} = 0.001*ones(S)

    # environmental transmission
    e::Vector{Float64} = 0.001*ones(S)

    # magnitudes of demographic stochasticity
    β::Matrix{Float64} = ones(S,K)

    # initial abundances
    N₀::Matrix{Float64} = 100*ones(S,K)

end

# deterministic component
function f(N, p, t)

    # # remove negative abundances (due to numerical round-off errors)
    ninds = findall(x->x<0,N)
    N[ninds] .= 0.0

    # unpack parameters
    @unpack_PARS p

    # turns the competition bit into a matrix
    αN::Matrix{Float64} = zeros(S,K)
    for k in 1:K
        αN[:,k] = α[k]*N[:,k] # the * here does matrix multiplication
    end

    # turns transmission bits into matrices
    sN::Matrix{Float64} = zeros(S,K)
    eN::Matrix{Float64} = zeros(S,K)
    for k in 1:K
        sN[:,k] = s .* N[:,k] # the .* does element-wise multiplication
        eN[:,k] = e .* N[:,k] # the .* does element-wise multiplication
    end

    # effects of within-host dynamics
    WTHN_HOST = (r - αN) .* N # the .* does element-wise multiplication

    # effects of between-host dynamics (ie, transmission)
    BTWN_HOST = sN*Ψ + sN*Ξ

    # net deterministic dynamics
    dN = WTHN_HOST + BTWN_HOST

    return dN

end

# stochastic component
function g(N, p, t)

    # remove negative abundances (due to numerical round-off errors)
    ninds = findall(x->x<0,N)
    N[ninds] .= 0.0

    # unpack model parameter
    @unpack_PARS p
 
    dN = β .* .√N

    return dN

end

function runSim(p::PARS)
    @unpack_PARS p

    # define timespan to simulate over
    tspan = (0.0,T)

    # run simulation
    prob = SDEProblem(f, g, N₀, tspan, p)
    out = solve(prob, SRIW1())

    return out

end

# saves parameter values as csv with filename fn
function savePARS(p::PARS;rbN0="rbN0",alpha="alpha",se="se",TKS="TKS",social="social",enviro="enviro")
    @unpack_PARS p

    # each column is a vectorized matrix of parameters
    # column lengths are S*K
    dfSK = DataFrame(r = vec(r), β = vec(β), N₀ = vec(N₀))

    # each column is a vectorized matrix of comp coeffs for host k
    # column lengths are S*S
    dfSS = DataFrame(α1 = vec(α[1]))
    for k in 2:K
        nm = "α"*string(k)
        dfSS[!,nm] = vec(α[k])
    end

    # each column is vector of microbe transmission "propensities"
    # so each column is length S
    dfS = DataFrame(s=s, e=e)

    # simulation time T, number of hosts K, and number of microbes s
    # ... each column length 1
    df1 = DataFrame(T=T,K=K,S=S)

    # social transmission matrix
    dfΨ = DataFrame(PSIx1 = Ψ[:,1])
    for k in 2:K
        nm = "PSIx"*string(k)
        dfΨ[!,nm] = Ψ[:,k]
    end

    # environmental transmission matrix
    dfΞ = DataFrame(XIx1 = Ξ[:,1])
    for k in 2:K
        nm = "XIx"*string(k)
        dfΞ[!,nm] = Ξ[:,k]
    end

    CSV.write(rbN0*".csv", dfSK)
    CSV.write(alpha*".csv", dfSS)
    CSV.write(se*".csv", dfS)
    CSV.write(TKS*".csv", df1)
    CSV.write(social*".csv", dfΨ)
    CSV.write(enviro*".csv", dfΞ)


end

# saves time series as csv with filename fn
function saveTS(out;abunds="abunds",times="times",pts=200)

    n = length(out.t)

    if pts > n
        println("\nMore time points requested than available
                 \nDefaulting to maximum time points available: pts=",n)
        pts = n
    end

    inds = round.(Int,LinRange(1,n,pts))

    # time points
    # single column of length pts
    dft = DataFrame(t=out.t[inds])

    # microbe abundances
    # each column is the S*K abundances at a time pt
    dfN = DataFrame(N_t1 = vec(out.u[1]))
    for k in 2:pts
        nm = "N_t"*string(k)
        dfN[!,nm] = vec(out.u[k])
    end
    
    CSV.write(times*".csv", dft)
    CSV.write(abunds*".csv", dfN)
    
end
    