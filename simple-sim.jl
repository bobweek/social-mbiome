include("sde-fcts.jl")

# model parameters
m = CSV.read("model-parameters.csv",DataFrame)

# can read binary or continuous valued edge matrix
Ψ = CSV.File("network.csv", header=false) |> Tables.matrix

# simulation parameters
s = CSV.read("sim-parameters.csv",DataFrame)

# p is the "vector" of parameters (where each entry can be any kind of object)
p = Vector{Float64}(undef,1)
p[1] = length(m.r) # = K, the number of hosts
append!(p,m.r)
append!(p,m.α)
append!(p,m.s)
append!(p,m.β)
append!(p,vec(Ψ))

tspan = (0.0, s.T[1])

# numerically solve SDE
prob = SDEProblem(f!, g!, m.N₀, tspan, p);
sol = solve(prob, SRIW1(), progress=true)

N = reduce(hcat,sol.u)'

pl = lineplot(sol.t,N,xlabel="time",ylabel="N",ylim=(0,maximum(N)),name=["N1","N2","N3","N4"])

print(pl)

#
# for reference see: https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/
#
