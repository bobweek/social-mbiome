using Distributions, DataFrames, CSV

#
# generate simulation parameters (modulo Ψ) CSV file
#

# random number of hosts
K = 4

# random growth rates across hosts
# r = rand(Exponential(1),K)
r = [1.0,2.0,3.0,4.0]

# random intraspecific competition coefficients across hosts
# α = rand(Exponential(0.01),K)
α = 0.01*ones(K)

# random propensity to transmit (should be small)
s = rand(Exponential(0.01))

# random levels of demographic stochasticity across hosts
#   (perhaps this can be constant across hosts?)
#β = rand(Exponential(1),K)
β = ones(K)

# random initial abundances across hosts
#N₀ = rand(Exponential(30),K)
N₀ = 200*ones(K)

mdf = DataFrame(r = r, α = α, s = s.*ones(K), β = β, N₀ = N₀)

CSV.write("model-parameters.csv", mdf)

#
# generate social network (Ψ) as edge matrix and writes to CSV file
#

# can make it binary
Ψ = rand(Bernoulli(0.9),K,K)

# or continuous
#Ψ = rand(Exponential(1),K,K)

# output is handled same either way
CSV.write("network.csv", Tables.table(Ψ), writeheader=false)

#
# generate simulation parameters
#

# duration of simulation
T = 10 #rand(Exponential(100))

# number of replicates
R = 1 + rand(Poisson(5))

# make dataframe
sdf = DataFrame(T=T, R=R)

# write as csv
CSV.write("sim-parameters.csv", sdf)
