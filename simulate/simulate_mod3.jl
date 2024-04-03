###### SIMULATE FROM SCALED EXPONENTIAL ALARM

###### SET-UP #######

# Load packages
using Random
using Distributions
using LinearAlgebra
using DataFrames
using CSV
using StatsBase
using Plots
using Plots.PlotMeasures

# External functions code inclusion
include("source/misc.jl")
include("source/models.jl")

# Constants
tmax = 50
rangex = [100 300]
rangey = [100 300]
num_ind = 1000
num_pop = 20
λ = 3
Random.seed!(1234)
seeds = sample(1000:10000, Dims([num_pop]), replace = false)
#init_inf = sample(1:num_ind, Dims([num_pop, 3]); replace=false)

αᵤ = 2.4 # Alpha for all spatially uniform populations
βᵤ = 2 # Beta for all spatially uniform populations


# Functions
function gen_unif_pops() # This function generates n=num_pop epidemics across spatially uniform populations
    unif_pops = DataFrame(x = Float64[], y = Float64[], status = Int64[], inftime = Int64[], rem_time = Int64[])
    gt = []
    for i in 1:num_pop
        x, y = spatial_uniform(num_ind, rangex, rangey, seeds[i])
        df, gt_temp = sim_data_mod3a(seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, [αᵤ βᵤ δ₁ δ₂])
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end

#### MODEL 3A

δ₁ = 0.8
δ₂ = 0.02

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))
CSV.write("data/mod3a_a24b2d80d02.txt", unif_pops)


δ₁ = 0.8
δ₂ = 0.03

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[3001:4000,:].inftime,tmax))
CSV.write("data/mod3a_a24b2d80d03.txt", unif_pops)


δ₁ = 0.8
δ₂ = 0.04

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[2001:3000,:].inftime,tmax))
CSV.write("data/mod3a_a24b2d80d04.txt", unif_pops)


function gen_unif_pops() # This function generates n=num_pop epidemics across spatially uniform populations
    unif_pops = DataFrame(x = Float64[], y = Float64[], status = Int64[], inftime = Int64[], rem_time = Int64[])
    gt = []
    for i in 1:num_pop
        x, y = spatial_uniform(num_ind, rangex, rangey, seeds[i])
        df, gt_temp = sim_data_mod3b(seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, [αᵤ βᵤ δ₁ δ₂])
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end

##### MODEL 3B

δ₁ = 0.4
δ₂ = 0.005

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))
CSV.write("data/mod3b_a24b2d04d005.txt", unif_pops)

δ₁ = 0.4
δ₂ = 0.007

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))
CSV.write("data/mod3b_a24b2d04d007.txt", unif_pops)


δ₁ = 0.4
δ₂ = 0.009

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))
CSV.write("data/mod3b_a24b2d04d009.txt", unif_pops)


