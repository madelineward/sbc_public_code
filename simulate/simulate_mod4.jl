###### SIMULATE FROM HILL-TYPE ALARM


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
        df, gt_temp = sim_data_mod4a(seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, [αᵤ βᵤ δ₁ δ₂])
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end


##### Model 4A (Behaviour change in alpha)

δ₁ = 3
δ₂ = 0.05

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4a_a24b2d3d05.txt", unif_pops)


δ₁ = 3
δ₂ = 0.075

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4a_a24b2d3d075.txt", unif_pops)



δ₁ = 3
δ₂ = 0.1

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4a_a24b2d3d10.txt", unif_pops)






function gen_unif_pops() # This function generates n=num_pop epidemics across spatially uniform populations
    unif_pops = DataFrame(x = Float64[], y = Float64[], status = Int64[], inftime = Int64[], rem_time = Int64[])
    gt = []
    for i in 1:num_pop
        x, y = spatial_uniform(num_ind, rangex, rangey, seeds[i])
        df, gt_temp = sim_data_mod5b(αᵤ, βᵤ, seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, δ₁, δ₂)
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end




##### Model 4B (Behaviour change in beta)

δ₁ = 3
δ₂ = 0.10

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4b_a24b2d3d10.txt", unif_pops)


δ₁ = 3
δ₂ = 0.15

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4b_a24b2d3d15.txt", unif_pops)


δ₁ = 3
δ₂ = 0.20

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))


CSV.write("data/mod4b_a24b2d3d20.txt", unif_pops)



