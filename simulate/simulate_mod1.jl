#### SIMULATE FROM BC MODEL WITH THRESHOLD ALARM


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

αᵤ = 2.2 # Alpha for all spatially uniform populations
βᵤ = 2 # Beta for all spatially uniform populations

# Functions
function gen_unif_pops() # This function generates n=num_pop epidemics across spatially uniform populations
    unif_pops = DataFrame(x = Float64[], y = Float64[], status = Int64[], inftime = Int64[], rem_time = Int64[])
    gt = []
    for i in 1:num_pop
        x, y = spatial_uniform(num_ind, rangex, rangey, seeds[i])
        df, gt_temp = sim_data_mod1a(seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, [αᵤ βᵤ δ H])
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end


##### Model 1A (Behaviour change in alpha)

# With α = 2.2 and β = 2

δ = 0
H = 0

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1a_a22b2d0h0.txt", unif_pops)

δ = 0.50
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1a_a22b2d50h40.txt", unif_pops)


δ = 0.65
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1a_a22b2d65h40.txt", unif_pops)


δ = 0.80
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1a_a22b2d80h40.txt", unif_pops)







function gen_unif_pops() # This function generates n=num_pop epidemics across spatially uniform populations
    unif_pops = DataFrame(x = Float64[], y = Float64[], status = Int64[], inftime = Int64[], rem_time = Int64[])
    gt = []
    for i in 1:num_pop
        x, y = spatial_uniform(num_ind, rangex, rangey, seeds[i])
        df, gt_temp = sim_data_mod1b(seeds[i], num_ind, x, y, λ, init_inf[i,:], tmax, [αᵤ βᵤ δ H])
        unif_pops = vcat(unif_pops, df)
        gt = vcat(gt, gt_temp)
    end
    insertcols!(unif_pops, 1, :pop => repeat(1:num_pop, inner = num_ind))
    return unif_pops, gt
end

##### Model 1B (Behaviour change in beta)

δ = 0.10
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1b_a22b2d10h40.txt", unif_pops)


δ = 0.15
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1b_a22b2d15h40.txt", unif_pops)



δ = 0.20
H = 40

unif_pops, gt = gen_unif_pops()
Plots.plot(count_inf(unif_pops[1:1000,:].inftime,tmax))

CSV.write("data/mod1b_a22b2d20h40.txt", unif_pops)




# Draw random values from priors and simulate curves

function draw_prior(n) 
    α = rand(Uniform(0,100),n)
    β = rand(Uniform(0,100),n)
    δ = rand(Beta(1,1), n)
    H = rand(Gamma(3,20),n)
    return [α β δ H]
end

tmax = 30
vals = draw_prior(100)
curves_t = zeros(100, tmax+1)
x, y  = spatial_uniform(num_ind, rangex, rangey, 123)

Random.seed!(1234)
seeds = sample(1000:10000, 100, replace = false)
for i in 1:100
    param = vals[i,:]
    unif_pops, gt = sim_data_mod1a(seeds[i], num_ind, x, y, λ, init_inf[1,:], tmax, param)
    curves_t[i,:] = count_inf(unif_pops.inftime,tmax)
end

p1 = Plots.plot(curves_t[1,:], legend = false)
Plots.xlabel!(p1, "Time")
Plots.ylabel!(p1, "Incident Infections")
Plots.title!(p1, "Model 1A")
for i in 2:100
    Plots.plot!(p1, curves_t[i,:], legend = false)
end


function draw_prior(n) 
    α = rand(Uniform(0,100),n)
    β = rand(Uniform(0,100),n)
    δ = rand(Beta(1,1), n)
    return [α β δ]
end

tmax = 30
vals = draw_prior(100)
curves_t2 = zeros(100, tmax+1)
x, y  = spatial_uniform(num_ind, rangex, rangey, 123)

Random.seed!(1234)
seeds = sample(1000:10000, 100, replace = false)
for i in 1:100
    param = vals[i,:]
    unif_pops, gt = sim_data_mod2a(seeds[i], num_ind, x, y, λ, init_inf[1,:], tmax, param)
    curves_t2[i,:] = count_inf(unif_pops.inftime,tmax)
end

p2 = Plots.plot(curves_t2[1,:], legend = false)
Plots.xlabel!(p2, "Time")
Plots.ylabel!(p2, "Incident Infections")
Plots.title!(p2, "Model 2A")
for i in 2:100
    Plots.plot!(p2, curves_t2[i,:], legend = false)
end



function draw_prior(n) 
    α = rand(Uniform(0,100),n)
    β = rand(Uniform(0,100),n)
    δ₁ = rand(Beta(1,1), n)
    δ₂ = rand(Beta(1,2), n)
    return [α β δ₁ δ₂]
end

tmax = 30
vals = draw_prior(100)
curves_t3 = zeros(100, tmax+1)
x, y  = spatial_uniform(num_ind, rangex, rangey, 123)

Random.seed!(1234)
seeds = sample(1000:10000, 100, replace = false)
for i in 1:100
    param = vals[i,:]
    unif_pops, gt = sim_data_mod3a(seeds[i], num_ind, x, y, λ, init_inf[1,:], tmax, param)
    curves_t3[i,:] = count_inf(unif_pops.inftime,tmax)
end

p3 = Plots.plot(curves_t3[1,:], legend = false)
Plots.xlabel!(p3, "Time")
Plots.ylabel!(p3, "Incident Infections")
Plots.title!(p3, "Model 3A")
for i in 2:100
    Plots.plot!(p3, curves_t3[i,:], legend = false)
end


function draw_prior(n) 
    α = rand(Uniform(0,100),n)
    β = rand(Uniform(0,100),n)
    δ₁ = rand(Gamma(2,4), n)
    δ₂ = rand(Beta(1,2), n)
    return [α β δ₁ δ₂]
end

tmax = 30
vals = draw_prior(100)
curves_t4 = zeros(100, tmax+1)
x, y  = spatial_uniform(num_ind, rangex, rangey, 123)

Random.seed!(1234)
seeds = sample(1000:10000, 100, replace = false)
for i in 1:100
    param = vals[i,:]
    unif_pops, gt = sim_data_mod4a(seeds[i], num_ind, x, y, λ, init_inf[1,:], tmax, param)
    curves_t4[i,:] = count_inf(unif_pops.inftime,tmax)
end

p4 = Plots.plot(curves_t4[1,:], legend = false)
Plots.xlabel!(p4, "Time")
Plots.ylabel!(p4, "Incident Infections")
Plots.title!(p4, "Model 4A")
for i in 2:100
    Plots.plot!(p4, curves_t4[i,:], legend = false)
end

savefig(Plots.plot(p1,p2,p3,p4, layout = 4, size = (1600, 1200), guidefont = (18), tickfont = (14), titlefont = (18), bottom_margin = 20mm, top_margin = 20mm, left_margin = 10mm, right_margin = 10mm), "epi_curves_under_priors.pdf")
