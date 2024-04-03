## CODE TO FIT MODEL 1A WHEN GENERATIVE MODEL WAS DIFFERENT (MISSPECIFIED ALARM)

# Run tasks go from 4-24

using Distributed
Nchains = 3
addprocs(Nchains)

# Load packages on all workers
@everywhere begin
    using Distributions
    using LinearAlgebra
    using DataFrames
    using CSV
    using StatsBase
    using DelimitedFiles
    using Random
    using MCMCDiagnosticTools
    using MCMCChains
end

# Load supporting code on all workers
@everywhere begin
    include("source/misc.jl")
    include("source/mcmc.jl")
    include("source/models.jl")
    pop_task_str = ENV["SLURM_ARRAY_TASK_ID"]
    pop_task = parse(Int, pop_task_str)
    run_task_str = ENV["A"]
    run_task = parse(Int, run_task_str)
end

println("Run task = ", run_task_str)

@everywhere begin
    tmax = 30
    tmin_f = 9
    tmax_f = 20
    num_ind = 1000
    num_pop = 1
    λ = 3
    n_iter = 100000
    n_burn = Int(0.1*n_iter)
    n_hpd = 500
end

########### FITTING MODEL 1A TO MODEL 1B DATA ##############

if run_task == 4
    # Weak BC
    @everywhere begin
        mod_gen = "mod1b"
        mod_fit = "mod1a"
        scen = "a22b2d10h40"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end


if run_task == 5
    @everywhere begin
       # Medium BC
        mod_gen = "mod1b"
        mod_fit = "mod1a"
        scen = "a22b2d15h40"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 6
    @everywhere begin
        # Strong BC
        mod_gen = "mod1b"
        mod_fit = "mod1a"
        scen = "a22b2d20h40"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end


########### MODEL 2A ##############

if run_task == 7
    @everywhere begin
        # Weak BC
        mod_gen = "mod2a"
        mod_fit = "mod1a"
        scen = "a24b2d005"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 8
    @everywhere begin
        # Medium BC
    mod_gen = "mod2a"
    mod_fit = "mod1a"
    scen = "a24b2d01"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 9
    @everywhere begin
        # Strong BC
    mod_gen = "mod2a"
    mod_fit = "mod1a"
    scen = "a24b2d015"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

########### MODEL 2B ##############

if run_task == 10
    @everywhere begin
        # Weak BC
    mod_gen = "mod2b"
    mod_fit = "mod1a"
    scen = "a24b2d001"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 11
    @everywhere begin
        # Medium BC
    mod_gen = "mod2b"
    mod_fit = "mod1a"
    scen = "a24b2d0015"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 12
    @everywhere begin
        # Strong BC
        mod_gen = "mod2b"
        mod_fit = "mod1a"
        scen = "a24b2d002"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

########### MODEL 3A ##############

if run_task == 13
    @everywhere begin
        # Weak BC
        mod_gen = "mod3a"
        mod_fit = "mod1a"
        scen = "a24b2d80d02"
        init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end


if run_task == 14
    @everywhere begin
        # Medium BC
    mod_gen = "mod3a"
    mod_fit = "mod1a"
    scen = "a24b2d80d03"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 15
    @everywhere begin
        # Strong BC
    mod_gen = "mod3a"
    mod_fit = "mod1a"
    scen = "a24b2d80d04"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

########### MODEL 3B ##############

if run_task == 16
    @everywhere begin
       # Weak BC
    mod_gen = "mod3b"
    mod_fit = "mod1a"
    scen = "a24b2d04d005"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 17
    @everywhere begin
       # Medium BC
    mod_gen = "mod3b"
    mod_fit = "mod1a"
    scen = "a24b2d04d007"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 18
    @everywhere begin
       # Strong BC
    mod_gen = "mod3b"
    mod_fit = "mod1a"
    scen = "a24b2d04d009"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        # Note the order of d1 and d2 is swapped in paper
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

########### MODEL 4A ##############

if run_task == 19
    @everywhere begin
       # Weak BC
    mod_gen = "mod4a"
    mod_fit = "mod1a"
    scen = "a24b2d3d05"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 20
    @everywhere begin
       # Medium BC
    mod_gen = "mod4a"
    mod_fit = "mod1a"
    scen = "a24b2d3d075"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 21
    @everywhere begin
      # Strong BC
    mod_gen = "mod4a"
    mod_fit = "mod1a"
    scen = "a24b2d3d10"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

########### MODEL 4B ##############

if run_task == 22
    @everywhere begin
      # Weak BC
    mod_gen = "mod4b"
    mod_fit = "mod1a"
    scen = "a24b2d3d10"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 23
    @everywhere begin
      # Medium BC
    mod_gen = "mod4b"
    mod_fit = "mod1a"
    scen = "a24b2d3d15"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end

if run_task == 24
    @everywhere begin
      # Strong BC
    mod_gen = "mod4b"
    mod_fit = "mod1a"
    scen = "a24b2d3d20"
    init_vals = [2.5 1.8 0.60 50]
        num_param = length(init_vals)

        epi_data = DataFrame(CSV.File("data/$(lpad(mod_gen,1))_$(lpad(scen,1)).txt"; delim = ','))
        popi = filter(:pop => ==(pop_task), epi_data)

        # Define prior
        function LogPrior(p) 
            α = p[1]
            β = p[2]
            δ₁ = p[3]
            δ₂ = p[4]
            prior_α = logpdf(Uniform(0,100), α)
            prior_β = logpdf(Uniform(0,100), β)
            prior_δ₁ = logpdf(Beta(1,1), δ₁)
            prior_δ₂ = logpdf(Gamma(3,20), δ₂)
            return prior_α + prior_β + prior_δ₁ + prior_δ₂
        end
        
        # Input population data into likelihood function
        function lik_pop(param2, tmax2) 
            likelihood_mod1a(popi.inftime, popi.rem_time, tmax2, param2, popi.x, popi.y, num_ind)
        end

        ## FUNCTIONS FOR FULL TIME PERIOD

        # Likelihood function for t=1-30
        function lik_pop_full(param2)
            lik_pop(param2, tmax)
        end
        
        function get_chain_full()
            RWMH_MCMC(lik_pop_full, LogPrior, num_param, init_vals, n_iter, tmax)
        end

        ## FUNCTIONS FOR FORECASTING

        # Likelihood function for t=1-8
        function lik_pop_forecast(param2)
            lik_pop(param2, tmin_f-1)
        end
        
        function get_chain_forecast()
            RWMH_MCMC(lik_pop_forecast, LogPrior, num_param, init_vals, n_iter, tmin_f-1)
        end

    end

    # RUN MCMC ON FULL STUDY PERIOD
    chains1, chains2, chains3 = pmap(x -> get_chain_full(),1:Nchains)

    process_chains(chains1, chains2, chains3, lik_pop_full, sim_data_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmax)


    # RUN MCMC ON FIRST 8 TIME POINTS THEN FORECAST
    chains1, chains2, chains3 = pmap(x -> get_chain_forecast(),1:Nchains)

    process_forecasts(chains1, chains2, chains3, lik_pop_forecast, sim_data_mod1a, forecast_mod1a, n_iter, n_burn, n_hpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)

end