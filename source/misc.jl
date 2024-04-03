###### MISCELLANEOUS FUNCTIONS #######

# Miscellaneous functions 

# function to count number of newly infectious individuals at each timepoint
# returns vector of length tmax + 1, with number of individuals entering infectious period at t=1,2,...tmax+1
function count_inf(inftimes, tmax) 
    x = collect(1:(tmax+1))           # Time points of the epidemic
    counts = repeat([0], outer = [length(x)])         # Placeholder vector of correct length
    for i in 1:length(x)              # For each timepoint, find the number of new infections
        counts[i] = sum(inftimes .== x[i])  
    end
    return counts
end

function count_inf2(inftimes, tmin, tmax) 
    x = collect(tmin:(tmax+1))           # Time points of the epidemic
    counts = repeat([0], outer = [length(x)])         # Placeholder vector of correct length
    for i in 1:length(x)              # For each timepoint, find the number of new infections
        counts[i] = sum(inftimes .== x[i])  
    end
    return counts
end

# function to calculate WAIC
function WAIC(x)
    ppd    = mapslices(mean, exp.(x); dims = 1)
    lppd   = sum(log.(ppd))
    pWAIC2 = sum(mapslices(var, x; dims = 1))
    WAIC   = -2*lppd + 2*pWAIC2
    return WAIC
end

# function to calculate DIC
function DIC(loglike, LLpostmean)
    Dbar    = -2 * mean(loglike)
    Dhat    = -2 * LLpostmean
    pD      = Dbar - Dhat
    dic     = pD + Dbar
    return dic
end

# In case seed is not consistent across machines, these are the initial infections for uniform pops:
init_inf = [860  840  245;
119  224  844;
478  688  125;
161  206  967;
652  926  756;
332  969   12;
519  901  930;
535  646   80;
299  846  446;
213  894  796;
298  854  892;
469  600   18;
 84  314  456;
802  434  843;
441  454  311;
569  244  308;
315  281  155;
326  811  784;
 97  640  187;
440  482  680]

# function that will generate spatial coordinates randomly generated from Uniform distributions
function spatial_uniform(num_ind, rangex, rangey, seed)
    Random.seed!(seed)
    N = num_ind
    xmin = rangex[1]; xmax = rangex[2]
    x = rand(Uniform(xmin, xmax), N)
    ymin = rangey[1]; ymax = rangey[2]
    y = rand(Uniform(ymin, ymax), N)
    return x, y
end


# Remove burn-in, calculate Gelman-Rubin statistic, obtain parameter estimates and CIs, calculate WAIC, calculate curve HPDs for full data set
function process_chains(chain1, chain2, chain3, likelihood, sim_fun, niter, burnin, nhpd, mod_gen, mod_fit, scen, pop_task, tmax)
    num_params = size(chain1[1])[2]
    param_names = [:alpha, :beta, :delta1, :delta2]
    all_chains = cat(chain1[1], chain2[1], chain3[1], dims = 3)

    # Remove burn-in iterations
    no_burnin = mapslices(x -> last(x, niter - burnin), all_chains; dims = 1)
    # Put all chains together for hpds etc.
    no_burnin_long = vcat(no_burnin[:,:,1], no_burnin[:,:,2], no_burnin[:,:,3])

    no_burnin_chobj = Chains(no_burnin, param_names[1:num_params])

    # Gelman-Rubin Statistic - converged if psrfci < 1.2
    r = DataFrame(gelmandiag(no_burnin_chobj))
    r.mod_gen = repeat([mod_gen],num_params)
    r.mod_fit = repeat([mod_fit],num_params)
    r.scen = repeat([scen],num_params)
    r.pop = repeat([pop_task],num_params)
    CSV.write("output/gr_full.txt", r; append = true)

    # Posterior parameter means and HPDs
    param_means = DataFrame(mean(no_burnin_chobj))
    param_hpds = DataFrame(hpd(no_burnin_chobj))
    param_means.lower = param_hpds.lower
    param_means.upper = param_hpds.upper
    param_means.mod_gen = repeat([mod_gen],num_params)
    param_means.mod_fit = repeat([mod_fit],num_params)
    param_means.scen = repeat([scen],num_params)
    param_means.pop = repeat([pop_task],num_params)

    #println(param_means)
    
    # Save parameter summaries to text file
    CSV.write("output/params_full.txt", param_means; append = true)

    # Format loglikelihoods and remove burnin for WAIC
    all_loglik = cat(chain1[2], chain2[2], chain3[2], dims = 3)
    all_loglik = mapslices(x -> last(x, niter - burnin), all_loglik; dims = 1)
    all_loglik = vcat(all_loglik[:,:,1], all_loglik[:,:,2], all_loglik[:,:,3])

    # Calculate WAIC and DIC (initial columns of loglik are for each timepoint, final column of loglik is overall)
    waic = WAIC(all_loglik[:, 1:(end-1)])
    
    post_means = mapslices(mean, no_burnin_long; dims = 1)
    LLpostmeans = likelihood(post_means)
    dic = DIC(all_loglik[:,end], LLpostmeans[end])

    # Put all metrics in same data frame
    metrics = DataFrame(WAIC = waic, DIC = dic, mod_gen = mod_gen, mod_fit = mod_fit, scen = scen, pop = pop_task)
    # Write metrics to text file
    CSV.write("output/waic_and_dic_full.txt", metrics; append = true)
    

    ### CALCULATE EPIDEMIC CURVE HPDS ###
    # Get indices
    ind_hpd = sample(1:Int(niter-burnin), nhpd; replace = false)
    # Initialize matrix to hold curves and alarms
    hpd_mat_c = zeros(nhpd,tmax+1)
    hpd_mat_a = zeros(nhpd,tmax)
    for i in 1:nhpd
        ind = ind_hpd[i]
        p = no_burnin_long[ind,:]
        d, g = sim_fun(sample(1000:10000,1)[1], 1000, popi.x, popi.y, 3, init_inf[pop_task,:], tmax, p)
        hold = count_inf(d.inftime, tmax)
        hpd_mat_c[i,:] = hold
        hpd_mat_a[i,:] = g
    end
    #println(hpd_mat_c)
    curve_med = mapslices(median, hpd_mat_c; dims = 1)
    #println(curve_med)
    curve_lower = mapslices(x -> quantile(x,0.025), hpd_mat_c; dims = 1)
    curve_upper = mapslices(x -> quantile(x,0.975), hpd_mat_c; dims = 1)

    alarm_med = mapslices(median, hpd_mat_a; dims = 1)
    alarm_lower = mapslices(x -> quantile(x,0.025), hpd_mat_a; dims = 1)
    alarm_upper = mapslices(x -> quantile(x,0.975), hpd_mat_a; dims = 1)

    curve_df = DataFrame(hcat(collect(1:(tmax+1)), transpose(curve_med), transpose(curve_lower), transpose(curve_upper), repeat([pop_task],tmax+1), repeat([mod_gen],tmax+1), repeat([scen],tmax+1), repeat([mod_fit],tmax+1)),:auto)
    CSV.write("output/curves_full.txt", curve_df; append = true)

    alarm_df = DataFrame(hcat(collect(1:tmax), transpose(alarm_med), transpose(alarm_lower), transpose(alarm_upper), repeat([pop_task],tmax), repeat([mod_gen],tmax), repeat([scen],tmax), repeat([mod_fit],tmax)),:auto)
    CSV.write("output/alarms_full.txt", alarm_df; append = true)

end

# Remove burn-in, calculate Gelman-Rubin statistic, obtain parameter estimates and CIs, calculate WAIC, calculate curve HPDs for forecasting
function process_forecasts(chain1, chain2, chain3, likelihood, sim_fun, forecast_fun, niter, burnin, nhpd, mod_gen, mod_fit, scen, pop_task, tmin_f, tmax_f)
    tmax = tmin_f - 1
    num_params = size(chain1[1])[2]
    param_names = [:alpha, :beta, :delta1, :delta2]
    all_chains = cat(chain1[1], chain2[1], chain3[1], dims = 3)

    # Remove burn-in iterations
    no_burnin = mapslices(x -> last(x, niter - burnin), all_chains; dims = 1)
    # Put all chains together for hpds etc.
    no_burnin_long = vcat(no_burnin[:,:,1], no_burnin[:,:,2], no_burnin[:,:,3])

    no_burnin_chobj = Chains(no_burnin, param_names[1:num_params])

    # Gelman-Rubin Statistic - converged if psrfci < 1.2
    r = DataFrame(gelmandiag(no_burnin_chobj))
    r.mod_gen = repeat([mod_gen],num_params)
    r.mod_fit = repeat([mod_fit],num_params)
    r.scen = repeat([scen],num_params)
    r.pop = repeat([pop_task],num_params)
    CSV.write("output/gr_forecast.txt", r; append = true)

    # Posterior parameter means and HPDs
    param_means = DataFrame(mean(no_burnin_chobj))
    param_hpds = DataFrame(hpd(no_burnin_chobj))
    param_means.lower = param_hpds.lower
    param_means.upper = param_hpds.upper
    param_means.mod_gen = repeat([mod_gen],num_params)
    param_means.mod_fit = repeat([mod_fit],num_params)
    param_means.scen = repeat([scen],num_params)
    param_means.pop = repeat([pop_task],num_params)

    #println(param_means)
    
    # Save parameter summaries to text file
    CSV.write("output/params_forecast.txt", param_means; append = true)

    # Format loglikelihoods and remove burnin for WAIC
    all_loglik = cat(chain1[2], chain2[2], chain3[2], dims = 3)
    all_loglik = mapslices(x -> last(x, niter - burnin), all_loglik; dims = 1)
    all_loglik = vcat(all_loglik[:,:,1], all_loglik[:,:,2], all_loglik[:,:,3])

    # Calculate WAIC and DIC (initial columns of loglik are for each timepoint, final column of loglik is overall)
    waic = WAIC(all_loglik[:, 1:(end-1)])
    
    post_means = mapslices(mean, no_burnin_long; dims = 1)
    LLpostmeans = likelihood(post_means)
    dic = DIC(all_loglik[:,end], LLpostmeans[end])

    # Put all metrics in same data frame
    metrics = DataFrame(WAIC = waic, DIC = dic, mod_gen = mod_gen, mod_fit = mod_fit, scen = scen, pop = pop_task)
    # Write metrics to text file
    CSV.write("output/waic_and_dic_forecast.txt", metrics; append = true)
    

    ### CALCULATE EPIDEMIC CURVE HPDS ###
    # Get indices
    ind_hpd = sample(1:Int(niter-burnin), nhpd; replace = false)
    # Initialize matrix to hold curves and alarms
    hpd_mat_c = zeros(nhpd,tmax+1)
    hpd_mat_a = zeros(nhpd,tmax)
    for i in 1:nhpd
        ind = ind_hpd[i]
        p = no_burnin_long[ind,:]
        d, g = sim_fun(sample(1000:10000,1)[1], 1000, popi.x, popi.y, 3, init_inf[pop_task,:], tmax, p)
        hold = count_inf(d.inftime, tmax)
        hpd_mat_c[i,:] = hold
        hpd_mat_a[i,:] = g
    end
    #println(hpd_mat_c)
    curve_med = mapslices(median, hpd_mat_c; dims = 1)
    #println(curve_med)
    curve_lower = mapslices(x -> quantile(x,0.025), hpd_mat_c; dims = 1)
    curve_upper = mapslices(x -> quantile(x,0.975), hpd_mat_c; dims = 1)

    alarm_med = mapslices(median, hpd_mat_a; dims = 1)
    alarm_lower = mapslices(x -> quantile(x,0.025), hpd_mat_a; dims = 1)
    alarm_upper = mapslices(x -> quantile(x,0.975), hpd_mat_a; dims = 1)

    curve_df = DataFrame(hcat(collect(1:(tmax+1)), transpose(curve_med), transpose(curve_lower), transpose(curve_upper), repeat([pop_task],tmax+1), repeat([mod_gen],tmax+1), repeat([scen],tmax+1), repeat([mod_fit],tmax+1)),:auto)
    CSV.write("output/curves_forecast.txt", curve_df; append = true)

    alarm_df = DataFrame(hcat(collect(1:tmax), transpose(alarm_med), transpose(alarm_lower), transpose(alarm_upper), repeat([pop_task],tmax), repeat([mod_gen],tmax), repeat([scen],tmax), repeat([mod_fit],tmax)),:auto)
    CSV.write("output/alarms_forecast.txt", alarm_df; append = true)

    init_inf2 = zeros(num_ind)
    for i in 1:num_ind
        if popi.inftime[i] < tmin_f
            init_inf2[i] = popi.inftime[i]
        end
    end

    ind_hpd = sample(Int(n_burn+1):n_iter, n_hpd; replace = false)
    hpd_mat = zeros(n_hpd,tmax_f+2-tmin_f)
    g_mat = zeros(n_hpd,tmax_f+1-tmin_f)
    for i in 1:n_hpd
        ind = ind_hpd[i]
        p = no_burnin_long[ind,:]
        #println(first(p, 6))
        d, g = forecast_fun(sample(1000:10000,1)[1], 1000, popi.x, popi.y, 3, init_inf2, tmin_f, tmax_f, p)
        g_mat[i,:] = g[tmin_f:end]
        hold = count_inf2(d.inftime, tmin_f, tmax_f)
        hpd_mat[i,:] = hold
    end

    #println(hpd_mat_c)
    curve_med = mapslices(median, hpd_mat; dims = 1)
    #println(curve_med)
    curve_lower = mapslices(x -> quantile(x,0.025), hpd_mat; dims = 1)
    curve_upper = mapslices(x -> quantile(x,0.975), hpd_mat; dims = 1)

    alarm_med = mapslices(median, g_mat; dims = 1)
    alarm_lower = mapslices(x -> quantile(x,0.025), g_mat; dims = 1)
    alarm_upper = mapslices(x -> quantile(x,0.975), g_mat; dims = 1)

    curve_df = DataFrame(hcat(collect(tmin_f:(tmax_f+1)), transpose(curve_med), transpose(curve_lower), transpose(curve_upper), repeat([pop_task],tmax_f+2-tmin_f), repeat([mod_gen],tmax_f+2-tmin_f), repeat([scen],tmax_f+2-tmin_f), repeat([mod_fit],tmax_f+2-tmin_f)),:auto)
    CSV.write("output/forecasted_curves.txt", curve_df; append = true)

    alarm_df = DataFrame(hcat(collect(tmin_f:tmax_f), transpose(alarm_med), transpose(alarm_lower), transpose(alarm_upper), repeat([pop_task],tmax_f+1-tmin_f), repeat([mod_gen],tmax_f+1-tmin_f), repeat([scen],tmax_f+1-tmin_f), repeat([mod_fit],tmax_f+1-tmin_f)),:auto)
    CSV.write("output/forecasted_alarms.txt", alarm_df; append = true)

end