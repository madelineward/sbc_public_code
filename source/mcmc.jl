####### MCMC FUNCTIONS #######

# Adaptive Metropolis MCMC algorithm of Roberts and Rosenthal 2006 (Roberts and Rosenthal Random Walk Metropolis Hastings)
# likelihood should be a function that takes only parameter values and returns log-likelihood for each time point to allow WAIC calculation
# LogPrior should be a function that takes parameter values and returns summed log priors
# Algorithm will be adaptive for first end_adapt iterations
# Adaptive Metropolis MCMC algorithm of Roberts and Rosenthal 2006/2008
# likelihood should be a function that takes only parameter values and returns log-likelihood for each time point to allow WAIC calculation
# LogPrior should be a function that takes parameter values and returns log priors added up

function RWMH_MCMC(likelihood, LogPrior, num_param, init_guess, n_iter, tmax)
    sigmaMat1 = (Matrix(I, num_param, num_param)*(0.1)^2)/num_param
    waicMat = zeros(n_iter, (tmax+1))
    loglik = zeros(n_iter, 1)
    waicMat[1,:] = likelihood(init_guess)
    loglik[1] = waicMat[1,(tmax+1)]
    x = zeros(Float64, n_iter, num_param)
    x[1, :] = init_guess
    for i in 2:n_iter
        old_p = x[i-1,:]
        if i <= 2*num_param
            prop = rand(MvNormal(old_p, sigmaMat1), 1)
        end
        if i > 2*num_param
            cov_dat = cov(x[1:(i-1),:]) + 0.0000001*I
            sigmaMat2 = (cov_dat*(2.38)^2)/num_param
            prop = (1-0.05)*rand(MvNormal(old_p, sigmaMat2), 1) + 0.05*rand(MvNormal(old_p, sigmaMat1), 1)
        end 
        waicold = waicMat[i-1,:]
        prior2 = LogPrior(prop)
        if prior2 == -Inf
            x[i,:] = x[i-1,:]
            waicMat[i,:] = waicold
            loglik[i] = loglik[i-1]
        else 
            prior1 = LogPrior(old_p)
            llprop = likelihood(prop)
            llold = loglik[i-1]
            ratio = exp(llprop[tmax+1] + prior2 - llold - prior1)
            if rand() < ratio
                x[i,:] = prop
                waicMat[i,:] = llprop
                loglik[i] = llprop[tmax+1]
            else
                x[i,:] = old_p
                waicMat[i,:] = waicold
                loglik[i] = llold
                #println("Rejected")
            end
            #println(waicMat[i,:])
            #println("Loglikelihood  = ", loglik[i])
        end
    end
    return x, waicMat
end