######### MODEL SIMULATION AND LIKELIHOOD FUNCTIONS ########

########## Baseline SILM ##########

## Simulate data from baseline SILM
# Param should vector [alpha beta]
function sim_data_base(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
	gt = repeat([0.0], outer = [tmax])
    for t in 1:tmax
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end       
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    #gt = 1
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_base(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status

	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    for t in (tmin-1):tmax
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end

# function to calculate the likelihood of a baseline SILM with only alpha and beta in an SIR framework
# returns log likelihood
function likelihood_base(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)
	for t in 1:tmax
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
    lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end


########## Model 1A ##########

## Simulate data from Model 1A (threshold alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod1a(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
    H = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        if prev_prev > H
            gt[t] = delta
        end
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod1a(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
    H = param[4]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status

	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        #println("t = ", t)
        if prev_prev > H
            gt[t] = delta
        end
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

# function to calculate the likelihood of Model 1A in an SIR framework
# returns log likelihood
function likelihood_mod1a(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    delta = param[3]
    H = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        if prev_prev > H
            gt[t] = delta
            #print("gt[t]: ", gt[t])
        end
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end


########## Model 1B ##########

## Simulate data from Model 1B (threshold alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod1b(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
    H = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        if prev_prev > H
            gt[t] = delta
        end
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 


function forecast_mod1b(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
    H = param[4]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        #println("t = ", t)
        if prev_prev > H
            gt[t] = delta
        end
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

# function to calculate the likelihood of Model 1B
# returns log likelihood
function likelihood_mod1b(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    delta = param[3]
    H = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        if prev_prev > H
            gt[t] = delta
            #print("gt[t]: ", gt[t])
        end
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end

########## Model 2A ##########

## Simulate data from Model 2A (exponential alarm)
# Param should vector [alpha beta delta1]
function sim_data_mod2a(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = (1-exp(-delta*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod2a(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status

	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        gt[t] = (1-exp(-delta*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 2A
# Returns Log-likelihood
function likelihood_mod2a(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d = param[3]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = (1-exp(-d*prev_prev))
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end



########## Model 2B ##########

## Simulate data from Model 2B (exponential alarm)
# Param should vector [alpha beta delta1]
function sim_data_mod2b(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = (1-exp(-delta*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end

function forecast_mod2b(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    delta = param[3]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        gt[t] = (1-exp(-delta*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 2B
# Returns Log-likelihood
function likelihood_mod2b(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d = param[3]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = (1-exp(-d*prev_prev))
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end

########## Model 3A ##########

## Simulate data from Model 3A (scaled exponential alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod3a(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod3a(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status

	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 3A
# Returns Log-likelihood
function likelihood_mod3a(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d1 = param[3]
    d2 = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end

########## Model 3B ##########

## Simulate data from Model 3B (scaled exponential alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod3b(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod3b(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[3]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)
    for t in (tmin-1):tmax
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 3B
# Returns Log-likelihood
function likelihood_mod3b(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d1 = param[3]
    d2 = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = d1*(1-exp(-d2*prev_prev))
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end

########## Model 4A ##########

## Simulate data from Model 4A (Hill-type alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod4a(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)/num_ind
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod4a(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)/num_ind
    for t in (tmin-1):tmax
        #println("t = ", t)
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^b
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*(1-gt[t])*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)/num_ind
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 4A
# Returns Log-likelihood
function likelihood_mod4a(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d1 = param[3]
    d2 = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b)
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*(1-gt[t])*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )/num_ind
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end

########## Model 4B ##########

## Simulate data from Model 4B (Hill-type alarm)
# Param should vector [alpha beta delta1 delta2]
function sim_data_mod4b(seed, num_ind, x, y, inf_per, init_inf, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	status[init_inf] .= 1                  	# This allows for initial infection of 1 location
	inftime = repeat([0], outer = [num_ind])
	inftime[init_inf] .= 1
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    for t in 1:tmax
        # println(prev_prev)
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)/num_ind
        
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(x = x, y = y, status = status, inftime = inftime, rem_time = rem_time), gt
end 

function forecast_mod4b(seed, num_ind, x, y, inf_per, init_inf, tmin, tmax, param)
	Random.seed!(seed)                      # Set random seed for reproducibilty
    alpha = param[1]
    beta = param[2]
    d1 = param[3]
    d2 = param[4]
	inftime = init_inf
    status = repeat([0], outer = [num_ind]) # Vector to keep track of infection status
	for i in 1:num_ind
        if inftime[i] <= tmin - inf_per && inftime[i] > 0
            status[i] = 2
        elseif inftime[i] > tmin - inf_per && inftime[i] < tmin
            status[i] = 1
        end
    end
	new_status = repeat([0], outer= [num_ind])     # Status update vector
    prev_status = repeat([0], outer= [num_ind])
    rem_time = repeat([0], outer = [num_ind])      # Removal time vector 
    a = -alpha
    b = -beta
    gt = repeat([0.0], outer = [tmax])
    prev_prev = sum(status .== 1)/num_ind
    for t in (tmin-1):tmax
        #println("t = ", t)
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if (status[i] == 0) 
				dist_count = 0
				for j in 1:num_ind
				    if j != i && status[j] ==1 
                        d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
                        dist_count = dist_count + d_ij
				    end
				end 
				p_it = 1-exp(a*dist_count)
				if p_it >= rand() 
					new_status[i] = 1
                    prev_status[i] = 0
					inftime[i] = t + 1
                else
					new_status[i] = 0
                    prev_status[i] = 0
					inftime[i] = 0
				end 
			end
			if status[i] == 1 && t >= inftime[i] + inf_per 
				new_status[i] = 2
                prev_status[i] = 2
				rem_time[i] = t
			elseif status[i] == 1 && t < inftime[i] + inf_per
				new_status[i] = 1
                prev_status[i] = 1
				rem_time[i] = 0
			end
            if status[i] == 2
                new_status[i] = 2
            end
        end
        #println(sum(new_status .== 1))
        for i in 1:num_ind
            status[i] = new_status[i]   
        end 
        prev_prev = sum(prev_status .== 1)/num_ind
        #println("prev = ", prev_prev)
    end
    for i in 1:num_ind
        if inftime[i] > 0
            rem_time[i] = inftime[i] + inf_per
        end
    end
    return DataFrame(inftime = inftime), gt
end 

#Calculates likelihood for Model 4B
# Returns Log-likelihood
function likelihood_mod4b(inftime, rem_time, tmax, param, x, y, num_ind)
	a = param[1]
    b = -param[2]
    d1 = param[3]
    d2 = param[4]
    gt = repeat([0.0], outer = [tmax])
    prev_prev = 0
    sump_it = 0
    sump_it2 = 0
    lnlike = zeros(tmax + 1)

	for t in 1:tmax
        gt[t] = prev_prev^d1/(d2^d1 + prev_prev^d1)
		for i in 1:num_ind
			if inftime[i] == t + 1
				dist_count = 0
				for j in 1:num_ind
				    if j != i && inftime[j] <= t && rem_time[j] > t 
					    d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					    dist_count = dist_count + d_ij
				    end 
				end 
				p_it = 1-exp(-a*dist_count)
				lp_it = log(p_it)
			else
				lp_it = 0
			end
			sump_it = sump_it + lp_it
		end 
		for i in 1:num_ind
			if inftime[i] > t + 1 || inftime[i] == 0
				dist_count = 0
				for j in 1:num_ind
				if j != i && inftime[j] <= t && rem_time[j] > t 
					d_ij = (sqrt((x[i]-x[j])^2 + (y[i] - y[j])^2)+1)^(b*((1-gt[t])^(-1)))
					dist_count = dist_count + d_ij
				end 
				end 

				p_it2 = 1-exp(-a*dist_count)
				lp_it2 = log(1 - p_it2)
			else
				lp_it2 = 0
			end 
			sump_it2 = sump_it2 + lp_it2
		end
	lnlike[t] = sump_it + sump_it2
	sump_it = 0
	sump_it2 = 0
    prev_prev = sum( (inftime .<= t) .& (rem_time .> t) )/num_ind
    #println(prev_prev)
	end
    lnlike[tmax + 1] = sum(lnlike[1:tmax])
    return lnlike
end