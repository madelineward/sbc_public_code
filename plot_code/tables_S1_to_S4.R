# Table S1

num_pop <- 20

param_estimates <- read.table("output/params_full.txt", sep = ",", fill = T)
names(param_estimates) <- c("param", "median", "lower", "upper", "mod_gen", "mod_fit", "scen", "pop")
param_estimates <- param_estimates[which(param_estimates$mod_gen == param_estimates$mod_fit),]
head(param_estimates)

true_param_values <- read.csv("data/true_param_values.csv", header = T)
head(true_param_values)

param_estimates <- merge(param_estimates, true_param_values, 
                         by = c("param", "mod_gen", "scen"))

param_estimates$captured <- ifelse(param_estimates$lower < param_estimates$true_value & param_estimates$upper > param_estimates$true_value, 1, 0)

param_estimates$degrees <- ifelse(param_estimates$scen == "a22b2d50h40" |
                                    param_estimates$scen == "a22b2d10h40" |
                                    param_estimates$scen == "a24b2d005" | 
                                    param_estimates$scen == "a24b2d001" |
                                    param_estimates$scen == "a24b2d80d02" |
                                    param_estimates$scen == "a24b2d04d005" |
                                    (param_estimates$scen == "a24b2d3d10" &
                                       param_estimates$mod_gen == "mod4a") |
                                    (param_estimates$scen == "a24b2d3d20" &
                                       param_estimates$mod_gen == "mod4b"), "Weak",
                                  ifelse(param_estimates$scen == "a22b2d65h40" |
                                           param_estimates$scen == "a22b2d15h40" |
                                           param_estimates$scen == "a24b2d01" |
                                           param_estimates$scen == "a24b2d0015" |
                                           param_estimates$scen == "a24b2d80d03" |
                                           param_estimates$scen == "a24b2d04d007" | 
                                           param_estimates$scen == "a24b2d3d15" |
                                           param_estimates$scen == "a24b2d3d075",
                                         "Medium", "Strong"))

mod1a <- param_estimates[which(param_estimates$mod_gen == "mod1a"),]
aggregate(captured ~ param  + degrees, data = mod1a, FUN = function(x) sum(x)/20)

# Table S2

mod2a <- param_estimates[which(param_estimates$mod_gen == "mod2a"),]
aggregate(captured ~ param  + degrees, data = mod2a, FUN = function(x) sum(x)/20)

# Table S3

mod3a <- param_estimates[which(param_estimates$mod_gen == "mod3a"),]
aggregate(captured ~ param  + degrees, data = mod3a, FUN = function(x) sum(x)/20)

# Table S4

mod4a <- param_estimates[which(param_estimates$mod_gen == "mod4a"),]
aggregate(captured ~ param  + degrees, data = mod4a, FUN = function(x) sum(x)/20)
