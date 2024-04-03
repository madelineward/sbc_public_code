# Figure S1

num_pop <- 20
n_ticks <- num_pop

param_estimates <- read.table("output/params_full.txt", sep = ",", fill = T)
names(param_estimates) <- c("param", "median", "lower", "upper", "mod_gen", "mod_fit", "scen", "pop")
param_estimates <- param_estimates[which(param_estimates$mod_gen == param_estimates$mod_fit),]
head(param_estimates)

true_param_values <- read.csv("data/true_param_values.csv", header = T)
head(true_param_values)

param_estimates <- merge(param_estimates, true_param_values, 
                         by = c("param", "mod_gen", "scen"))
head(param_estimates)

scens_mod1a <- unique(param_estimates$scen[which(param_estimates$mod_gen == "mod1a")])
mod1a_w <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[1]),]
mod1a_m <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[2]),]
mod1a_s <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[3]),]

titles <- c(", Weak BC", ", Medium BC", ", Strong BC")

par(mfrow = c(4,3))

for(i in 1:3) {
  mod1a <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[i] &
                                   param_estimates$param == "alpha"),]
  t <- titles[i]
  plot(1:n_ticks, mod1a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("α", t), 
       col = "darkblue", 
       ylim = c(min(mod1a$lower),
                max(mod1a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod1a$lower, 
         1:n_ticks, mod1a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod1a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod1a <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[i] &
                                   param_estimates$param == "beta"),]
  t <- titles[i]
  plot(1:n_ticks, mod1a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("β", t), 
       col = "darkblue", 
       ylim = c(min(mod1a$lower),
                max(mod1a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod1a$lower, 
         1:n_ticks, mod1a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod1a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod1a <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[i] &
                                   param_estimates$param == "delta1"),]
  t <- titles[i]
  plot(1:n_ticks, mod1a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ1", t), 
       col = "darkblue", 
       ylim = c(min(mod1a$lower),
                max(mod1a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod1a$lower, 
         1:n_ticks, mod1a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod1a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod1a <- param_estimates[which(param_estimates$mod_gen == "mod1a" &
                                   param_estimates$scen == scens_mod1a[i] &
                                   param_estimates$param == "delta2"),]
  t <- titles[i]
  plot(1:n_ticks, mod1a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ2", t), 
       col = "darkblue", 
       ylim = c(min(mod1a$lower),
                max(mod1a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod1a$lower, 
         1:n_ticks, mod1a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod1a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}


# Figure S2

scens_mod2a <- unique(param_estimates$scen[which(param_estimates$mod_gen == "mod2a")])
mod2a_w <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[1]),]
mod2a_m <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[2]),]
mod2a_s <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[3]),]

titles <- c(", Weak BC", ", Medium BC", ", Strong BC")

par(mfrow = c(3,3))

for(i in 1:3) {
  mod2a <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[i] &
                                   param_estimates$param == "alpha"),]
  t <- titles[i]
  plot(1:n_ticks, mod2a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("α", t), 
       col = "darkblue", 
       ylim = c(min(mod2a$lower),
                max(mod2a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod2a$lower, 
         1:n_ticks, mod2a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod2a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod2a <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[i] &
                                   param_estimates$param == "beta"),]
  t <- titles[i]
  plot(1:n_ticks, mod2a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("β", t), 
       col = "darkblue", 
       ylim = c(min(mod2a$lower),
                max(mod2a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod2a$lower, 
         1:n_ticks, mod2a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod2a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod2a <- param_estimates[which(param_estimates$mod_gen == "mod2a" &
                                   param_estimates$scen == scens_mod2a[i] &
                                   param_estimates$param == "delta1"),]
  t <- titles[i]
  plot(1:n_ticks, mod2a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ1", t), 
       col = "darkblue", 
       ylim = c(min(mod2a$lower),
                max(mod2a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod2a$lower, 
         1:n_ticks, mod2a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod2a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}


# Figure S3

scens_mod3a <- unique(param_estimates$scen[which(param_estimates$mod_gen == "mod3a")])
mod3a_w <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[1]),]
mod3a_m <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[2]),]
mod3a_s <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[3]),]

titles <- c(", Weak BC", ", Medium BC", ", Strong BC")

par(mfrow = c(4,3))

for(i in 1:3) {
  mod3a <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[i] &
                                   param_estimates$param == "alpha"),]
  t <- titles[i]
  plot(1:n_ticks, mod3a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("α", t), 
       col = "darkblue", 
       ylim = c(min(mod3a$lower),
                max(mod3a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod3a$lower, 
         1:n_ticks, mod3a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod3a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod3a <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[i] &
                                   param_estimates$param == "beta"),]
  t <- titles[i]
  plot(1:n_ticks, mod3a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("β", t), 
       col = "darkblue", 
       ylim = c(min(mod3a$lower),
                max(mod3a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod3a$lower, 
         1:n_ticks, mod3a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod3a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod3a <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[i] &
                                   param_estimates$param == "delta2"),]
  t <- titles[i]
  plot(1:n_ticks, mod3a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ1", t), 
       col = "darkblue", 
       ylim = c(min(mod3a$lower),
                max(mod3a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod3a$lower, 
         1:n_ticks, mod3a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod3a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod3a <- param_estimates[which(param_estimates$mod_gen == "mod3a" &
                                   param_estimates$scen == scens_mod3a[i] &
                                   param_estimates$param == "delta1"),]
  t <- titles[i]
  plot(1:n_ticks, mod3a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ2", t), 
       col = "darkblue", 
       ylim = c(min(mod3a$lower),
                max(mod3a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod3a$lower, 
         1:n_ticks, mod3a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod3a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

# Figure S4

scens_mod4a <- rev(unique(param_estimates$scen[which(param_estimates$mod_gen == "mod4a")]))

titles <- c(", Weak BC", ", Medium BC", ", Strong BC")

par(mfrow = c(4,3))

for(i in 1:3) {
  mod4a <- param_estimates[which(param_estimates$mod_gen == "mod4a" &
                                   param_estimates$scen == scens_mod4a[i] &
                                   param_estimates$param == "alpha"),]
  t <- titles[i]
  plot(1:n_ticks, mod4a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("α", t), 
       col = "darkblue", 
       ylim = c(min(mod4a$lower),
                max(mod4a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod4a$lower, 
         1:n_ticks, mod4a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod4a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod4a <- param_estimates[which(param_estimates$mod_gen == "mod4a" &
                                   param_estimates$scen == scens_mod4a[i] &
                                   param_estimates$param == "beta"),]
  t <- titles[i]
  plot(1:n_ticks, mod4a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("β", t), 
       col = "darkblue", 
       ylim = c(min(mod4a$lower),
                max(mod4a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod4a$lower, 
         1:n_ticks, mod4a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod4a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod4a <- param_estimates[which(param_estimates$mod_gen == "mod4a" &
                                   param_estimates$scen == scens_mod4a[i] &
                                   param_estimates$param == "delta2"),]
  t <- titles[i]
  plot(1:n_ticks, mod4a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ1", t), 
       col = "darkblue", 
       ylim = c(min(mod4a$lower),
                max(mod4a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod4a$lower, 
         1:n_ticks, mod4a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod4a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}

for(i in 1:3) {
  mod4a <- param_estimates[which(param_estimates$mod_gen == "mod4a" &
                                   param_estimates$scen == scens_mod4a[i] &
                                   param_estimates$param == "delta2"),]
  t <- titles[i]
  plot(1:n_ticks, mod4a$median, 
       xlab = "Population", ylab= "Value", 
       main= paste0("δ2", t), 
       col = "darkblue", 
       ylim = c(min(mod4a$lower),
                max(mod4a$upper)), 
       pch = 19,xaxt = "n")
  arrows(1:n_ticks, mod4a$lower, 
         1:n_ticks, mod4a$upper, 
         length=0.05, angle=90, code=3, lwd = 1)
  abline(h = mod4a$true_value[1],
         col = "#9c0534", lwd = 1)
  axis(1, at=1:n_ticks, labels=1:20)
}
