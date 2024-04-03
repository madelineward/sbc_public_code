true_bc_curves <- read.table("data/true_curves.txt", sep = ",")
names(true_bc_curves) <- c("pop", "model", "scen", "deg", paste("t", 1:31, sep = ""))
head(true_bc_curves)

true_param_values <- read.csv("data/true_param_values.csv", header = T)
head(true_param_values)

true_params_wide <- reshape(true_param_values, 
                            idvar = c("mod_gen","scen"),
                            timevar = "param",
                            direction = "wide")
names(true_params_wide) <- c("model", "scen", "alpha", "beta", "delta1", "delta2")
head(true_params_wide)

true_bc_curves <- merge(true_bc_curves, true_params_wide, by = c("model", "scen"))

my_cols <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')

lty1 <- "solid"
lty2 <- "dotted"
lty3 <- "dashed"
lty4 <- "dotdash"
lty5 <- "longdash"
lty6 <- "twodash"

### Model 1 Alarm Function (Threshold Alarm)

model1 <- function(x, params) {
  delta <- params[3]
  H <- params[4]
  return(ifelse(x < H, 0, delta))
}

### Model 2 Alarm Function (Exponential Alarm)

model2 <- function(x, params) {
  delta <- params[3]
  return(1-exp(-delta*x))
}

### Model 3 Alarm Function (Scaled Exponential Alarm)

model3 <- function(x, params) {
  d1 <- params[3]
  d2 <- params[4]
  return(d1*(1-exp(-d2*x)))
}

### Model 4 Alarm Function (Hill-Type Alarm)

model4 <- function(x, params) {
  d1 <- params[1]
  d2 <- params[2]
  return(x^d1/(d2^d1 + x^d1))
}

all_alarm_funs <- list(mod1a = model1, mod2a = model2, mod3a = model3, mod4a = model4)

true_alarms <- matrix(nrow = nrow(true_bc_curves), ncol = 30)

for(i in 1:nrow(true_bc_curves)) {
  params <- true_bc_curves[i,37:40]
  mod_ind <- ifelse(true_bc_curves$model[i] == "mod1a" | true_bc_curves$model[i] == "mod1b", 1, ifelse(true_bc_curves$model[i] == "mod2a" | true_bc_curves$model[i] == "mod2b", 2, ifelse(true_bc_curves$model[i] == "mod3a" | true_bc_curves$model[i] == "mod3b",3,4)))
  alarm_fun <- all_alarm_funs[[mod_ind]]
  true_alarms[i,] <- alarm_fun(true_bc_curves[i,5:34], params)
}
