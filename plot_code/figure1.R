# Figure 1 from paper

library(ggplot2)
library(tidyr)
library(gridExtra)

my_cols <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')

lty1 <- "solid"
lty2 <- "dotted"
lty3 <- "dashed"
lty4 <- "dotdash"
lty5 <- "longdash"
lty6 <- "twodash"

### Model 1 Alarm Function

model1 <- function(x, delta, H) {
  return(ifelse(x < H, 0, delta))
}

### Model 1A Alarm Functions (Threshold Alarm)

x1 <- seq(0,500, by = 0.1)
# Weak: Delta1 = 0.5, Delta2/H = 40
yA1a <- model1(x1, 0.50, 100)
# Medium: Delta1 = 0.65, Delta2/H = 40
yA1b <- model1(x1, 0.65, 100)
# Strong: Delta1 = 0.8, Delta2/H = 40
yA1c <- model1(x1, 0.80, 100)

### Model 1B Alarm Functions (Threshold Alarm)

x1 <- seq(0,500, by = 0.1)
# Weak: Delta1 = 0.1, Delta2/H = 40
yB1a <- model1(x1, 0.10, 40)
# Medium: Delta1 = 0.15, Delta2/H = 40
yB1b <- model1(x1, 0.15, 40)
# Strong: Delta1 = 0.2, Delta2/H = 40
yB1c <- model1(x1, 0.20, 40)

mod1_alarms = data.frame(x1, yA1a, yA1b, yA1c, yB1a, yB1b, yB1c)
mod1_alarms <- mod1_alarms %>% 
  pivot_longer(-x1, names_to = "variable", values_to = "value")

labels_mod1 <- c(expression(paste(delta[1], "= 0.50")), 
                 expression(paste(delta[1], "= 0.65")), 
                 expression(paste(delta[1], "= 0.80")), 
                 expression(paste(delta[1], "= 0.10")), 
                 expression(paste(delta[1], "= 0.15")), 
                 expression(paste(delta[1], "= 0.20")))

mod1_alarms_plot <- ggplot(data = mod1_alarms, 
                           aes(x1, value, colour = variable, linetype = variable)) +
  geom_line(size = 1.5) +
  scale_color_manual(name = "Scenario", labels = labels_mod1, 
                     values = c(my_cols[1], my_cols[2], my_cols[3], 
                                my_cols[4], my_cols[5], my_cols[6])) +
  scale_linetype_manual(name = "Scenario", labels = labels_mod1, 
                        values = c(lty1, lty2, lty3, lty4, lty5, lty6)) +
  theme_bw() +
  ylim(0,1) +
  labs(title = "Threshold Alarm (Model 1)",
       x = "Prevalence at t-1",
       y = "Alarm Value") +
  theme( 
    plot.title = element_text(size=25),
    axis.title = element_text(size=22),
    axis.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))
mod1_alarms_plot 

### Model 2 Alarm Function (Exponential Alarm)

model2 <- function(x, delta) {
  return(1-exp(-delta*x))
}

### Model 2A Alarm Functions (Exponential Alarm)

x2 <- seq(0,500, by = 0.1)
# Weak: Delta = 0.005
yA2a <- model2(x2, 0.005)
# Medium: Delta = 0.01
yA2b <- model2(x2, 0.01)
# Strong: Delta = 0.015
yA2c <- model2(x2, 0.015)


### Model 2B Alarm Functions (Exponential Alarm)

x2 <- seq(0,500, by = 0.1)
# Weak: Delta = 0.001
yB2a <- model2(x2, 0.001)
# Medium: Delta = 0.0015
yB2b <- model2(x2, 0.0015)
# Strong: Delta = 0.002
yB2c <- model2(x2, 0.002)


mod2_alarms = data.frame(x2, yA2a, yA2b, yA2c, yB2a, yB2b, yB2c)
mod2_alarms <- mod2_alarms %>% 
  pivot_longer(-x2, names_to = "variable", values_to = "value")

labels_mod2 <- c(expression(paste(delta[1], "= 0.005")), 
                 expression(paste(delta[1], "= 0.01")), 
                 expression(paste(delta[1], "= 0.015")), 
                 expression(paste(delta[1], "= 0.001")), 
                 expression(paste(delta[1], "= 0.0015")), 
                 expression(paste(delta[1], "= 0.002")))

mod2_alarms_plot <- ggplot(data = mod2_alarms, 
                           aes(x2, value, colour = variable, linetype = variable)) +
  geom_line(size = 1.5) +
  scale_color_manual(name = "Scenario", labels = labels_mod2, 
                     values = c(my_cols[1], my_cols[2], my_cols[3], 
                                my_cols[4], my_cols[5], my_cols[6])) +
  scale_linetype_manual(name = "Scenario", labels = labels_mod2, 
                        values = c(lty1, lty2, lty3, lty4, lty5, lty6)) +
  theme_bw() +
  labs(title = "Exponential Alarm (Model 2)",
       x = "Prevalence at t-1",
       y = "Alarm Value") +
  theme( 
    plot.title = element_text(size=25),
    axis.title = element_text(size=22),
    axis.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))
mod2_alarms_plot 


### Model 3 Alarm Function (Scaled Exponential Alarm)

model3 <- function(x, d1, d2) {
  return(d1*(1-exp(-d2*x)))
}

### Model 3A Alarm Functions (Scaled Exponential Alarm)

x3 <- seq(0,500, by = 0.1)
# Weak: Delta = 0.001
yA3a <- model3(x3, 0.80, 0.02)
# Medium: Delta = 0.0015
yA3b <- model3(x3, 0.80, 0.03)
# Strong: Delta = 0.002
yA3c <- model3(x3, 0.80, 0.04)

### Model 3B Alarm Functions (Scaled Exponential Alarm)

x3 <- seq(0,500, by = 0.1)
# Weak: Delta = 0.001
yB3a <- model3(x3, 0.40, 0.005)
# Medium: Delta = 0.0015
yB3b <- model3(x3, 0.40, 0.007)
# Strong: Delta = 0.002
yB3c <- model3(x3, 0.40, 0.009)

mod3_alarms = data.frame(x3, yA3a, yA3b, yA3c, yB3a, yB3b, yB3c)
mod3_alarms <- mod3_alarms %>% 
  pivot_longer(-x3, names_to = "variable", values_to = "value")

labels_mod3 <- c(expression(paste(delta[1], "= 0.02")), 
                 expression(paste(delta[1], "= 0.03")), 
                 expression(paste(delta[1], "= 0.04")), 
                 expression(paste(delta[1], "= 0.005")), 
                 expression(paste(delta[1], "= 0.007")), 
                 expression(paste(delta[1], "= 0.009")))

mod3_alarms_plot <- ggplot(data = mod3_alarms, 
                           aes(x3, value, colour = variable, linetype = variable)) +
  geom_line(size = 1.5) +
  scale_color_manual(name = "Scenario", labels = labels_mod3, 
                     values = c(my_cols[1], my_cols[2], my_cols[3], 
                                my_cols[4], my_cols[5], my_cols[6])) +
  scale_linetype_manual(name = "Scenario", labels = labels_mod3, 
                        values = c(lty1, lty2, lty3, lty4, lty5, lty6)) +
  theme_bw() +
  ylim(0,1) +
  labs(title = "Scaled Exponential Alarm (Model 3)",
       x = "Prevalence at t-1",
       y = "Alarm Value") +
  theme( 
    plot.title = element_text(size=25),
    axis.title = element_text(size=22),
    axis.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))
mod3_alarms_plot 

### Model 4 Alarm Function (Hill-Type Alarm)

model4 <- function(x, d1, d2) {
  return(x^d1/(d2^d1 + x^d1))
}

### Model 4B Alarm Functions 
x4 <- seq(0,0.5, by = 0.001)
# Weak: Delta = 0.001
yA4a <- model4(x4, 3, 0.10)
# Medium: Delta = 0.0015
yA4b <- model4(x4, 3, 0.075)
# Strong: Delta = 0.002
yA4c <- model4(x4, 3, 0.05)

### Model 4B Alarm Functions 
x4 <- seq(0,0.5, by = 0.001)
# Weak: Delta = 0.001
yB4a <- model4(x4, 3, 0.20)
# Medium: Delta = 0.0015
yB4b <- model4(x4, 3, 0.15)
# Strong: Delta = 0.002
yB4c <- model4(x4, 3, 0.10)


mod4_alarms = data.frame(x4, yA4a, yA4b, yA4c, yB4a, yB4b, yB4c)
mod4_alarms <- mod4_alarms %>% 
  pivot_longer(-x4, names_to = "variable", values_to = "value")

labels_mod4 <- c(expression(paste(delta[1], "= 0.10")), 
                 expression(paste(delta[1], "= 0.075")), 
                 expression(paste(delta[1], "= 0.05")), 
                 expression(paste(delta[1], "= 0.20")), 
                 expression(paste(delta[1], "= 0.15")), 
                 expression(paste(delta[1], "= 0.10")))

mod4_alarms_plot <- ggplot(data = mod4_alarms, 
                           aes(x4, value, colour = variable, linetype = variable)) +
  geom_line(size = 1.5) +
  scale_color_manual(name = "Scenario", labels = labels_mod4, 
                     values = c(my_cols[1], my_cols[2], my_cols[3], 
                                my_cols[4], my_cols[5], my_cols[6])) +
  scale_linetype_manual(name = "Scenario", labels = labels_mod4, 
                        values = c(lty1, lty2, lty3, lty4, lty5, lty6)) +
  theme_bw() +
  ylim(0,1) +
  labs(title = "Hill-Type Alarm (Model 4)",
       x = "Proportion Infectious at t-1",
       y = "Alarm Value") +
  theme( 
    plot.title = element_text(size=25),
    axis.title = element_text(size=22),
    axis.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))
mod4_alarms_plot 

grid.arrange(mod1_alarms_plot, mod2_alarms_plot, mod3_alarms_plot, mod4_alarms_plot)
