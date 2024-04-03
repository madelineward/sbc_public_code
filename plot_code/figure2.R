# Figure 2a)

library(ggplot2)
library(tibble)

# Note need to run combine_true_curves.R first

true_bc_curves <- read.table("data/true_curves.txt", sep = ",")
names(true_bc_curves) <- c("pop", "model", "scen", "deg", paste("t", 1:31, sep = ""))
head(true_bc_curves)

true_nobc_curves <- read.table("data/true_curves_nobc.txt", sep = ",")
names(true_nobc_curves) <- c("pop", "model", paste("t", 1:31, sep = ""))
head(true_nobc_curves)
par(mfrow = c(2,2))
popi = 8

# MODEL 1A
plot(1:31, true_nobc_curves[which(true_nobc_curves$model == "mod1a" & true_nobc_curves$pop == popi),3:33], type = "l", col ='#EE7733', lwd = 3, main = "Model 1A (Threshold Alarm)", xlab = "Time", ylab = "New Infections", ylim = c(0,160))
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod1a" & true_bc_curves$deg == "weak" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#0077BB')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod1a" & true_bc_curves$deg == "med" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#33BBEE')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod1a" & true_bc_curves$deg == "str" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#EE3377')
legend(16, 140, c("No BC", expression(paste(delta[1], " = 0.50, ", delta[2], " = 40")), expression(paste(delta[1], " = 0.65, ", delta[2], " = 40")), expression(paste(delta[1], " = 0.80, ", delta[2], " = 40"))), col = c( '#EE7733', '#0077BB', '#33BBEE', '#EE3377'), lty = c(1,1,1,1), lwd = c(3,3,3,3), bty = "n")
# MODEL 2A
plot(1:31, true_nobc_curves[which(true_nobc_curves$model == "mod2a" & true_nobc_curves$pop == popi),3:33], type = "l", col ='#EE7733', lwd = 3, main = "Model 2A (Exponential Alarm)", xlab = "Time", ylab = "New Infections", ylim = c(0,160))
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod2a" & true_bc_curves$deg == "weak" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#0077BB')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod2a" & true_bc_curves$deg == "med" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#33BBEE')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod2a" & true_bc_curves$deg == "str" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#EE3377')
legend(16, 140, c("No BC", expression(paste(delta[1], " = 0.005")), expression(paste(delta[1], " = 0.010")), expression(paste(delta[1], " = 0.015"))), col = c( '#EE7733', '#0077BB', '#33BBEE', '#EE3377'), lty = c(1,1,1,1), lwd = c(3,3,3,3), bty = "n")
# MODEL 3A
plot(1:31, true_nobc_curves[which(true_nobc_curves$model == "mod3a" & true_nobc_curves$pop == popi),3:33], type = "l", col ='#EE7733', lwd = 3, main = "Model 3A (Scaled Exponential Alarm)", xlab = "Time", ylab = "New Infections", ylim = c(0,160))
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod3a" & true_bc_curves$deg == "weak" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#0077BB')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod3a" & true_bc_curves$deg == "med" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#33BBEE')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod3a" & true_bc_curves$deg == "str" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#EE3377')
legend(16, 140, c("No BC", expression(paste(delta[1], " = 0.02, ", delta[2], " = 0.80")), expression(paste(delta[1], " = 0.03, ", delta[2], " = 0.80")), expression(paste(delta[1], " = 0.04, ", delta[2], " = 0.80"))), col = c( '#EE7733', '#0077BB', '#33BBEE', '#EE3377'), lty = c(1,1,1,1), lwd = c(3,3,3,3), bty = "n")
# MODEL 4A
plot(1:31, true_nobc_curves[which(true_nobc_curves$model == "mod4a" & true_nobc_curves$pop == popi),3:33], type = "l", col ='#EE7733', lwd = 3, main = "Model 4A (Hill-Type Alarm)", xlab = "Time", ylab = "New Infections", ylim = c(0,160))
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod4a" & true_bc_curves$deg == "weak" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#0077BB')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod4a" & true_bc_curves$deg == "med" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#33BBEE')
lines(1:31, true_bc_curves[which(true_bc_curves$model == "mod4a" & true_bc_curves$deg == "str" & true_bc_curves$pop == popi),5:35], lwd = 3, col = '#EE3377')
legend(16, 140, c("No BC", expression(paste(delta[1], " = 0.10, ", delta[2], " = 3")), expression(paste(delta[1], " = 0.075, ", delta[2], " = 3")), expression(paste(delta[1], " = 0.050, ", delta[2], " = 3"))), col = c( '#EE7733', '#0077BB', '#33BBEE', '#EE3377'), lty = c(1,1,1,1), lwd = c(3,3,3,3), bty = "n")


# Figure 2b)

true_bc_curves$size <- rowSums(true_bc_curves[,5:ncol(true_bc_curves)]) 
true_nobc_curves$size <- rowSums(true_nobc_curves[,3:ncol(true_nobc_curves)])

bc_size_means <- aggregate(true_bc_curves$size, list(true_bc_curves$model, true_bc_curves$deg), FUN=mean) 
names(bc_size_means) <- c("model", "degree", "mean_size")

nobc_size_means <- aggregate(true_nobc_curves$size, list(true_nobc_curves$model), FUN=mean) 
nobc_size_means <- add_column(nobc_size_means, "degree" = rep("no",8), .after = "Group.1")
names(nobc_size_means) <- c("model", "degree", "mean_size")

bc_size_sds <- aggregate(true_bc_curves$size, list(true_bc_curves$model, true_bc_curves$deg), FUN=sd) 
names(bc_size_sds) <- c("model", "degree", "size_sd")
bc_size_means$size_sd <- bc_size_sds$size_sd

nobc_size_sds <- aggregate(true_nobc_curves$size, list(true_nobc_curves$model), FUN=sd) 
nobc_size_sds <- add_column(nobc_size_sds, "degree" = rep("no",8), .after = "Group.1")
names(nobc_size_sds) <- c("model", "degree", "size_sd")
nobc_size_means$size_sd <- nobc_size_sds$size_sd

size_means <- rbind.data.frame(bc_size_means, nobc_size_means)
size_means$size_se <- size_means$size_sd/sqrt(20)
size_means$degree <- factor(size_means$degree, levels = c("str", "med", "weak", "no"))

deg_labels <- c("Strong", "Medium", "Weak", "None")

ggplot(size_means, aes(fill=degree, y=mean_size, x=model)) + 
  geom_bar(position="dodge", stat="identity", alpha = 0.75) + 
  scale_fill_manual(name = "Degree of BC", values=rev(c('#EE7733', '#0077BB', '#33BBEE', '#EE3377')), labels = deg_labels) +
  theme_minimal() + 
  coord_flip() + 
  ylab("Mean Epidemic Size (+/-SD)") + 
  xlab("Model") +
  geom_errorbar(aes(ymin=mean_size-size_sd, ymax=mean_size+size_sd), width=.2,
                position=position_dodge(.9)) + 
  theme( 
    axis.title = element_text(size=14),
    axis.text = element_text(size=12),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "top") +
  scale_x_discrete(labels=c("mod1a" = "1A", "mod1b" = "1B",
                            "mod2a" = "2A", "mod2b" = "2B",
                            "mod3a" = "3A", "mod3b" = "3B",
                            "mod4a" = "4A", "mod4b" = "4B"))
