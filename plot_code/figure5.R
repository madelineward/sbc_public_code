# Figure 5 

tmax <- 30

full_curves <- read.table("output/curves_full.txt", sep = ",", fill = T)
names(full_curves) <- c("time", "median", "lower", "upper", "pop", "mod_gen", "scen", "mod_fit")
head(full_curves)
full_curves <- full_curves[-61474,]
full_curves[2:4] <- sapply(full_curves[2:4],as.numeric)
summary(full_curves)

# Remove models/populations that didn't converge
not_converged <- read.table("output/non_converged_models.txt", sep = ",")
names(not_converged) <- c("mod_gen", "mod_fit", "scen", "pop")
not_converged$converge <- 1

full_curves <- merge(full_curves, not_converged, by = c("mod_gen", "mod_fit", "scen", "pop"), all.x = T)
head(full_curves)

full_curves <- full_curves[-which(full_curves$converge == 1),1:8]

# Read in true curves
true_bc_curves <- read.table("data/true_curves.txt", sep = ",")
names(true_bc_curves) <- c("pop", "model", "scen", "deg", paste("t", 1:31, sep = ""))
head(true_bc_curves)

agg_bc_curves <- aggregate(true_bc_curves[,5:35], 
                           list(true_bc_curves$model,
                                true_bc_curves$scen), mean)

names(agg_bc_curves)[names(agg_bc_curves) == "Group.2"] <- "scen"


# Convert data to wide format

medians <- reshape(full_curves[,c(1:6)], 
                   idvar = c("pop","mod_gen","scen","mod_fit"), 
                   direction = "wide")

lowers <- reshape(full_curves[,c(1:5,7)], 
                  idvar = c("pop","mod_gen","scen","mod_fit"), 
                  direction = "wide")

uppers <- reshape(full_curves[,c(1:5,8)], 
                  idvar = c("pop","mod_gen","scen","mod_fit"), 
                  direction = "wide")

# Obtain mean curve HPDs and merge with true curves

agg_medians <- aggregate(medians[,5:35], list(medians$mod_gen, medians$mod_fit, medians$scen), function(x) mean(x, na.rm = T))
names(agg_medians)[names(agg_medians) == "Group.3"] <- "scen"
agg_medians <- merge(agg_medians, agg_bc_curves, by = c("Group.1", "scen"))

agg_lowers <- aggregate(lowers[,5:35], list(lowers$mod_gen, lowers$mod_fit, lowers$scen), function(x) mean(x, na.rm = T))

agg_uppers <- aggregate(uppers[,5:35], list(uppers$mod_gen, uppers$mod_fit, uppers$scen), function(x) mean(x, na.rm = T))


par(mfrow = c(4,5))
#Model 1A, fit with Baseline SILM
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40" & agg_uppers$Group.2 == "base"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Baseline SILM", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40" & agg_lowers$Group.2 == "base"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"  & agg_lowers$Group.2 == "base"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "base"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "base"),4:34], lty = 2, col = "navy")

#Model 1A, fit with Model 1A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40" & agg_uppers$Group.2 == "mod1a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 1A*", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40" & agg_lowers$Group.2 == "mod1a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"  & agg_lowers$Group.2 == "mod1a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod1a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod1a"),4:34], lty = 2, col = "navy")

#Model 1A, fit with Model 2A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40" & agg_uppers$Group.2 == "mod2a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 2A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40" & agg_lowers$Group.2 == "mod2a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"  & agg_lowers$Group.2 == "mod2a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod2a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod2a"),4:34], lty = 2, col = "navy")

#Model 1A, fit with Model 3A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40" & agg_uppers$Group.2 == "mod3a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 3A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40" & agg_lowers$Group.2 == "mod3a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"  & agg_lowers$Group.2 == "mod3a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod3a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod3a"),4:34], lty = 2, col = "navy")


#Model 1A, fit with Model 4A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40" & agg_uppers$Group.2 == "mod4a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 4A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40" & agg_lowers$Group.2 == "mod4a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"  & agg_lowers$Group.2 == "mod4a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod4a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40" & agg_medians$Group.2 == "mod4a"),4:34], lty = 2, col = "navy")


#Model 2A, fit with Baseline SILM
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01" & agg_uppers$Group.2 == "base"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Baseline SILM", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01" & agg_lowers$Group.2 == "base"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"  & agg_lowers$Group.2 == "base"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "base"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "base"),4:34], lty = 2, col = "navy")

#Model 2A, fit with Model 1A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01" & agg_uppers$Group.2 == "mod1a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 1A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01" & agg_lowers$Group.2 == "mod1a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"  & agg_lowers$Group.2 == "mod1a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod1a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod1a"),4:34], lty = 2, col = "navy")

#Model 2A, fit with Model 2A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01" & agg_uppers$Group.2 == "mod2a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 2A*", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01" & agg_lowers$Group.2 == "mod2a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"  & agg_lowers$Group.2 == "mod2a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod2a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod2a"),4:34], lty = 2, col = "navy")

#Model 2A, fit with Model 3A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01" & agg_uppers$Group.2 == "mod3a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 3A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01" & agg_lowers$Group.2 == "mod3a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"  & agg_lowers$Group.2 == "mod3a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod3a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod3a"),4:34], lty = 2, col = "navy")

#Model 2A, fit with Model 4A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01" & agg_uppers$Group.2 == "mod4a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 4A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01" & agg_lowers$Group.2 == "mod4a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"  & agg_lowers$Group.2 == "mod4a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod4a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01" & agg_medians$Group.2 == "mod4a"),4:34], lty = 2, col = "navy")

#Model 3A, fit with Baseline SILM
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03" & agg_uppers$Group.2 == "base"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Baseline SILM", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03" & agg_lowers$Group.2 == "base"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"  & agg_lowers$Group.2 == "base"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "base"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "base"),4:34], lty = 2, col = "navy")

#Model 3A, fit with Model 1A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03" & agg_uppers$Group.2 == "mod1a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 1A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03" & agg_lowers$Group.2 == "mod1a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"  & agg_lowers$Group.2 == "mod1a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod1a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod1a"),4:34], lty = 2, col = "navy")

#Model 3A, fit with Model 2A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03" & agg_uppers$Group.2 == "mod2a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 2A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03" & agg_lowers$Group.2 == "mod2a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"  & agg_lowers$Group.2 == "mod2a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod2a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod2a"),4:34], lty = 2, col = "navy")

#Model 3A, fit with Model 3A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03" & agg_uppers$Group.2 == "mod3a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 3A*", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03" & agg_lowers$Group.2 == "mod3a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"  & agg_lowers$Group.2 == "mod3a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod3a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod3a"),4:34], lty = 2, col = "navy")

#Model 3A, fit with Model 4A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03" & agg_uppers$Group.2 == "mod4a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 4A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03" & agg_lowers$Group.2 == "mod4a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"  & agg_lowers$Group.2 == "mod4a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod4a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03" & agg_medians$Group.2 == "mod4a"),4:34], lty = 2, col = "navy")

#Model 4A, fit with Baseline SILM
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075" & agg_uppers$Group.2 == "base"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Baseline SILM", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075" & agg_lowers$Group.2 == "base"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"  & agg_lowers$Group.2 == "base"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "base"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "base"),4:34], lty = 2, col = "navy")

#Model 4A, fit with Model 1A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075" & agg_uppers$Group.2 == "mod1a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 1A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075" & agg_lowers$Group.2 == "mod1a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"  & agg_lowers$Group.2 == "mod1a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod1a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod1a"),4:34], lty = 2, col = "navy")

#Model 4A, fit with Model 2A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075" & agg_uppers$Group.2 == "mod2a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 2A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075" & agg_lowers$Group.2 == "mod2a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"  & agg_lowers$Group.2 == "mod2a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod2a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod2a"),4:34], lty = 2, col = "navy")

#Model 4A, fit with Model 3A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075" & agg_uppers$Group.2 == "mod3a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 3A", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075" & agg_lowers$Group.2 == "mod3a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"  & agg_lowers$Group.2 == "mod3a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod3a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod3a"),4:34], lty = 2, col = "navy")

#Model 4A, fit with Model 4A
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075" & agg_uppers$Group.2 == "mod4a"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 4A*", col = "white", cex.main = 1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075" & agg_lowers$Group.2 == "mod4a"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"  & agg_lowers$Group.2 == "mod4a"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod4a"),35:65], lwd = 2)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075" & agg_medians$Group.2 == "mod4a"),4:34], lty = 2, col = "navy")