# Figure 3 

tmax <- 30

full_curves <- read.table("output/curves_full.txt", sep = ",", fill = T)
names(full_curves) <- c("time", "median", "lower", "upper", "pop", "mod_gen", "scen", "mod_fit")
head(full_curves)
full_curves <- full_curves[-61474,]
full_curves[2:4] <- sapply(full_curves[2:4],as.numeric)

# Only need to keep curves where true model was fit
full_curves <- full_curves[which(full_curves$mod_gen == full_curves$mod_fit),]

true_bc_curves <- read.table("data/true_curves.txt", sep = ",")
names(true_bc_curves) <- c("pop", "model", "scen", "deg", paste("t", 1:31, sep = ""))
head(true_bc_curves)

agg_bc_curves <- aggregate(true_bc_curves[,5:35], 
                           list(true_bc_curves$model,
                                true_bc_curves$scen), mean)

names(agg_bc_curves)[names(agg_bc_curves) == "Group.2"] <- "scen"

# Convert data to wide format

medians <- reshape(full_curves[,c(1,2,5:8)], 
                   idvar = c("pop","mod_gen","scen","mod_fit"), 
                   direction = "wide")

lowers <- reshape(full_curves[,c(1,3,5:8)], 
                   idvar = c("pop","mod_gen","scen","mod_fit"), 
                   direction = "wide")

uppers <- reshape(full_curves[,c(1,4,5:8)], 
                   idvar = c("pop","mod_gen","scen","mod_fit"), 
                   direction = "wide")

# Obtain mean curve HPDs and merge with true curves

agg_medians <- aggregate(medians[,5:35], list(medians$mod_gen, medians$mod_fit, medians$scen), mean)
names(agg_medians)[names(agg_medians) == "Group.3"] <- "scen"
agg_medians <- merge(agg_medians, agg_bc_curves, by = c("Group.1", "scen"))

agg_lowers <- aggregate(lowers[,5:35], list(lowers$mod_gen, lowers$mod_fit, lowers$scen), mean)

agg_uppers <- aggregate(uppers[,5:35], list(uppers$mod_gen, uppers$mod_fit, uppers$scen), mean)

# Plot curves for every simulation scenario
# par(mfrow = c(2,3))
# for(i in 1:nrow(agg_medians)) {
#   plot(1:(tmax+1), agg_uppers[i,4:34], type = "l",
#        xlab = "Time", ylab = "Incidence",
#        main = agg_uppers[i,"scen"], col = "white")
#   polygon(c(1:(tmax+1), (tmax+1):1),
#           c(agg_lowers[i,4:34], rev(agg_uppers[i,4:34])),
#           col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
#   lines(1:(tmax+1), agg_medians[i,35:65], lwd = 2)
#   lines(1:(tmax+1), agg_medians[i,4:34], lty = 2, col = "navy")
# }


par(mfrow = c(2,2))
#Model 1A, medium BC
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 1A", col = "white", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod1a" & agg_lowers$Group.3 == "a22b2d65h40"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod1a" & agg_uppers$Group.3 == "a22b2d65h40"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40"),35:65], lwd = 3)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod1a" & agg_medians$scen == "a22b2d65h40"),4:34], lty = 2, lwd = 2, col = "navy")

#Model 2A, medium BC
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 2A", col = "white", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod2a" & agg_lowers$Group.3 == "a24b2d01"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod2a" & agg_uppers$Group.3 == "a24b2d01"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01"),35:65], lwd = 3)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod2a" & agg_medians$scen == "a24b2d01"),4:34], lty = 2, lwd = 2, col = "navy")

#Model 3A, medium BC
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 3A", col = "white", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod3a" & agg_lowers$Group.3 == "a24b2d80d03"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod3a" & agg_uppers$Group.3 == "a24b2d80d03"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03"),35:65], lwd = 3)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod3a" & agg_medians$scen == "a24b2d80d03"),4:34], lty = 2, lwd = 2, col = "navy")

#Model 4A, medium BC
plot(1:(tmax+1), agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"),4:34], type = "l",
     xlab = "Time", ylab = "Incidence",
     main = "Model 4A", col = "white", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
polygon(c(1:(tmax+1), (tmax+1):1),
        c(agg_lowers[which(agg_lowers$Group.1 == "mod4a" & agg_lowers$Group.3 == "a24b2d3d075"),4:34], rev(agg_uppers[which(agg_uppers$Group.1 == "mod4a" & agg_uppers$Group.3 == "a24b2d3d075"),4:34])),
        col = adjustcolor('lightblue',alpha.f=0.5), lty = 0)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075"),35:65], lwd = 3)
lines(1:(tmax+1), agg_medians[which(agg_medians$Group.1 == "mod4a" & agg_medians$scen == "a24b2d3d075"),4:34], lty = 2, lwd = 2, col = "navy")