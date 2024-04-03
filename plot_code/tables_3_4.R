# Tables 3 and 4

diagnos <- read.table("output/waic_and_dic_full.txt", header = F, sep = ",", fill = T)
head(diagnos)
names(diagnos) <- c("waic", "dic", "mod_gen","mod_fit", "scen", "pop")

# Remove models/populations that didn't converge
not_converged <- read.table("output/non_converged_models.txt", sep = ",")
names(not_converged) <- c("mod_gen", "mod_fit", "scen", "pop")
not_converged$converge <- 1

diagnos <- merge(diagnos, not_converged, by = c("mod_gen", "mod_fit", "scen", "pop"), all.x = T)
head(diagnos)
diagnos <- diagnos[-which(diagnos$converge == 1),1:6]

# Remove duplicates (ran some pops twice and both converged)

diagnos <- diagnos[!duplicated(diagnos[,1:4]),]

library(xtable)


### MODEL 1A
mod_1ai <- diagnos[which(diagnos$mod_gen == "mod1a" & 
                           diagnos$scen == "a22b2d50h40"),]
waic_1ai <- xtabs(waic ~ pop + mod_fit, mod_1ai)

mod_1aii <- diagnos[which(diagnos$mod_gen == "mod1a" & 
                            diagnos$scen == "a22b2d65h40"),]
waic_1aii <- xtabs(waic ~ pop + mod_fit, mod_1aii)

mod_1aiii <- diagnos[which(diagnos$mod_gen == "mod1a" & 
                             diagnos$scen == "a22b2d80h40"),]
waic_1aiii <- xtabs(waic ~ pop + mod_fit, mod_1aiii)


### MODEL 1B
mod_1bi <- diagnos[which(diagnos$mod_gen == "mod1b" & 
                           diagnos$scen == "a22b2d10h40"),]
waic_1bi <- xtabs(waic ~ pop + mod_fit, mod_1bi)

mod_1bii <- diagnos[which(diagnos$mod_gen == "mod1b" & 
                            diagnos$scen == "a22b2d15h40"),]
waic_1bii <- xtabs(waic ~ pop + mod_fit, mod_1bii)

mod_1biii <- diagnos[which(diagnos$mod_gen == "mod1b" & 
                             diagnos$scen == "a22b2d20h40"),]
waic_1biii <- xtabs(waic ~ pop + mod_fit, mod_1biii)


### MODEL 2A
mod_2ai <- diagnos[which(diagnos$mod_gen == "mod2a" & 
                           diagnos$scen == "a24b2d005"),]
waic_2ai <- xtabs(waic ~ pop + mod_fit, mod_2ai)

mod_2aii <- diagnos[which(diagnos$mod_gen == "mod2a" & 
                            diagnos$scen == "a24b2d01"),]
waic_2aii <- xtabs(waic ~ pop + mod_fit, mod_2aii)

mod_2aiii <- diagnos[which(diagnos$mod_gen == "mod2a" & 
                             diagnos$scen == "a24b2d015"),]
waic_2aiii <- xtabs(waic ~ pop + mod_fit, mod_2aiii)



### MODEL 2B

mod_2bi <- diagnos[which(diagnos$mod_gen == "mod2b" & 
                           diagnos$scen == "a24b2d001"),]
waic_2bi <- xtabs(waic ~ pop + mod_fit, mod_2bi)
dic_2bi <- xtabs(dic ~ pop + mod_fit, mod_2bi)

mod_2bii <- diagnos[which(diagnos$mod_gen == "mod2b" & 
                            diagnos$scen == "a24b2d0015"),]
waic_2bii <- xtabs(waic ~ pop + mod_fit, mod_2bii)
dic_2bii <- xtabs(dic ~ pop + mod_fit, mod_2bii)

mod_2biii <- diagnos[which(diagnos$mod_gen == "mod2b" & 
                             diagnos$scen == "a24b2d002"),]
waic_2biii <- xtabs(waic ~ pop + mod_fit, mod_2biii)
dic_2biii <- xtabs(dic ~ pop + mod_fit, mod_2biii)


### MODEL 3A
mod_3ai <- diagnos[which(diagnos$mod_gen == "mod3a" & 
                           diagnos$scen == "a24b2d80d02"),]
waic_3ai <- xtabs(waic ~ pop + mod_fit, mod_3ai)
dic_3ai <- xtabs(dic ~ pop + mod_fit, mod_3ai)

mod_3aii <- diagnos[which(diagnos$mod_gen == "mod3a" & 
                            diagnos$scen == "a24b2d80d03"),]
waic_3aii <- xtabs(waic ~ pop + mod_fit, mod_3aii)
dic_3aii <- xtabs(dic ~ pop + mod_fit, mod_3aii)

mod_3aiii <- diagnos[which(diagnos$mod_gen == "mod3a" & 
                             diagnos$scen == "a24b2d80d04"),]
waic_3aiii <- xtabs(waic ~ pop + mod_fit, mod_3aiii)
dic_3aiii <- xtabs(dic ~ pop + mod_fit, mod_3aiii)


### MODEL 3B
mod_3bi <- diagnos[which(diagnos$mod_gen == "mod3b" & 
                           diagnos$scen == "a24b2d04d005"),]
waic_3bi <- xtabs(waic ~ pop + mod_fit, mod_3bi)
dic_3bi <- xtabs(dic ~ pop + mod_fit, mod_3bi)

mod_3bii <- diagnos[which(diagnos$mod_gen == "mod3b" & 
                            diagnos$scen == "a24b2d04d007"),]
waic_3bii <- xtabs(waic ~ pop + mod_fit, mod_3bii)
dic_3bii <- xtabs(dic ~ pop + mod_fit, mod_3bii)

mod_3biii <- diagnos[which(diagnos$mod_gen == "mod3b" & 
                             diagnos$scen == "a24b2d04d009"),]
waic_3biii <- xtabs(waic ~ pop + mod_fit, mod_3biii)
dic_3biii <- xtabs(dic ~ pop + mod_fit, mod_3biii)


### MODEL 4A
mod_4ai <- diagnos[which(diagnos$mod_gen == "mod4a" & 
                           diagnos$scen == "a24b2d3d05"),]
waic_4ai <- xtabs(waic ~ pop + mod_fit, mod_4ai)
dic_4ai <- xtabs(dic ~ pop + mod_fit, mod_4ai)

mod_4aii <- diagnos[which(diagnos$mod_gen == "mod4a" & 
                            diagnos$scen == "a24b2d3d075"),]
waic_4aii <- xtabs(waic ~ pop + mod_fit, mod_4aii)
dic_4aii <- xtabs(dic ~ pop + mod_fit, mod_4aii)

mod_4aiii <- diagnos[which(diagnos$mod_gen == "mod4a" & 
                             diagnos$scen == "a24b2d3d10"),]
waic_4aiii <- xtabs(waic ~ pop + mod_fit, mod_4aiii)
dic_4aiii <- xtabs(dic ~ pop + mod_fit, mod_4aiii)





### MODEL 4B
mod_4bi <- diagnos[which(diagnos$mod_gen == "mod4b" & 
                           diagnos$scen == "a24b2d3d10"),]
waic_4bi <- xtabs(waic ~ pop + mod_fit, mod_4bi)
dic_4bi <- xtabs(dic ~ pop + mod_fit, mod_4bi)

mod_4bii <- diagnos[which(diagnos$mod_gen == "mod4b" & 
                            diagnos$scen == "a24b2d3d15"),]
waic_4bii <- xtabs(waic ~ pop + mod_fit, mod_4bii)
dic_4bii <- xtabs(dic ~ pop + mod_fit, mod_4bii)

mod_4biii <- diagnos[which(diagnos$mod_gen == "mod4b" & 
                             diagnos$scen == "a24b2d3d20"),]
waic_4biii <- xtabs(waic ~ pop + mod_fit, mod_4biii)
dic_4biii <- xtabs(dic ~ pop + mod_fit, mod_4biii)



#### Make WAIC Tables

mod_names <-c("Base", "1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B")

# Proportion of time model is selected 
min_g0 <- function(x) {
  return(min(x[which(x>0)]))
}
mean_n0 <- function(x) {
  remove_bad <- function(y) {
    z <- y[which(y < 4000 & y > -4000)]
    return(mean(z))
  }
  return(apply(x, 2, remove_bad))
}
prop_select <- function(x) {
  y <- apply(x, 1, min_g0)
  z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  z[which(x == y)] = 1
  return(z)
}

waic_weak <- rbind(colMeans(prop_select(waic_1ai)),
                   colMeans(prop_select(waic_1bi)),
                   colMeans(prop_select(waic_2ai)),
                   colMeans(prop_select(waic_2bi)),
                   colMeans(prop_select(waic_3ai)),
                   colMeans(prop_select(waic_3bi)),
                   colMeans(prop_select(waic_4ai)),
                   colMeans(prop_select(waic_4bi))
                   )
rownames(waic_weak) <- mod_names[-1]
colnames(waic_weak) <- mod_names
waic_weak

waic_med <- rbind(colMeans(prop_select(waic_1aii)),
                  colMeans(prop_select(waic_1bii)),
                  colMeans(prop_select(waic_2aii)),
                  colMeans(prop_select(waic_2bii)),
                  colMeans(prop_select(waic_3aii)),
                  colMeans(prop_select(waic_3bii)),
                  colMeans(prop_select(waic_4aii)),
                  colMeans(prop_select(waic_4bii)))
rownames(waic_med) <- mod_names[-1]
colnames(waic_med) <- mod_names
xtable(waic_med, type = "latex")

waic_str <- rbind(colMeans(prop_select(waic_1aiii)),
                  colMeans(prop_select(waic_1biii)),
                  colMeans(prop_select(waic_2aiii)),
                  colMeans(prop_select(waic_2biii)),
                  colMeans(prop_select(waic_3aiii)),
                  colMeans(prop_select(waic_3biii)),
                  colMeans(prop_select(waic_4aiii)),
                  colMeans(prop_select(waic_4biii)))
rownames(waic_str) <- mod_names[-1]
colnames(waic_str) <- mod_names
waic_str

# Mean difference from true model
wdif_weak <- rbind(round(mean_n0(apply(waic_1ai,2, function(x) x - waic_1ai[,2])),2),
                   round(mean_n0(apply(waic_1bi,2, function(x) x - waic_1bi[,3])),2),
                   round(mean_n0(apply(waic_2ai,2, function(x) x - waic_2ai[,4])),2),
                   round(mean_n0(apply(waic_2bi,2, function(x) x - waic_2bi[,5])),2),
                   round(mean_n0(apply(waic_3ai,2, function(x) x - waic_3ai[,6])),2),
                   round(mean_n0(apply(waic_3bi,2, function(x) x - waic_3bi[,7])),2),
                   round(mean_n0(apply(waic_4ai,2, function(x) x - waic_4ai[,8])),2),
                   round(mean_n0(apply(waic_4bi,2, function(x) x - waic_4bi[,9])),2))
rownames(wdif_weak) <- mod_names[-1]
colnames(wdif_weak) <- mod_names
xtable(wdif_weak, type = "latex")

wdif_med <- rbind(round(mean_n0(apply(waic_1aii,2, function(x) x - waic_1aii[,2])),2),
                  round(mean_n0(apply(waic_1bii,2, function(x) x - waic_1bii[,3])),2),
                  round(mean_n0(apply(waic_2aii,2, function(x) x - waic_2aii[,4])),2),
                  round(mean_n0(apply(waic_2bii,2, function(x) x - waic_2bii[,5])),2),
                  round(mean_n0(apply(waic_3aii,2, function(x) x - waic_3aii[,6])),2),
                  round(mean_n0(apply(waic_3bii,2, function(x) x - waic_3bii[,7])),2),
                  round(mean_n0(apply(waic_4aii,2, function(x) x - waic_4aii[,8])),2),
                  round(mean_n0(apply(waic_4bii,2, function(x) x - waic_4bii[,9])),2))
rownames(wdif_med) <- mod_names[-1]
colnames(wdif_med) <- mod_names
xtable(wdif_med, type = "latex")

rbind(round(mean_n0(apply(waic_1aiii,2, function(x) x - waic_1aiii[,2])),2),
      round(mean_n0(apply(waic_1biii,2, function(x) x - waic_1biii[,3])),2),
      round(mean_n0(apply(waic_2aiii,2, function(x) x - waic_2aiii[,4])),2),
      round(mean_n0(apply(waic_2biii,2, function(x) x - waic_2biii[,5])),2),
      round(mean_n0(apply(waic_3aiii,2, function(x) x - waic_3aiii[,6])),2),
      round(mean_n0(apply(waic_3biii,2, function(x) x - waic_3biii[,7])),2),
      round(mean_n0(apply(waic_4aiii,2, function(x) x - waic_4aiii[,8])),2),
      round(mean_n0(apply(waic_4biii,2, function(x) x - waic_4biii[,9])),2))

#### Make DIC Tables

# Proportion of time model is selected 

rbind(colMeans(prop_select(dic_1ai)),
      colMeans(prop_select(dic_1bi)),
      colMeans(prop_select(dic_2ai)),
      colMeans(prop_select(dic_2bi)),
      colMeans(prop_select(dic_3ai)),
      colMeans(prop_select(dic_3bi)),
      colMeans(prop_select(dic_4ai)),
      colMeans(prop_select(dic_4bi)))

rbind(colMeans(prop_select(dic_1aii)),
      colMeans(prop_select(dic_1bii)),
      colMeans(prop_select(dic_2aii)),
      colMeans(prop_select(dic_2bii)),
      colMeans(prop_select(dic_3aii)),
      colMeans(prop_select(dic_3bii)),
      colMeans(prop_select(dic_4aii)),
      colMeans(prop_select(dic_4bii)))

rbind(colMeans(prop_select(dic_1aiii)),
      colMeans(prop_select(dic_1biii)),
      colMeans(prop_select(dic_2aiii)),
      colMeans(prop_select(dic_2biii)),
      colMeans(prop_select(dic_3aiii)),
      colMeans(prop_select(dic_3biii)),
      colMeans(prop_select(dic_4aiii)),
      colMeans(prop_select(dic_4biii)))

# Mean difference from true model
rbind(round(colMeans(apply(dic_1ai,2, function(x) x - dic_1ai[,2])),2),
      round(colMeans(apply(dic_1bi,2, function(x) x - dic_1bi[,3])),2),
      round(colMeans(apply(dic_2ai,2, function(x) x - dic_2ai[,4])),2),
      round(colMeans(apply(dic_2bi,2, function(x) x - dic_2bi[,5])),2),
      round(colMeans(apply(dic_3ai,2, function(x) x - dic_3ai[,6])),2),
      round(colMeans(apply(dic_3bi,2, function(x) x - dic_3bi[,7])),2),
      round(colMeans(apply(dic_4ai,2, function(x) x - dic_4ai[,8])),2),
      round(colMeans(apply(dic_4bi,2, function(x) x - dic_4bi[,9])),2))

rbind(round(colMeans(apply(dic_1aii,2, function(x) x - dic_1aii[,2])),2),
      round(colMeans(apply(dic_1bii,2, function(x) x - dic_1bii[,3])),2),
      round(colMeans(apply(dic_2aii,2, function(x) x - dic_2aii[,4])),2),
      round(colMeans(apply(dic_2bii,2, function(x) x - dic_2bii[,5])),2),
      round(colMeans(apply(dic_3aii,2, function(x) x - dic_3aii[,6])),2),
      round(colMeans(apply(dic_3bii,2, function(x) x - dic_3bii[,7])),2),
      round(colMeans(apply(dic_4aii,2, function(x) x - dic_4aii[,8])),2),
      round(colMeans(apply(dic_4bii,2, function(x) x - dic_4bii[,9])),2))

rbind(round(colMeans(apply(dic_1aiii,2, function(x) x - dic_1aiii[,2])),2),
      round(colMeans(apply(dic_1biii,2, function(x) x - dic_1biii[,3])),2),
      round(colMeans(apply(dic_2aiii,2, function(x) x - dic_2aiii[,4])),2),
      round(colMeans(apply(dic_2biii,2, function(x) x - dic_2biii[,5])),2),
      round(colMeans(apply(dic_3aiii,2, function(x) x - dic_3aiii[,6])),2),
      round(colMeans(apply(dic_3biii,2, function(x) x - dic_3biii[,7])),2),
      round(colMeans(apply(dic_4aiii,2, function(x) x - dic_4aiii[,8])),2),
      round(colMeans(apply(dic_4biii,2, function(x) x - dic_4biii[,9])),2))
