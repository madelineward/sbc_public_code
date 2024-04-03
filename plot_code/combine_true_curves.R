## This file makes a new data set containing the true epidemic curves for each population, model, and simulation scenario

num_pop = 20
tmax <- 30

# Function that calculates epidemic curve
count_inf <- function(lattimes, tmin, tmax) {
  x = tmin:(tmax+1)                        # Time points of the epidemic
  counts = rep(0,length(x))             # Placeholder vector of correct length
  for(i in 1:length(x)) {               # For each timepoint, find the number of new infections
    counts[i] = sum(lattimes %in% x[i])  
  }
  counts
}

mod_list <- c("mod1a", "mod1b", "mod2a", "mod2b", "mod3a", "mod3b", "mod4a", "mod4b")
all_mods <- rep(mod_list, each = 3)
scens <- c("a22b2d50h40", "a22b2d65h40", "a22b2d80h40",
           "a22b2d10h40","a22b2d15h40", "a22b2d20h40", 
           "a24b2d005","a24b2d01", "a24b2d015", 
           "a24b2d001", "a24b2d0015", "a24b2d002",
           "a24b2d80d02","a24b2d80d03", "a24b2d80d04",
           "a24b2d04d005", "a24b2d04d007", "a24b2d04d009", 
           "a24b2d3d05",  "a24b2d3d075", "a24b2d3d10",
           "a24b2d3d15", "a24b2d3d20", "a24b2d3d10")
strength_scen <- c(rep(c("weak", "med", "str"),6), c("str", "med", "weak", "med", "weak", "str"))

true_curves <- matrix(nrow = num_pop*length(scens), ncol = 31)
count <- 1
for(i in seq_along(all_mods)) {
  file_name <- paste("data/", all_mods[i], "_", scens[i], ".txt", sep = "")
  epi_df <- read.table(file_name, sep = ",", header = T)
  for(j in 1:num_pop) {
    prevs <- epi_df[which(epi_df$pop == j),"inftime"]
    true_curves[count,] <- count_inf(prevs, 1, 30)
    count <- count + 1
  }
}
true_curves <- as.data.frame(true_curves)
true_curves <- data.frame("pop" = rep(1:20, length(scens)), "model" = rep(all_mods, each = 20), "scen" = rep(scens, each = 20), "deg" = rep(strength_scen, each = 20), true_curves)

write.table(true_curves,file = "data/true_curves.txt", row.names = F, col.names = F, sep = ",")


## True curves, no BC effect 

mod_list <- c("mod1a", "mod1b", "mod2a", "mod2b", "mod3a", "mod3b", "mod4a", "mod4b")
scens <- c("a22b2d0h0", "a22b2d0h0", 
           "a24b2d0","a24b2d0", 
           "a24b2d0d0","a24b2d0d0", 
           "a24b2d0d0", "a24b2d0d0")

true_curves <- matrix(nrow = 20*length(scens), ncol = 31)
count <- 1
for(i in seq_along(mod_list)) {
  file_name <- paste("data/", mod_list[i], "_", scens[i], ".txt", sep = "")
  epi_df <- read.table(file_name, sep = ",", header = T)
  for(j in 1:num_pop) {
    prevs <- epi_df[which(epi_df$pop == j),"inftime"]
    true_curves[count,] <- count_inf(prevs, 1, 30)
    count <- count + 1
  }
}
true_curves <- as.data.frame(true_curves)
true_curves <- data.frame("pop" = rep(1:20, length(scens)), "model" = rep(mod_list, each = 20), true_curves)

write.table(true_curves,file = "data/true_curves_nobc.txt", row.names = F, col.names = F, sep = ",")
