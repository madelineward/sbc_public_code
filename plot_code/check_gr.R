full_gr <- read.table("output/gr_full.txt", sep = ",", fill = T)
head(full_gr)
no_converge <- full_gr[which(full_gr$V2 > 1.2),]
names(no_converge) <- c("param", "psrf", "psrfci", "mod_gen", "mod_fit", "scen", "pop")
#no_converge[!duplicated(no_converge[,4:7]),]

no_converge[order(no_converge[,5]),]

no_converge_b <- no_converge[grep("b",no_converge$mod_gen),]
no_converge_b[order(no_converge_b[,5]),]

non_converge_mods <- no_converge[!duplicated(no_converge[,4:7]), 4:7]

write.table(non_converge_mods,file = "output/non_converged_models.txt", row.names = F, col.names = F, sep = ",")



