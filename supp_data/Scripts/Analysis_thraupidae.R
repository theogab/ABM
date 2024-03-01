### READ Thraupidae data
lib <- c("ape", "bite", "coda")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

data <- read.csv("supp_data/Data/Individual beak data (Gaboriau).csv")
tree <- read.tree("supp_data/Data/Coerebinae.tre")

data <- data[data$Juv %in% c("", "Adult", "-"),]
data <- data[data$Jetz_taxo %in% gsub("_", " ", tree$tip.label),c(3, 9:12)]
data_brut <- list(mean = sapply(tree$tip.label, function(sp) sapply(data[data$Jetz_taxo == gsub("_", " ", sp),2:5], mean, na.rm = T)),
                     range = sapply(tree$tip.label, function(sp) sapply(data[data$Jetz_taxo == gsub("_", " ", sp),2:5], function(x) diff(range(x, na.rm = T)))),
                     mid = sapply(tree$tip.label, function(sp) sapply(data[data$Jetz_taxo == gsub("_", " ", sp),2:5], function(x) mean(range(x, na.rm = T)))),
                     count = sapply(tree$tip.label, function(sp) sapply(data[data$Jetz_taxo == gsub("_", " ", sp),2:5], function(x) sum(!is.na(x)))),
                     sd = sapply(tree$tip.label, function(sp) sapply(data[data$Jetz_taxo == gsub("_", " ", sp),2:5], function(x) sd(x, na.rm = T))))


####### MODEL TESTING ########
## JIVE ##
data_jive <- data
data_jive$Jetz_taxo <- gsub(" ", "_", data_jive$Jetz_taxo)
#total culmen
jive.BM <- make_jive(tree, data_jive[, 1:2], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_total_culmen_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
jive.OU <- make_jive(tree, data_jive[, 1:2], model.priors = list(mean = "BM", logvar = "OU"), scale = F)
mcmc_bite(jive.OU, log.file = "Application/Results/jive_M-BM_V-OU_total_culmen_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
#bill nares
jive.BM <- make_jive(tree, data_jive[, c(1,3)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_nares_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
jive.OU <- make_jive(tree, data_jive[, c(1,3)], model.priors = list(mean = "BM", logvar = "OU"), scale = F)
mcmc_bite(jive.OU, log.file = "Application/Results/jive_M-BM_V-OU_bill_nares_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
#bill width
jive.BM <- make_jive(tree, data_jive[, c(1,4)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_width_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
jive.OU <- make_jive(tree, data_jive[, c(1,4)], model.priors = list(mean = "BM", logvar = "OU"), scale = F)
mcmc_bite(jive.OU, log.file = "Application/Results/jive_M-BM_V-OU_bill_width_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
#bill depth
jive.BM <- make_jive(tree, data_jive[, c(1,5)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_depth_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)
jive.OU <- make_jive(tree, data_jive[, c(1,5)], model.priors = list(mean = "BM", logvar = "OU"), scale = F)
mcmc_bite(jive.OU, log.file = "Application/Results/jive_M-BM_V-OU_bill_depth_TI_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 10, ngen = 2000000)

#total culmen
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_total_culmen_TI_mcmc.log", header = TRUE, sep = "\t")
res_OU <- read.csv("Application/Results/jive_M-BM_V-OU_total_culmen_TI_mcmc.log", header = TRUE, sep = "\t")
mlik_BM <- marginal_lik(res_BM, burnin = 0.1, method = "SS")
mlik_OU <- marginal_lik(res_OU, burnin = 0.1, method = "SS")
#bill nares
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_nares_TI_mcmc.log", header = TRUE, sep = "\t")
res_OU <- read.csv("Application/Results/jive_M-BM_V-OU_bill_nares_TI_mcmc.log", header = TRUE, sep = "\t")
mlik_BM <- marginal_lik(res_BM, burnin = 0.1, method = "SS")
mlik_OU <- marginal_lik(res_OU, burnin = 0.1, method = "SS")
#bill width
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_width_TI_mcmc.log", header = TRUE, sep = "\t")
res_OU <- read.csv("Application/Results/jive_M-BM_V-OU_bill_width_TI_mcmc.log", header = TRUE, sep = "\t")
mlik_BM <- marginal_lik(res_BM, burnin = 0.1, method = "SS")
mlik_OU <- marginal_lik(res_OU, burnin = 0.1, method = "SS")
#bill depth
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_depth_TI_mcmc.log", header = TRUE, sep = "\t")
res_OU <- read.csv("Application/Results/jive_M-BM_V-OU_bill_depth_TI_mcmc.log", header = TRUE, sep = "\t")
mlik_BM <- marginal_lik(res_BM, burnin = 0.1, method = "SS")
mlik_OU <- marginal_lik(res_OU, burnin = 0.1, method = "SS")

## Parameter Estimation
#total culmen
jive.BM <- make_jive(tree, data_jive[, 1:2], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_total_culmen_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 1, ngen = 1000000)
#bill nares
jive.BM <- make_jive(tree, data_jive[, c(1,3)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_nares_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 1, ngen = 1000000)
#bill width
jive.BM <- make_jive(tree, data_jive[, c(1,4)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_width_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 1, ngen = 1000000)
#bill depth
jive.BM <- make_jive(tree, data_jive[, c(1,5)], model.priors = list(mean = "BM", logvar = "BM"), scale = F)
mcmc_bite(jive.BM, log.file = "Application/Results/jive_M-BM_V-BM_bill_depth_mcmc.log", sampling.freq = 1000, print.freq = 1000, ncat = 1, ngen = 1000000)

res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_total_culmen_mcmc.log", header = TRUE, sep = "\t")
colMeans(res_BM)
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_nares_mcmc.log", header = TRUE, sep = "\t")
colMeans(res_BM)
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_width_TI_mcmc.log", header = TRUE, sep = "\t")
colMeans(res_BM)
res_BM <- read.csv("Application/Results/jive_M-BM_V-BM_bill_depth_TI_mcmc.log", header = TRUE, sep = "\t")
colMeans(res_BM)


### ABM ###
ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[1,], data_brut$sd[1,]), scale = F, switch = T)
ADIM.obj$init$pars[1:2] <- 0.01
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_total_culmen_Ind1.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_total_culmen_Ind2.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
res_tot_cul1 <- read.csv("Application/Results/coerebinae_total_culmen_Ind1.log", header = T, sep = "\t")
res_tot_cul2 <- read.csv("Application/Results/coerebinae_total_culmen_Ind2.log", header = T, sep = "\t")
res_tot_cul1 <- res_tot_cul1[-(1:500),]
res_tot_cul2 <- res_tot_cul2[-(1:500),]
res_tot_cul_BF <- rbind(res_tot_cul1, res_tot_cul2)
BF_nu <- 2 * log((mean(res_tot_cul_BF$In)/(1 - mean(res_tot_cul_BF$In)))/((0.05)/(0.95)))
BF_omega <- 2 * log((mean(res_tot_cul_BF$Io)/(1 - mean(res_tot_cul_BF$Io)))/((0.05)/(0.95)))

ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[1,], data_brut$sd[1,]), scale = F, switch = F)
ADIM.obj$update.freq[c(3,7,8)] <- 0
ADIM.obj$init$pars[1:2] <- 0.01
ADIM.obj$init$I[1] <- 0
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_total_culmen_par.log", sampling.freq = 100, print.freq = 100, ngen = 500000)
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_total_culmen_par2.log", sampling.freq = 100, print.freq = 100, ngen = 500000)
res1 <- read.csv("Application/Results/coerebinae_total_culmen_par.log", header = T, sep = "\t")
res2 <- read.csv("Application/Results/coerebinae_total_culmen_par2.log", header = T, sep = "\t")
res1[,grepl("logv_", colnames(res1))] <- exp(res1[,grepl("logv_", colnames(res1))])
res2[,grepl("logv_", colnames(res2))] <- exp(res2[,grepl("logv_", colnames(res2))])
res_tot_cul <- as.mcmc(rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)]))
summary(res_tot_cul)
res_tot_cul <- rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)])

png("Application/Figures/coerebinae_total_culmen_BF.png", width = 10, height = 20, unit = "cm", res = 400)
par(mfrow = c(2,1), las = 1)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", nu)), xlab = expression(bar(nu)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_tot_cul_BF$nu*res_tot_cul_BF$In, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, sprintf("BF = %s", round(BF_nu, 2)), pos = 4)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", omega)), xlab = expression(bar(omega)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_tot_cul_BF$omega*res_tot_cul_BF$Io, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, sprintf("BF = %s", round(BF_omega, 2)), pos = 4)
dev.off()

png("Application/Figures/coerebinae_total_culmen.png", width = 20, height = 30, unit = "cm", res = 400)
par(mfrow = c(3,2), las = 1)
#hist(res_tot_cul$nu[res_tot_cul2$In == 1] , ylab = "Posterior density", xlab = expression(bar(nu)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_tot_cul$omega, ylab = "Posterior density", xlab = expression(bar(omega)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_tot_cul$sigma.2_m, ylab = "Posterior density", xlab = expression(bar(sigma[m]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_tot_cul$sigma.2_logv, ylab = "Posterior density", xlab = expression(bar(sigma[zeta]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_tot_cul$m_root, ylab = "Posterior density", xlab = expression(bar(m[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_tot_cul$logv_root, ylab = "Posterior density", xlab = expression(bar(zeta[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
dev.off()



#####  BILL NARES #####
ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[2,], data_brut$sd[2,]), scale = F, switch = T)
ADIM.obj$init$pars[1:2] <- 0.01
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_nares_Ind1.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_nares_Ind2.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
res_bill_nar1 <- read.csv("Application/Results/coerebinae_bill_nares_Ind1.log", header = T, sep = "\t")
res_bill_nar2 <- read.csv("Application/Results/coerebinae_bill_nares_Ind2.log", header = T, sep = "\t")
res_bill_nar1 <- res_bill_nar1[-(1:500),]
res_bill_nar2 <- res_bill_nar2[-(1:500),]
res_bill_nar_BF <- rbind(res_bill_nar1, res_bill_nar2)
BF_nu <- 2 * log((mean(res_bill_nar_BF$In)/(1 - mean(res_bill_nar_BF$In)))/((0.05)/(0.95)))
BF_omega <- 2 * log((mean(res_bill_nar_BF$Io)/(1 - mean(res_bill_nar_BF$Io)))/((0.05)/(0.95)))

ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[2,], data_brut$sd[2,]), scale = F, switch = F)
ADIM.obj$update.freq[c(3,7,8)] <- 0
ADIM.obj$init$pars[1:2] <- 0.01
ADIM.obj$init$I[1] <- 0
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_nares_par1.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_nares_par2.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
res1 <- read.csv("Application/Results/coerebinae_bill_nares_par1.log", header = T, sep = "\t")
res1[,grepl("logv_", colnames(res1))] <- exp(res1[,grepl("logv_", colnames(res1))])
res2 <- read.csv("Application/Results/coerebinae_bill_nares_par.log", header = T, sep = "\t")
res2[,grepl("logv_", colnames(res2))] <- exp(res2[,grepl("logv_", colnames(res2))])
res_bill_nar <- as.mcmc(rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)]))
summary(res_bill_nar)
res_bill_nar <- rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)])

png("Application/Figures/coerebinae_bill_nares_BF.png", width = 10, height = 20, unit = "cm", res = 400)
par(mfrow = c(2,1), las = 1)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", nu)), xlab = expression(bar(nu)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_nar$nu*res_bill_nar$In, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, sprintf("BF = %s", round(BF_nu, 2)), pos = 4)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", omega)), xlab = expression(bar(omega)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_nar$omega*res_bill_nar$Io, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, sprintf("BF = %s", round(BF_omega, 2)), pos = 4)
dev.off()

png("Application/Figures/coerebinae_billl_nares.png", width = 20, height = 30, unit = "cm", res = 400)
par(mfrow = c(3,2), las = 1)
#hist(res_bill_nar$nu[res_bill_nar2$In == 1] , ylab = "Posterior density", xlab = expression(bar(nu)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_nar$omega, ylab = "Posterior density", xlab = expression(bar(omega)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_nar$sigma.2_m, ylab = "Posterior density", xlab = expression(bar(sigma[m]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_nar$sigma.2_logv, ylab = "Posterior density", xlab = expression(bar(sigma[zeta]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_nar$m_root, ylab = "Posterior density", xlab = expression(bar(m[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_nar$logv_root, ylab = "Posterior density", xlab = expression(bar(zeta[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
dev.off()



##### BILL WIDTH #####
ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[3,], data_brut$sd[3,]), scale = F, switch = T)
ADIM.obj$init$pars[1:2] <- 0.01
mcmc_ADIM(ADIM.obj, log.file = "Application/Results/ADIM_coerebinae_bill_width1.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)
mcmc_ADIM(ADIM.obj, log.file = "Application/Results/ADIM_coerebinae_bill_width2.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)
res_bill_width1 <- read.csv("Application/Results/ADIM_coerebinae_bill_width1.log", header = T, sep = "\t")
res_bill_width2 <- read.csv("Application/Results/ADIM_coerebinae_bill_width2.log", header = T, sep = "\t")
res_bill_width1 <- res_bill_width1[-(1:500),]
res_bill_width2 <- res_bill_width2[-(1:500),]
res_bill_width_BF <- rbind(res_bill_width1, res_bill_width2)
BF_nu <- 2 * log((mean(res_bill_width_BF$In)/(1 - mean(res_bill_width_BF$In)))/((0.05)/(0.95)))
BF_omega <- 2 * log((mean(res_bill_width_BF$Io)/(1 - mean(res_bill_width_BF$Io)))/((0.05)/(0.95)))

ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[3,], data_brut$sd[3,]), scale = F, switch = F)
ADIM.obj$update.freq[c(3,4,7,8)] <- 0
ADIM.obj$init$pars[1:2] <- 0.01
ADIM.obj$init$I <- c(0,0)
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_width_par.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_width_par2.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
res1 <- read.csv("Application/Results/coerebinae_bill_width_par.log", header = T, sep = "\t")
res1[,grepl("logv_", colnames(res1))] <- exp(res1[,grepl("logv_", colnames(res1))])
res2 <- read.csv("Application/Results/coerebinae_bill_width_par2.log", header = T, sep = "\t")
res2[,grepl("logv_", colnames(res2))] <- exp(res2[,grepl("logv_", colnames(res2))])
res_bill_width <- as.mcmc(rbind(res1[-(1:500),-ncol(res1)],res2[-(1:500),-ncol(res2)]))
summary(res_bill_width)
res_bill_width <- rbind(res1[-(1:500),-ncol(res1)],res2[-(1:500),-ncol(res2)])

png("Application/Figures/coerebinae_bill_width_BF.png", width = 10, height = 20, unit = "cm", res = 400)
par(mfrow = c(2,1), las = 1)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", nu)), xlab = expression(bar(nu)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_width$nu*res_bill_width$In, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, "BF = -6.19", pos = 4)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", omega)), xlab = expression(bar(omega)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_width$omega*res_bill_width$Io, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, "BF = 1.10", pos = 4)
dev.off()

png("Application/Figures/coerebinae_bill_width.png", width = 20, height = 30, unit = "cm", res = 400)
par(mfrow = c(3,2), las = 1)
#hist(res_bill_width2$nu[res_bill_width2$In == 1] , ylab = "Posterior density", xlab = expression(bar(nu)), col =  adjustcolor("#3943b7", 0.8), main = "")
#hist(res_bill_width$omega, ylab = "Posterior density", xlab = expression(bar(omega)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_width$sigma.2_m, ylab = "Posterior density", xlab = expression(bar(sigma[m]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_width$sigma.2_logv, ylab = "Posterior density", xlab = expression(bar(sigma[zeta]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_width$m_root, ylab = "Posterior density", xlab = expression(bar(m[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_width$logv_root, ylab = "Posterior density", xlab = expression(bar(zeta[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
dev.off()


### BILL DEPTH ###
ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[4,], data_brut$sd[4,]), scale = F, switch = T)
ADIM.obj$init$pars[1:2] <- 0.01
mcmc_ADIM(ADIM.obj, log.file = "Application/Results/ADIM_coerebinae_bill_depth1.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
mcmc_ADIM(ADIM.obj, log.file = "Application/Results/ADIM_coerebinae_bill_depth2.log", sampling.freq = 100, print.freq = 100, ngen = 1000000)  
res_bill_depth1 <- read.csv("Application/Results/ADIM_coerebinae_bill_depth1.log", header = T, sep = "\t")
res_bill_depth2 <- read.csv("Application/Results/ADIM_coerebinae_bill_depth2.log", header = T, sep = "\t")
res_bill_depth1 <- res_bill_depth1[-(1:500),]
res_bill_depth2 <- res_bill_depth2[-(1:500),]
res_bill_depth_BF <- rbind(res_bill_depth1, res_bill_depth2)
BF_nu <- 2 * log((mean(res_bill_depth_BF$In)/(1 - mean(res_bill_depth_BF$In)))/((0.05)/(0.95)))
BF_omega <- 2 * log((mean(res_bill_depth_BF$Io)/(1 - mean(res_bill_depth_BF$Io)))/((0.05)/(0.95)))

ADIM.obj <- make_ADIM(phy = tree, traits = cbind(data_brut$mean[4,], data_brut$sd[4,]), scale = F, switch = F)
ADIM.obj$update.freq[c(3,7,8)] <- 0
ADIM.obj$init$pars[1:2] <- 0.01
ADIM.obj$init$I[1] <- 0
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_depth_par1.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
mcmc_ADIM(ADIM = ADIM.obj, log.file = "Application/Results/coerebinae_bill_depth_par2.log", sampling.freq = 100, print.freq = 100, ngen = 500000)  
res1 <- read.csv("Application/Results/coerebinae_bill_depth_par1.log", header = T, sep = "\t")
res1[,grepl("logv_", colnames(res1))] <- exp(res1[,grepl("logv_", colnames(res1))])# * var(data_nt[,5], na.rm = T) + mean(data_nt[,5], na.rm = T)
res2 <- read.csv("Application/Results/coerebinae_bill_depth_par2.log", header = T, sep = "\t")
res2[,grepl("logv_", colnames(res2))] <- exp(res2[,grepl("logv_", colnames(res2))])# * var(data_nt[,5], na.rm = T) + mean(data_nt[,5], na.rm = T)

res_bill_depth <- as.mcmc(rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)]))
summary(res_bill_depth)
res_bill_depth <- rbind(res1[-(1:500),-ncol(res1)], res2[-(1:500),-ncol(res2)])

png("Application/Figures/coerebinae_bill_depth_BF.png", width = 10, height = 20, unit = "cm", res = 400)
par(mfrow = c(2,1), las = 1)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", nu)), xlab = expression(bar(nu)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_depth_BF$nu*res_bill_depth_BF$In, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, "BF = -5.25", pos = 4)
prior <- hist(rbinom(1000,1,0.05)*runif(1000,0,.99), col =  adjustcolor("#3943b7", 0.6), main = expression(paste("Posterior odds for ", omega)), xlab = expression(bar(omega)), border = NA, breaks = 20, freq = F)
par(new = T)
hist(res_bill_depth_BF$omega*res_bill_depth_BF$Io, col =  adjustcolor("#4EA699", 0.6), main = "", xlab = "", border = NA, breaks = 20, axes = F, ylab = "", xlim = c(0,1), ylim = c(0,max(prior$density)), freq = F)
legend(0.6,14,c("prior", "posterior"), fill = adjustcolor(c("#3943b7","#4EA699"), 0.6), bty = "n")
text(0.6,15, "BF = 6.32*", pos = 4)
dev.off()

png("Application/Figures/coerebinae_bill_depth.png", width = 20, height = 30, unit = "cm", res = 400)
par(mfrow = c(3,2), las = 1)
#hist(res_bill_depth2$nu[res_bill_depth2$In == 1] , ylab = "Posterior density", xlab = expression(bar(nu)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_depth$omega, ylab = "Posterior density", xlab = expression(bar(omega)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_depth$sigma.2_m, ylab = "Posterior density", xlab = expression(bar(sigma[m]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_depth$sigma.2_logv, ylab = "Posterior density", xlab = expression(bar(sigma[zeta]^2)), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_depth$m_root, ylab = "Posterior density", xlab = expression(bar(m[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
hist(res_bill_depth$logv_root, ylab = "Posterior density", xlab = expression(bar(zeta[root])), col =  adjustcolor("#3943b7", 0.8), main = "")
dev.off()

####################################################################################
############################     PLOTS    ##########################################
####################################################################################


png("Application/Figures/coerebinae_hist_traits.png", width = 20, height = 24, unit = "cm", res = 400)
par(mfrow = c(2,2), las = 1)
hist(data_nt$Bill_TotalCulmen, main = "Total culmen (mm)", xlab = "Total culmen", col = adjustcolor("#2E86AB", .8))
hist(data_nt$Bill_Nares, main = "Beak Nares (mm)", xlab = "Beak Nares", col = adjustcolor("#BDD358", .8))
hist(data_nt$Bill_Width, main = "Beak width (mm)", xlab = "Beak width", col = adjustcolor("#484041", .8))
hist(data_nt$Bill_Depth, main = "Beak depth (mm)", xlab = "Beak depth", col = adjustcolor("#EA9010", .8))
dev.off()

source("~/Dropbox/Documents/R/Graphs/add_geoscale.R")
plot.pvo <- function(n, lim, trait, lab, color){
  rge <- seq(lim[1], lim[2], length = 1e4)
  plot(md$yy[1:n]~data_brut$mean[2,], yaxt = "n", ylab = "", xlab = lab, pch = 16, col = "#FFFFFF", xlim = lim, bty = "n")
  for(i in 1:n){
    polygon(c(50,0,0,50), c(i-0.2,i-0.2,i+0.8,i+0.8), border = NA, col= adjustcolor(c("lightgrey", "darkgrey")[i %% 2 + 1], .4), xpd = T)
  }
  for(i in 1:n){
    d <- dnorm(rge, data_brut$mean[trait,tree$tip.label[i]], data_brut$sd[trait,tree$tip.label[i]])
    d <- d/max(d)*0.8
    polygon(rge, md$yy[i]-0.2+d , border = NA, col = color, xpd = T)
    
    pt <- data_nt[data_nt[,1] == gsub("_", " ", tree$tip.label[i]),trait + 1]
    pt <- pt[!is.na(pt)]
    points(pt,rep(md$yy[i]-0.1, length(pt)), pch = 16, cex = 0.8)
    for(j in 1:length(pt)){
      lines(rep(pt[j],2), md$yy[i] - c(0.1, 0.2), lwd = 0.5)
    }
  } 
}

pdf("Application/Figures/coerebinae_phy_traits.pdf", width = 20, height = 12.5)
layout(matrix(c(1,1,2,3,4,5,6), nrow = 1))
par(mar = c(5.1, 1.1, 4.1, 0))
plot(tree,label.offset = .1, x.lim = c(0,12.2), cex = 1.2, plot = F)
for(i in 1:14){
  polygon(c(6.38,20,20,6.38), c(i-0.2,i-0.2,i+0.8,i+0.8), border = NA, col= adjustcolor(c("lightgrey", "darkgrey")[i %% 2 + 1], .4), xpd = T)
}
par(new = T)
plot(tree,label.offset = .1, x.lim = c(0,12.2), cex = 1.2)
axisGeoX(c(0,6.28), ymax = 0, posbar = 0.5, unit = 'epoch', cex = 1.2, ages = T, posage = -0.7,
         col = c("#fef2e0","#fff2a2","#ffff99","#ffff00", "#fdc07a","#fdb46c","#fda75f"))
usr <- par()$usr[1:2]
md <- get('last_plot.phylo', envir = .PlotPhyloEnv)
par(mar = c(5.1, 0, 4.1, 0), cex.lab = 1.2)
plot.pvo(14, c(4,28), 1,  "Total culmen (mm)", adjustcolor("#2E86AB", .8))
plot.pvo(14, c(2,26), 2, "Beak nares (mm)", adjustcolor("#BDD358", .8))
plot.pvo(14, c(2,18), 3, "Beak width (mm)", adjustcolor("#484041", .8))
plot.pvo(14,  c(2,22), 4, "Beak depth (mm)", adjustcolor("#EA9010", .8))

dev.off()

plot.pvo <- function(n, lim, trait, lab, color){
  rge <- seq(lim[1], lim[2], length = 1e4)
  plot(md$yy[1:n]~data_grouped$mean[2,], yaxt = "n", ylab = "", xlab = lab, pch = 16, col = "#FFFFFF", xlim = lim, bty = "n")
  for(i in 1:n){
    polygon(c(50,-50,-50,50), c(i-0.2,i-0.2,i+0.8,i+0.8), border = NA, col= adjustcolor(c("lightgrey", "darkgrey")[i %% 2 + 1], .4), xpd = T)
  }
  for(i in 1:n){
    d <- dnorm(rge, data_grouped$mean[trait,tree$tip.label[i]], data_grouped$sd[trait,tree$tip.label[i]])
    d <- d/max(d)*0.8
    polygon(rge, md$yy[i]-0.2+d , border = NA, col = color, xpd = T)
    
    pt <- data[data[,1] == gsub("_", " ", tree$tip.label[i]),trait + 1]
    pt <- pt[!is.na(pt)]
    points(pt,rep(md$yy[i]-0.1, length(pt)), pch = 16, cex = 0.8)
    for(j in 1:length(pt)){
      lines(rep(pt[j],2), md$yy[i] - c(0.1, 0.2), lwd = 0.5)
    }
  } 
}

pdf("Application/Figures/coerebinae_phy_traits_std.pdf", width = 20, height = 12.5)
layout(matrix(c(1,1,2,3,4,5,6), nrow = 1))
par(mar = c(5.1, 1.1, 4.1, 0))
plot(tree,label.offset = .1, x.lim = c(0,12.2), cex = 1.2, plot = F)
for(i in 1:14){
  polygon(c(6.38,20,20,6.38), c(i-0.2,i-0.2,i+0.8,i+0.8), border = NA, col= adjustcolor(c("lightgrey", "darkgrey")[i %% 2 + 1], .4), xpd = T)
}
par(new = T)
plot(tree,label.offset = .1, x.lim = c(0,12.2), cex = 1.2)
axisGeoX(c(0,6.28), ymax = 0, posbar = 0.5, unit = 'epoch', cex = 1.2, ages = T, posage = -0.7,
         col = c("#fef2e0","#fff2a2","#ffff99","#ffff00", "#fdc07a","#fdb46c","#fda75f"))
usr <- par()$usr[1:2]
md <- get('last_plot.phylo', envir = .PlotPhyloEnv)
par(mar = c(5.1, 0, 4.1, 0), cex.lab = 1.2)
plot.pvo(14, c(-0.6,1.2), 1,  "Total culmen (mm)", adjustcolor("#2E86AB", .8))
plot.pvo(14, c(-1.2,2.2), 2, "Beak nares (mm)", adjustcolor("#BDD358", .8))
plot.pvo(14, c(-2,3), 3, "Beak width (mm)", adjustcolor("#484041", .8))
plot.pvo(14, c(-1,2), 4, "Beak depth (mm)", adjustcolor("#EA9010", .8))

dev.off()