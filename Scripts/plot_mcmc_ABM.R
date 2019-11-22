### ABM fit interpret results
lib <- c("ape", "phytools", "coda")
sapply(lib, library, character.only = T)

## input of the initial simulation
rho <- 0.1
delta <- 0.1
sig2 <- 0.002
load(sprintf("Results/sim_tree_rho-%s_delta-%s_sig-%s.Rdata", rho, delta, sig2))
eta <- rho/(rho*(delta-1)+1)
if(eta > 1) eta <- delta-1+1/rho
omega <- 1-delta*rho
set.seed(300)
phy <- pbtree(n = 10, scale = 100)
## all the simulated values ordered as in log.mcmc
sim <- c(NA, NA, sim_tree[[2]][(length(phy$tip.label)+1):(length(phy$tip.label) + phy$Nnode)], sig2, rho, delta, NA, NA, eta, omega)

## Results of the mcmc
simn <- 2
tree <- 1
log.mcmc <- read.csv(sprintf("Results/ABM_mcmc_rho-%s_delta-%s_sig-%s_tree-%s_simn-%s.log", rho, delta, sig2, tree, simn), sep = "\t")
log.mcmc <- log.mcmc[-ncol(log.mcmc)]
names(sim) <- colnames(log.mcmc)
eta <- log.mcmc$rho/(log.mcmc$rho*(log.mcmc$delta-1)+1)
eta[eta>1] <- log.mcmc$delta[eta>1] - 1 + 1/log.mcmc$rho[eta>1]
omega <- 1 - log.mcmc$rho*log.mcmc$delta
log.mcmc <- cbind(log.mcmc, eta, omega)


burnin <- 1:80
layout(t(c(1,2)))

for(i in 2:ncol(log.mcmc)){
  plot(log.mcmc[-burnin,i], log.mcmc$iter[-burnin], type = "l", xlab = "Iterations", ylab = parse(text = colnames(log.mcmc)[i])[[1]])
  hist(log.mcmc[-burnin,i], xlab = parse(text = colnames(log.mcmc)[i])[[1]], main = "")
  lines(rep(sim[i],2), c(0, par()$usr[4]), col ="red")
}

