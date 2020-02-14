## Test fct° mcmc_ABM
lib <- c("ape", "phytools", "bite")
sapply(lib, library, character.only = T
       #, lib.loc = "/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/R"
       )
source("Scripts/mcmc_ABM.R")
source("Scripts/sim_ABM.R")

#args <- as.numeric(commandArgs(trailingOnly = T))

n <- 40
sm <- 0.001
sx <- 0.001
nu <- 0.9
om <- 0
sim_nb <- 1

tree <- pbtree(n = n, scale = 1)
sim <- sim_ABM(tree, pars_logX = c(sx, 0.7, -Inf, Inf), pars_m = c(sm, 10, -Inf, Inf), nu = nu, omega = om)[[1]]
traits <- cbind(sim$m, exp(sim$logX))[1:n,]
rownames(traits) <- tree$tip.label
ABM.obj <- make_ABM(phy = tree, traits = traits, scale = F)
ABM.obj$init$Sm <- sim$Sm
ABM.obj$init$Sx <- sim$Sx
#ABM.obj$init$pars[1:2] <- c(0.001, 0.001)
#ABM.obj$init$m <- sim$m[-(1:n)]
#ABM.obj$init$X <- sim$X[-(1:n)]
ABM.obj$update.freq <- c(0.2, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
lik_ABM(tree = ABM.obj$data$tree, logX = sim$logX, m = sim$m, pars = c(sm, sx, nu, om), Sx = ABM.obj$init$Sx, Sm = ABM.obj$init$Sm)

mcmc_ABM(ABM.obj, log.file = sprintf("Results/abm_gibbs_n-%s_sm-%s_sx-%s_nu-%s_om-%s_sim-%s.log", n, sm, sx, nu, om, sim_nb), sampling.freq = 1000, print.freq = 1000, ngen = 5000000)


m0 <- c(ABM.obj$data$m, ABM.obj$init$m)
logX0 <- c(ABM.obj$data$logX, ABM.obj$init$logX)
Sm0 <- ABM.obj$init$Sm
Sx0 <- ABM.obj$init$Sx
pars0 <- ABM.obj$init$pars
lik0 <- ABM.obj$lik(ABM.obj$data$tree, logX0, m0, pars0, Sx0, Sm0)
prior.m0 <- unlist(mapply(do.call, ABM.obj$priors$m, lapply(c(ABM.obj$init$m), list))[1,])
prior.logX0 <- unlist(mapply(do.call, ABM.obj$priors$logX, lapply(c(ABM.obj$init$logX), list))[1,])
prior.pars0 <- unlist(mapply(do.call, ABM.obj$priors$pars, lapply(c(ABM.obj$init$pars), list))[1,])
cat(lik0, "\n")

for(sX in c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000)){
  if(sX == 0.0001){
    cat("#### Wrong ancestral states ####\n")
    cat("sX\tlik\n")
  }
  pars1 <- pars0
  pars1[2] <- sX
  st <- gibbs_move(tree = ABM.obj$data$tree, n = ABM.obj$data$n, m = m0, logX = logX0, Sm = Sm0, SX = Sx0, pars = pars0, priors = ABM.obj$priors)
  logX1 <- st$logX0
  m1 <- st$m0
  lik1 <- ABM.obj$lik(ABM.obj$data$tree, logX1, m1, pars1, Sx0, Sm0)
  cat(sX, ":", lik1, "\n")
}

for(sX in c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000)){
  if(sX == 0.0001){
    cat("#### Good ancestral states ####\n")
    cat("sX\tlik\n")
  }
  pars1 <- pars0
  pars1[2] <- sX
  lik1 <- ABM.obj$lik(ABM.obj$data$tree, sim$logX, sim$m, pars1, Sx0, Sm0)
  cat(sX, ":", lik1, "\n")
}



for(sm in c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000)){
  if(sm == 0.0001){
    cat("#### Wrong ancestral states ####\n")
    cat("sm\tlik\n")
  }
  pars1 <- pars0
  pars1[1] <- sm
  lik1 <- ABM.obj$lik(ABM.obj$data$tree, logX0, m0, pars1, Sx0, Sm0)
  cat(sm, ":", lik1, "\n")
}

for(sm in c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000)){
  if(sm == 0.0001){
    cat("#### Good ancestral states ####\n")
    cat("sm\tlik\n")
  }
  pars1 <- pars0
  pars1[1] <- sm
  lik1 <- ABM.obj$lik(ABM.obj$data$tree, sim$logX, sim$m, pars1, Sx0, Sm0)
  cat(sm, ":", lik1, "\n")
}

