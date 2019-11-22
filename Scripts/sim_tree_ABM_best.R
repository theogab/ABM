### Jive-relaxed simus 
args <- as.numeric(commandArgs(trailingOnly = T))
cat(args, "\n")
sig2_rge <- args[1]
sig2_mdp <- args[2]
nu <- args[3]
omega <- args[4]
nsim <- args[5]


source('Scripts/mcmc_ABM_best.R')
lib <- c("ape","phytools","jive","geiger")
sapply(lib, library, character.only = T, lib.loc = "/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/R")

# Simple assymetry
# Assymetric Brownian Motion (ABM)
tree1 <- pbtree(n = 20, scale = 100)
tree2 <- pbtree(n = 100, scale = 100)
tree3 <- pbtree(n = 500, scale = 100)

## simulations
sim1 <- sim_ABM(tree1, nu = nu, omega = omega, pars_rge = c(sig2_rge, 2, 0, Inf), pars_mdp = c(sig2_mdp, 10, -Inf, Inf), nsim = 1)

sim2 <- sim_ABM(tree2, nu = nu, omega = omega, pars_rge = c(sig2_rge, 2, 0, Inf), pars_mdp = c(sig2_mdp, 10, -Inf, Inf), nsim = 1)

sim3 <- sim_ABM(tree3, nu = nu, omega = omega, pars_rge = c(sig2_rge, 2, 0, Inf), pars_mdp = c(sig2_mdp, 10, -Inf, Inf), nsim = 1)

sim_tree <- list(sim1, sim2, sim3)

save(sim_tree, file = sprintf("Results/sim_tree_best_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_sim-%s.rda", nu, omega, sig2_rge, sig2_mdp, nsim))

## mcmc chains

for(tree in 1:3){
    eval(parse(text = sprintf("phy <- tree%s", tree))) 
    traits <- sim_tree[[tree]][[1]][1:length(phy$tip.label),]
    rownames(traits) <- phy$tip.label
    my_ABM <- make_ABM(phy, traits, scale = F)
    mcmc_ABM(my_ABM, log.file = sprintf("Results/ABM_mcmc_best_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, omega, sig2_rge, sig2_mdp, tree, nsim),
             sampling.freq = 500, print.freq = 500, ngen = 200000) 
}

