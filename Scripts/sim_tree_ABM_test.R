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

## Simulate data
sim_ABM <- function(tree, pars_rge = c(.1, 2, -Inf, Inf), pars_mdp = c(.1, 10, -Inf, Inf), nu = 0.5, omega = 0.5, nsim = 1){
  
  tree <- reorder(tree, "cladewise")
  b.len <- tree$edge.length
  
  # simulate data
  X <- lapply(1:nsim, function(i){
    
    # standard BM variation
    var_rge <- rnorm(length(b.len), 0, sqrt(pars_rge[1]*b.len))
    var_mdp <- rnorm(length(b.len), 0, sqrt(pars_mdp[1]*b.len))
    
    # inheritance assymetry
    as <- rep(0, nrow(tree$edge))
    n <- 1
    while (any(as == 0)){
      sis <- which(tree$edge[,1] == tree$edge[n,1])
      as[sis] <- sample(1:2, 2, replace = F)
      n <- which(as == 0)[1]
    }
    as_rge <- as.logical(as-1)
    
    as <- rep(0, nrow(tree$edge))
    n <- 1
    while (any(as == 0)){
      sis <- which(tree$edge[,1] == tree$edge[n,1])
      as[sis] <- sample(1:2, 2, replace = F)
      n <- which(as == 0)[1]
    }
    as_mdp <- as.logical(as-1)
    
    ## function to respect bounds
    reflect <- function(yy, bounds) {
      while (yy < bounds[1] || yy > bounds[2]) {
        if (yy < bounds[1]) 
          yy <- 2 * bounds[1] - yy
        if (yy > bounds[2]) 
          yy <- 2 * bounds[2] - yy
      }
      return(yy)
    }
    
    # go from roots to tips
    val <- cbind(c(rep(0, length(tree$tip.label)), pars_mdp[2], rep(0, tree$Nnode-1)), c(rep(0, length(tree$tip.label)), pars_rge[2], rep(0, tree$Nnode-1)))
    for(k in 1:nrow(tree$edge)){
      anc <- tree$edge[k, 1]
      desc <- tree$edge[k, 2]
      val[desc,1] <- reflect(midpt(omega, nu, t(as.matrix(val[anc,])), as_rge[k], as_mdp[k]) + var_mdp[k], pars_mdp[3:4])
      val[desc,2] <- reflect(rnge(omega, nu, val[anc,2], as_rge[k]) + var_rge[k], pars_rge[3:4])
    }
    return(val)
  })
  
  return(X)
}

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

