## Simulations ABM
lib <- c("ape", "phytools", "bite","matrixStats", "TreeSim")
sapply(lib, library, character.only = T
       #, lib.loc = "/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/R"
       )
source("Scripts/mcmc_ADIM.R")
source("Scripts/sim_ADIM.R")
source("Scripts/plot_states.R")

set.seed(32)
for(n in c(20,50,100)){
  tree <- pbtree(n = n, scale = 1, nsim = 100)
  write.tree(tree, file = sprintf("Results/Simulations/tree_n%s.tre", n))
}

### With extinction
set.seed(493)
for(mu in c(0.2,0.5,0.8)){
  tree <- sim.bd.taxa(20,100,1,mu)
  tree <- lapply(tree, function(tr){
    tr$edge.length <- tr$edge.length/max(vcv(tr))
    class(tr) <- "phylo"
    return(tr)
  })
  class(tree)<- "multiPhylo"
  write.tree(tree, file = sprintf("Results/Simulations/tree_mu%s.tre", mu))
}

### With ITS
set.seed(543)
for(n in c(25,30,40)){
  tree <- sim.bd.taxa(n,100,1,0)
  tree <- lapply(tree, function(tr){
    tr$edge.length <- tr$edge.length/max(vcv(tr))
    class(tr) <- "phylo"
    return(tr)
  })
  class(tree)<- "multiPhylo"
  write.tree(tree, file = sprintf("Results/Simulations/tree_n%s.tre", n))
}


### Full birth simulations
seeds <- c(45,98,984,49479,49,974,23189,84984,6849)
for(n in c(20,50,100)){
  tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
  s <- 1
	for(sig_m in c(0.1, 0.5, 1)){
 	 	for(sig_logv in c(0.1, 0.5, 1)){
 	 	  seed <- seeds[s]
 	 	  s <- s + 1
    	for(nu in c(0, 0.2, 0.5)){
      		for(omega in c(0, 0.2, 0.5)){
        			cat("#############################\n", n, sig_m, sig_logv, nu, omega, "\n#############################\n")
        			set.seed(seed)
        			sim <- sim_ADIM(tree, pars_logv = c(sig_logv, 0, -Inf, Inf), pars_m = c(sig_m, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 100)
        			saveRDS(sim, file = sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_nu%s_omega%s.rds", n,sig_m, sig_logv, nu, omega))
      		}
    		}
  		}
	}
}

## With extinction
set.seed(8745)
for(mu in c(0.2,0.5, 0.8)){
  tree <- read.tree(sprintf("Results/Simulations/tree_mu%s.tre", mu))
  for(nu in c(0, 0.2, 0.5)){
        for(omega in c(0, 0.2, 0.5)){
          cat("#############################\n", mu, nu, omega, "\n#############################\n")
          sim <- sim_ADIM(tree, pars_logv = c(0.1, 0, -Inf, Inf), pars_m = c(0.1, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 100)
          saveRDS(sim, file = sprintf("Results/Simulations/sim_mu%s_nu%s_omega%s.rds", mu, nu, omega))
    }
  }
}


#### DEBUGING
tree <- read.tree(sprintf("Results/Simulations/tree_mu%s.tre", 0))
for(nu in c(0.5)){
  for(omega in c(0)){
    cat("#############################\n", mu, nu, omega, "\n#############################\n")
    sim <- sim_ADIM(tree, pars_logv = c(0.1, 0, -Inf, Inf), pars_m = c(0.1, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 10)
    saveRDS(sim, file = sprintf("Results/Simulations/sim_mu%s_nu%s_omega%s.rds", 0, nu, omega))
  }
}

tree <- read.tree(sprintf("Results/Simulations/test_tree_n%s.tre", 20))
for(nu in c(0.5)){
  for(omega in c(0)){
    cat("#############################\n", mu, nu, omega, "\n#############################\n")
    sim <- sim_ADIM(tree, pars_logv = c(0.1, 0, -Inf, Inf), pars_m = c(0.1, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 10)
    saveRDS(sim, file = sprintf("Results/Simulations/test_sim_nu%s_omega%s.rds", nu, omega))
  }
}
#####################

## With ITS
set.seed(1028)
for(n in c(25,30,40)){
  tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
  for(nu in c(0, 0.2, 0.5)){
    for(omega in c(0, 0.2, 0.5)){
      cat("#############################\n", n, nu, omega, "\n#############################\n")
      sim <- sim_ADIM(tree, pars_logv = c(0.1, 0, -Inf, Inf), pars_m = c(0.1, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 100)
      saveRDS(sim, file = sprintf("Results/Simulations/sim_n%s_nu%s_omega%s.rds", n, nu, omega))
    }
  }
}

## heterogeneity in the process
tree <- read.tree("Results/Simulations/tree_n20.tre")
for(prop in c(0.2,0.5,0.8)){
  set.seed(300)
  tree_ran <- lapply(tree, function(x){
    x$node.label <- as.character(rbinom(x$Nnode,1,prop) + 1)
    return(x)
  })
  class(tree_ran) <- "multiPhylo"
  saveRDS(tree_ran, file = sprintf("Results/Simulations/tree_prop%s_ran.rds", prop))
  sig_m <- sig_logv <- 0.1
  for(nu in c(0, 0.2, 0.5)){
    for(omega in c(0, 0.2, 0.5)){
      cat("#############################\n", prop, nu, omega, "\n#############################\n")
      sim <- sim_ADIM(tree_ran, pars_logv = c(sig_logv, 0, -Inf, Inf), pars_m = c(sig_m, 10, -Inf, Inf), nu = nu, omega = omega, nsim = 100, nodemap = T)
      saveRDS(sim, file = sprintf("Results/Simulations/sim_nu%s_omega%s_prop%s_ran.rds", nu, omega, prop))
    }
  }
}

cat(file = "Scripts/par.set", append = F)
for(n in c(20,50,100)){
  for(sig in c(0.1, 0.5, 1)){
    for(sim in 1:100){
      cat(n, sig, sim, "\n", file = "Scripts/par.set", append = T)
    }
  }
}

cat(file = "Scripts/par.set.ext", append = F)
for(mu in c(0.2,0.5,0.8)){
  nu <- 0
  omega <- 0.5
  for(sim in 1:100){
    cat(mu, nu, omega, sim, "\n", file = "Scripts/par.set.ext", append = T)
  }
  nu <- 0.2
  omega <- 0.2
  for(sim in 1:100){
    cat(mu, nu, omega, sim, "\n", file = "Scripts/par.set.ext", append = T)
  }
  nu <- 0.5
  omega <- 0.5
  for(sim in 1:100){
    cat(mu, nu, omega, sim, "\n", file = "Scripts/par.set.ext", append = T)
  }
  nu <- 0.5
  omega <- 0
  for(sim in 1:100){
    cat(mu, nu, omega, sim, "\n", file = "Scripts/par.set.ext", append = T)
  }
}

cat(file = "Scripts/par.set.its", append = F)
for(n in c(25,30,40)){
  nu <- 0
  omega <- 0.5
  for(sim in 1:100){
    cat(n, nu, omega, sim, "\n", file = "Scripts/par.set.its", append = T)
  }
  nu <- 0.2
  omega <- 0.2
  for(sim in 1:100){
    cat(n, nu, omega, sim, "\n", file = "Scripts/par.set.its", append = T)
  }
  nu <- 0.5
  omega <- 0.5
  for(sim in 1:100){
    cat(n, nu, omega, sim, "\n", file = "Scripts/par.set.its", append = T)
  }
  nu <- 0.5
  omega <- 0
  for(sim in 1:100){
    cat(n, nu, omega, sim, "\n", file = "Scripts/par.set.its", append = T)
  }
}  

cat(file = "Scripts/par.set.h", append = F)
for(prop in c(0.2,0.5,0.8)){
  nu <- 0
  omega <- 0.5
  for(sim in 1:100){
    cat(prop, nu, omega, sim, "\n", file = "Scripts/par.set.h", append = T)
  }
  nu <- 0.2
  omega <- 0.2
  for(sim in 1:100){
    cat(prop, nu, omega, sim, "\n", file = "Scripts/par.set.h", append = T)
  }
  nu <- 0.5
  omega <- 0.5
  for(sim in 1:100){
    cat(prop, nu, omega, sim, "\n", file = "Scripts/par.set.h", append = T)
  }
  nu <- 0.5
  omega <- 0
  for(sim in 1:100){
    cat(prop, nu, omega, sim, "\n", file = "Scripts/par.set.h", append = T)
  }
}  

