## Test fct° mcmc_ABM
lib <- c("codetools", "ape", "maps", "geiger", "phytools", "bite", "phangorn")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

args <- as.numeric(commandArgs(trailingOnly = T))
prop <- args[1]
nu <- args[2]
omega <- args[3]
sim <- args[4]

cat("#############################\n", prop, nu, omega, sim, "\n#############################\n")

tree <- readRDS(sprintf("Results/Simulations/tree_prop%s_ran.rds", prop))
data <- readRDS(file = sprintf("Results/Simulations/sim_nu%s_omega%s_prop%s_ran.rds", nu, omega, prop))
if(!file.exists(sprintf("Results/adim_h_prop%s_nu%s_omega-%s_sim%s.log", prop, nu, omega, sim))){
  tr <- tree[[sim]]
  traits <- cbind(data[[sim]]$m, exp(data[[sim]]$logv))[1:20,]
  rownames(traits) <- tr$tip.label
  ADIM.obj <- make_ADIM(phy = tr, traits = traits, scale = F, switch = T)
  mcmc_ADIM(ADIM.obj, log.file = sprintf("Results/adim_h_prop%s_nu%s_omega-%s_sim%s.log", prop, nu, omega, sim), sampling.freq = 100, print.freq = 100, ngen = 500000)  
} else {cat("No overwriting of an existing file\n")}
