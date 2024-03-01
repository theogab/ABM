## Test fct° mcmc_ABM
lib <- c("ape", "phytools", "bite", "phangorn")
sapply(lib, library, character.only = T, lib.loc = "/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/R")
source("Scripts/mcmc_ABM.R")
source("Scripts/sim_ABM.R")
source("Scripts/plot_states.R")

args <- as.numeric(commandArgs(trailingOnly = T))
n <- args[1]
sigm <- args[2]
siglogX <- args[3]
nu <- args[4]
omega <- args[5]
sim <- args[6]

cat("#############################\n", sigm, siglogX, nu, omega, sim, "\n#############################\n")

tree <- read.tree(sprintf("Results/Simulations/tree_n%s_sigm%s_siglogX%s_nu%s_omega%s.tre", n, sigm, siglogX, nu, omega))
data <- readRDS(file = sprintf("Results/Simulations/sim_n%s_sigm%s_siglogX%s_nu%s_omega%s.rds", n,sigm, siglogX, nu, omega))
if(!file.exists(sprintf("Results/abm_gibbs_n%s_sigm%s_slogX%s_nu%s_omega-%s_sim%s.log", n, sigm, siglogX, nu, omega, sim))){
  traits <- cbind(data[[sim]]$m, exp(data[[sim]]$logX))[1:n,]
  rownames(traits) <- tree[[sim]]$tip.label
  ABM.obj <- make_ABM(phy = tree[[sim]], traits = traits, scale = F, switch = T)
  mcmc_ABM(ABM.obj, log.file = sprintf("Results/abm_gibbs_n%s_sigm%s_slogX%s_nu%s_omega-%s_sim%s.log", n, sigm, siglogX, nu, omega, sim), sampling.freq = 100, print.freq = 100, ngen = 20000, psml = F)  
}

