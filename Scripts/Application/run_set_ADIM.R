## Test fct° mcmc_ABM
lib <- c("codetools", "ape", "maps", "geiger", "phytools", "bite", "phangorn")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

args <- as.numeric(commandArgs(trailingOnly = T))
n <- args[1]
sigm <- args[2]
siglogv <- args[2]
nu <- args[4]
omega <- args[5]
sim <- args[3]

cat("#############################\n", sigm, siglogv, nu, omega, sim, "\n#############################\n")

tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
data <- readRDS(file = sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_nu%s_omega%s.rds", n,sigm, siglogv, nu, omega))
if(!file.exists(sprintf("Results/adim_n%s_sigm%s_siglogv%s_nu%s_omega-%s_sim%s.log", n, sigm, siglogv, nu, omega, sim))){
  traits <- cbind(data[[sim]]$m, exp(data[[sim]]$logv))[1:n,]
  rownames(traits) <- tree[[sim]]$tip.label
  ADIM.obj <- make_ADIM(phy = tree[[sim]], traits = traits, scale = F, switch = T)
  mcmc_ADIM(ADIM.obj, log.file = sprintf("Results/adim_n%s_sigm%s_siglogv%s_nu%s_omega-%s_sim%s.log", n, sigm, siglogv, nu, omega, sim), sampling.freq = 100, print.freq = 100, ngen = 500000)  
} else {cat("No overwriting of an existing file\n")}

