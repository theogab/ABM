## Test fct° mcmc_ABM
lib <- c("codetools", "ape", "maps", "geiger", "phytools", "bite", "phangorn")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

args <- as.numeric(commandArgs(trailingOnly = T))
mu <- args[1]
nu <- args[2]
omega <- args[3]
sim <- args[4]

cat("#############################\n", mu, nu, omega, sim, "\n#############################\n")

tree <- read.tree(sprintf("Results/Simulations/tree_mu%s.tre", mu))
data <- readRDS(file = sprintf("Results/Simulations/sim_mu%s_nu%s_omega%s.rds", mu, nu, omega))
if(!file.exists(sprintf("Results/adim_ext_mu%s_nu%s_omega-%s_sim%s.log", mu, nu, omega, sim))){
  obs_tips <- which(round(diag(vcv(tree[[sim]])),9)==1)
  tr <- keep.tip(tree[[sim]], obs_tips)
  traits <- cbind(data[[sim]]$m, exp(data[[sim]]$logv))[obs_tips,]
  rownames(traits) <- tr$tip.label
  ADIM.obj <- make_ADIM(phy = tr, traits = traits, scale = F, switch = T)
  mcmc_ADIM(ADIM.obj, log.file = sprintf("Results/adim_ext_mu%s_nu%s_omega-%s_sim%s.log", mu, nu, omega, sim), sampling.freq = 100, print.freq = 100, ngen = 500000)  
} else {cat("No overwriting of an existing file\n")}

