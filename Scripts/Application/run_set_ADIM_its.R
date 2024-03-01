## Test fct° mcmc_ABM
lib <- c("codetools", "ape", "maps", "geiger", "phytools", "bite", "phangorn")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

args <- as.numeric(commandArgs(trailingOnly = T))
n <- args[1]
nu <- args[2]
omega <- args[3]
sim <- args[4]

cat("#############################\n", n, nu, omega, sim, "\n#############################\n")

tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
data <- readRDS(file = sprintf("Results/Simulations/sim_n%s_nu%s_omega%s.rds", n, nu, omega))
if(!file.exists(sprintf("Results/adim_its_n%s_nu%s_omega-%s_sim%s.log", n, nu, omega, sim))){
  obs_tips <- sample(1:n, 20)
  ## save kept tips
  cat(obs_tips, file = sprintf("Results/kept_tips_its_n%s_nu%s_omega-%s_sim%s.log", n, nu, omega, sim))
  tr <- keep.tip(tree[[sim]], obs_tips)
  traits <- cbind(data[[sim]]$m, exp(data[[sim]]$logv))[tr$tip.label,]
  ADIM.obj <- make_ADIM(phy = tr, traits = traits, scale = F, switch = T)
  mcmc_ADIM(ADIM.obj, log.file = sprintf("Results/adim_its_n%s_nu%s_omega-%s_sim%s.log", n, nu, omega, sim), sampling.freq = 100, print.freq = 100, ngen = 500000)  
} else {cat("No overwriting of an existing file\n")}

