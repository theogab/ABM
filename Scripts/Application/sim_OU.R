### Simulate OU models for variance.

## Simulations ABM
lib <- c("ape", "phytools", "bite")
sapply(lib, library, character.only = T
       #, lib.loc = "/scratch/axiom/FAC/FBM/DBC/nsalamin/default/tgaboria/R"
)

for(n in c(20, 50, 100)){
  for(sig_m in c(0.1, 0.5, 1)){
    for(sig_logX in c(0.1, 0.5, 1)){
          cat("#############################\n", n, sig_m, sig_logX, nu, omega, "\n#############################\n")
          set.seed(300)
          trees <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
          unif.f <- function(n, pars){
            runif(n, pars[1] - exp(pars[2])/2, pars[1] + exp(pars[2])/2)
          }
          Q <- cbind(c(-.2, .2), c(.2, -.2))
          sim <- lapply(trees, function(tree){
            out <- sim_jive(tree, map = sim.history(tree, Q)$maps, models = list(mid = c("BM"), logrange = c("OU", "theta")), pars = list(mid = c(root = 10, sigma_sq = sig_m), logrange = c(root = 0, theta1 = 0, theta2 = log(0.5), sigma_sq = sig_logX, alpha = 0.1)))$pars
            list(m = out[,1], logX = out[,2])
          } )
          saveRDS(sim, file = sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_OU.rds", n,sig_m, sig_logX))
    }
  }
}
 