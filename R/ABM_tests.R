### Jive-relaxed simus 
args <- as.numeric(commandArgs(trailingOnly = T))
cat(args, "\n")
rho <- args[1]
h <- args[2]
sig2 <- args[3]
nsim <- args[4]

lib <- c("ape","phytools")
sapply(lib, library, character.only = T)

# Simple assymetry
# Assymetric Brownian Motion (ABM) Gamma /[1,2] controls overlap and rho /[0,0.5] controls assymetry 
set.seed(300)
tree1 <- pbtree(n = 10, scale = 100)
tree2 <- pbtree(n = 50, scale = 100)
tree3 <- pbtree(n = 200, scale = 100)

## Simulate data
sim_ABM <- function(tree, sig2 = .1, theta0 = 2, rho = 0.25, h = .5, bounds = c(0, Inf), nsim = 1){
  tree <- reorder(tree, "cladewise")
  b.len <- tree$edge.length
  rho.val <- c(rho, 1 - rho)
  
  # simulate data
  X <- sapply(1:nsim, function(i){
    
    # standard BM variation
    var <- rnorm(length(b.len), 0, sqrt(sig2*b.len))
    
    # inheritance assymetry
    as <- rep(0, nrow(tree$edge))
    n <- 1
    while (any(as == 0)){
      sis <- which(tree$edge[,1] == tree$edge[n,1])
      as[sis] <- sample(1:2, 2, replace = F)
      n <- which(as == 0)[1]
    }
    
    ## function to respect boubds
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
    val <- c(rep(0, tree$Nnode+1), theta0, rep(0, tree$Nnode-1))
    for(k in 1:nrow(tree$edge)){
      anc <- tree$edge[k, 1]
      desc <- tree$edge[k, 2]
      val[desc] <- h * val[anc] + (1 - h) * rho.val[as[k]] * val[anc] + var[k]
      val[desc] <- reflect(val[desc], bounds)
    }
    val
  })
  
  return(X)
  
}

## Likelihood under an ABM model
# pars : vector with [1] sig2 /[2] theta0 /[3] rho /[4] h
lik_ABM <- function(tree, x, pars){
  tree  <- reorder(tree, "cladewise")
  tips  <- x[1:length(tree$tip.label)]
  root  <- pars[2]
  nodes <- x[tree$edge[tree$edge[,2] %in% 1:length(tree$tip.label), 1]-1]
  sig2  <- pars[1]
  rho   <- pars[3]
  h     <- pars[4]
  b.len <- tree$edge.length[tree$edge[,2] %in% 1:length(tree$tip.label)]
  
  left.prob <- dnorm(tips, h *  nodes + (1 - h) * rho * nodes, sqrt(sig2 * b.len))
  right.prob <- dnorm(tips, h * nodes + (1 - h) * (1 - rho) * nodes, sqrt(sig2 * b.len))
  
  lik <- sum(log(1/2*(left.prob + right.prob)))
  
  return(lik)
} 

## simulations
set.seed(300)
sim1 <- sim_ABM(tree1, rho = rho, h = h, sig2 = sig2, nsim = nsim)

set.seed(300)
sim2 <- sim_ABM(tree2, rho = rho, h = h, sig2 = sig2, nsim = nsim)

set.seed(300)
sim3 <- sim_ABM(tree3, rho = rho, h = h, sig2 = sig2, nsim = nsim)


## likelihoods
sig.parse <- c(0.01,.1,5)
rho.parse <- seq(0.1,0.5, by = 0.1)
comb <- cbind(sig = rep(sig.parse, each = length(rho.parse)), rho = rep(rho.parse, length(sig.parse)))
lik1 <- apply(sim1, 2, function(x){
  apply(comb, 1, function(y){
    lik_ABM(tree1, x[-(length(tree1$tip.label)+1)], c(y[1], 2, y[2], 0.5))
  })
})
lik1 <- cbind(comb,lik1)

lik2 <- apply(sim2, 2, function(x){
  apply(comb, 1, function(y){
    lik_ABM(tree2, x[-(length(tree2$tip.label)+1)], c(y[1], 2, y[2], 0.5))
  })
})
lik2 <- cbind(comb,lik2)

lik3 <- apply(sim3, 2, function(x){
  apply(comb, 1, function(y){
    lik_ABM(tree3, x[-(length(tree3$tip.label)+1)], c(y[1], 2, y[2], 0.5))
  })
})
lik3 <- cbind(comb,lik3)

lik <- list(tree1 = list(sim1, lik1), tree2 = list(sim2, lik2), tree3 = list(sim3, lik3))

save(lik, file = sprintf("Results/lik_rho-%s_h-%s_sig-%s.Rdata", rho, h, sig2))