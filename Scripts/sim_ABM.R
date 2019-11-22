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
    return(cbind(val[,1], exp(val[,2])))
  })
  
  return(X)
}