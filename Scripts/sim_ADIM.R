## Simulate data with normal ADIM
sim_ADIM <- function(tree, pars_logv = c(.01, 2, -Inf, Inf), pars_m = c(.1, 10, -Inf, Inf), nu = 0.5, omega = 0.5, nsim = 1, nodemap = F){
  
  # simulate data
  dat <- lapply(1:nsim, function(i){
    
    if(class(tree) == "multiPhylo"){
      if(length(tree)!=nsim) error("tree must be of the same length as nsim")
      phy <- tree[[i]]
    } else phy <- tree
    
    b.len <- phy$edge.length
    e1 <- phy$edge[,1]
    e2 <- phy$edge[,2]
    n <- length(phy$tip.label)
    
    # standard BM variation
    var_logv <- rnorm(length(b.len), 0, sqrt(pars_logv[1]*b.len))
    var_m <- rnorm(length(b.len), 0, sqrt(pars_m[1]*b.len))
    
    # inheritance assymetry
    Sn <- sample(c(-1,1), phy$Nnode, replace = T)
    So <- sample(c(-1,1), phy$Nnode, replace = T)
  
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
    m <- c(rep(0, n), pars_m[2], rep(0, phy$Nnode-1))
    logv <- c(rep(0, n), pars_logv[2], rep(0, phy$Nnode-1))
    order <- n + (1:phy$Nnode)
    for(anc in order){
      desc <- e2[e1 == anc]
      bn <- which(e1 == anc)
      if(nodemap){
        nup <- ifelse(phy$node.label[anc - n] == "1", 0, nu)
        omegap <- ifelse(phy$node.label[anc - n] == "1", 0, omega)
      } else {
        nup <- nu
        omegap <- omega
      }
      m[desc] <- sapply(mud(nup, omegap, exp(logv[anc]), m[anc], Sn[anc - n], So[anc - n], c(-1,1)) + var_m[bn],reflect, pars_m[3:4])
      logv[desc] <- sapply(log(sigd(nup, omegap, exp(logv[anc]), Sn[anc - n], So[anc - n], c(-1,1))) + var_logv[bn],reflect, pars_logv[3:4])
    }
    return(list(m = m, logv = logv, Sn = Sn, So = So, n.map = phy$node.label))
  })
  
  return(dat)
}
