## Simulate data
sim_ABM <- function(tree, pars_logX = c(.1, 2, -Inf, Inf), pars_m = c(.1, 10, -Inf, Inf), nu = 0.5, omega = 0.5, nsim = 1){
  
  b.len <- tree$edge.length
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  n <- length(tree$tip.label)
  
  # simulate data
  dat <- lapply(1:nsim, function(i){
    
    # standard BM variation
    var_logX <- rnorm(length(b.len), 0, sqrt(pars_logX[1]*b.len))
    var_m <- rnorm(length(b.len), 0, sqrt(pars_m[1]*b.len))
    
    # inheritance assymetry
    Sm <- rep(0, n+tree$Nnode)
    Sx <- rep(0, n+tree$Nnode)
    for(j in 1:(n+tree$Nnode)){
      if(j %in% e2 & Sm[j] == 0 | Sx[j] == 0){
        anc <- e1[e2 == j]
        sis <- e2[e1 == anc]
        Sm[sis] <- sample(c(-1,1), 2, replace = F)
        Sx[sis] <- sample(c(-1,1), 2, replace = F)  
      }
    }
    
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
    m <- c(rep(0, n), pars_m[2], rep(0, tree$Nnode-1))
    logX <- c(rep(0, n), pars_logX[2], rep(0, tree$Nnode-1))
    order <- c(n + (2:tree$Nnode), 1:n)
    for(desc in order){
      anc <- e1[e2 == desc]
      bn <- which(e2 == desc)
      m[desc] <- reflect(midpt(omega, nu, logX[anc], m[anc], Sx[desc], Sm[desc]) + var_m[bn], pars_m[3:4])
      logX[desc] <- reflect(rnge(omega, nu, logX[anc], Sx[desc]) + var_logX[bn], pars_logX[3:4])
    }
    return(list(m = m, logX = logX, Sm = Sm, Sx = Sx))
  })
  
  return(dat)
}