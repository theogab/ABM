## Likelihood under an ABM model
# pars : vector with /[1] sig2_mdp /[2] sig2_rge /[3] nu /[4] omega
lik_ABM <- function(tab, x, pars){
  sig2  <- pars[1:2]
  nu   <- pars[3]
  omega <- pars[4]

  # scenario 1 -> sp1 inherits high midpoint and narrow range
  rg1 <- dnorm(x[tab[,2],2], rnge(omega, nu, x[tab[,1],2],F), sqrt(sig2[2] * tab[,4]), log = T) + dnorm(x[tab[,3],2], rnge(omega, nu, x[tab[,1],2],T), sqrt(sig2[2] * tab[,5]), log = T)
  mp1 <- dnorm(x[tab[,2],1], midpt(omega, nu, x[tab[,1],],F,F), sqrt(sig2[1] * tab[,4]), log = T) + dnorm(x[tab[,3],1], midpt(omega, nu, x[tab[,1],],T,T), sqrt(sig2[1] * tab[,5]), log = T)
  # scenario 2 -> sp1 inherits high midpoint and wide range
  rg2 <- dnorm(x[tab[,2],2], rnge(omega, nu, x[tab[,1],2],T), sqrt(sig2[2] * tab[,4]), log = T) + dnorm(x[tab[,3],2], rnge(omega, nu, x[tab[,1],2],F), sqrt(sig2[2] * tab[,5]), log = T)
  mp2 <- dnorm(x[tab[,2],1], midpt(omega, nu, x[tab[,1],],T,F), sqrt(sig2[1] * tab[,4]), log = T) + dnorm(x[tab[,3],1], midpt(omega, nu, x[tab[,1],],F,T), sqrt(sig2[1] * tab[,5]), log = T)
  # scenario 3 -> sp1 inherits low midpoint and narrow range
  mp3 <- dnorm(x[tab[,2],1], midpt(omega, nu, x[tab[,1],],F,T), sqrt(sig2[1] * tab[,4]), log = T) + dnorm(x[tab[,3],1], midpt(omega, nu, x[tab[,1],],T,F), sqrt(sig2[1] * tab[,5]), log = T)
  # scenario 4 -> sp1 inherits low midpoint and wide range
  mp4 <- dnorm(x[tab[,2],1], midpt(omega, nu, x[tab[,1],],T,T), sqrt(sig2[1] * tab[,4]), log = T) + dnorm(x[tab[,3],1], midpt(omega, nu, x[tab[,1],],F,F), sqrt(sig2[1] * tab[,5]), log = T)
  
  #cat("narrow range:", rnge(omega, nu, x[tab[,1],2],F), "\n")
  #cat("low mdp:", midpt(omega, nu, x[tab[,1],],F,T), "\n")
  #cat("high mdp:", midpt(omega, nu, x[tab[,1],],F,F), "\n")
  #cat("wide range:", rnge(omega, nu, x[tab[,1],2],T), "\n")
  #cat("low mdp:", midpt(omega, nu, x[tab[,1],],T,T), "\n")
  #cat("high mdp:", midpt(omega, nu, x[tab[,1],],T,F), "\n")
  
  minilik <- cbind(mp1+rg1,mp2+rg2,mp3+rg1,mp4+rg2)
  best <- apply(minilik, 1, max)
  scaleminilik <- minilik-best
  lik <- sum(log(apply(exp(scaleminilik), 1, sum)) + best)
  return(lik)
} 

# Range inheritance function
rnge <- function(omega, nu, X0, wide = F){
  if(wide) return(log(exp(X0)/2*(omega*(1-nu^2)/(1+nu)+1+nu)))
  else return(log(exp(X0)/2*(omega*(1-nu^2)/(1+nu)+1-nu)))
}

# Midpoint inheritance function
midpt <- function(omega, nu, x0, wide = F, low = F){
  up <- x0[,1] + exp(x0[,2])/2
  down <- x0[,1] - exp(x0[,2])/2
  rg <- rnge(omega, nu, x0[,2], wide)
  if(low) return(down + exp(rg)/2)
  else return(up - exp(rg)/2)
}

build_table <- function (tree){
  
  tree  <- reorder(tree, "cladewise")
  ntips <- length(tree$tip.label)
  root <- ntips+1
  
  # the table with ancestor, the two descendants and the two branch length
  table.tree<-matrix(NA,tree$Nnode,5)
  colnames(table.tree)<-c("anc","desc1","desc2","b.len1", "b.len2")
  
  #contains the root and all the nodes 
  table.tree[,1] <- unique(tree$edge[,1])
  
  for(i in 1:dim(table.tree)[1]){
    table.tree[i,c("desc1","desc2")] <- tree$edge[which(tree$edge[,1]==table.tree[i,1]),2]
    table.tree[i,c("b.len1","b.len2")] <- tree$edge.length[which(tree$edge[,1]==table.tree[i,1])]
  }
  return(table.tree)
}


## Make a mapping object to parse to mcmc ABM
make_ABM <- function(phy, traits, scale = F){
  
  ### validity test ###
  if (geiger::name.check(phy, traits) != "OK") {
    stop("Species do not match in tree and traits")
  }
  
  
  ### scale height ###
  if(scale){
    t.len <- max(branching.times(phy))
    phy$edge.length <- phy$edge.length/t.len
  }
  
  ABM <- list()
  
  ### Global variables ###
  ABM$data$traits <- cbind(traits[phy$tip.label,1], log(traits[phy$tip.label,2]))
  ABM$data$n      <- length(phy$tip.label)
  ABM$data$tree   <- phy
  ABM$data$nodes  <- (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
  ABM$data$scale  <- scale
  ABM$data$tab <- build_table(phy)
  
  
  ### Likelihood parameters ###
  ABM$lik <- lik_ABM
  ABM$ws <- c(0.2, 1.2, 0.1, 0.1, 0.1, 0.1)
  ABM$nodes <- matrix(0,nrow = ABM$data$tree$Nnode, ncol = 2)
  ABM$nodes[1:ABM$data$tree$Nnode,1] <- sapply(prop.part(ABM$data$tree), function(n){
    mean(ABM$data$traits[n,1])
  })
  ABM$nodes[1:ABM$data$tree$Nnode,2] <- sapply(prop.part(ABM$data$tree), function(n){
    mean(ABM$data$traits[n,2])/(.5 + 1/(2*length(n)))
  })
  ABM$init <- c(runif(2, 0.01, 0.99), 0.5, 0.5)
  ABM$proposals <- list()
  ABM$priors <- list()
  for (i in 1:ABM$data$tree$Nnode) {
    ABM$priors[[i]] <- jive:::hpfun("Normal", c(ABM$nodes[i,1],.5))
    ABM$priors[[i+ABM$data$tree$Nnode]] <- jive:::hpfun("Lognormal", c(ABM$nodes[i,2],.5))
  }
  ABM$proposals[[1]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[2]] <- jive:::proposal("slidingWin")
  ABM$proposals[[3]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[4]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[5]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[6]] <- jive:::proposal("slidingWinAbs")
  ABM$priors[[2*ABM$data$tree$Nnode+1]] <- jive:::hpfun("Gamma", c(1.1, 5))
  ABM$priors[[2*ABM$data$tree$Nnode+2]] <- jive:::hpfun("Gamma", c(1.1, 5))
  ABM$priors[[2*ABM$data$tree$Nnode+3]] <- jive:::hpfun("Uniform", c(0,0.99))
  ABM$priors[[2*ABM$data$tree$Nnode+4]] <- jive:::hpfun("Uniform", c(0,1))
  ABM$update.freq <- c(0.3,0.3,0.1,0.1,0.1,0.1)
  
  #### Prepare headers of log file ####
  
  ABM$header <- c("iter", "posterior", sprintf("range_node_n%s", 1:ABM$data$tree$Nnode), sprintf("midpoint_node_n%s", 1:ABM$data$tree$Nnode),"sigma^2_rge", "sigma^2_mdp", "nu", "omega", "log.lik", "acc")			
  
  return(ABM)
  
}


## mcmc chain
mcmc_ABM <- function(ABM, log.file = "my_ABM_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                     ngen = 5000000, psml = F){
  
  # General syntax
  # 0 : used for claculations
  # 1 : not used for calculations
  
  ## initial conditions
  cat("setting initial conditions\n")
  
  # likelihood level
  x0 <- rbind(ABM$data$traits, ABM$nodes)
  pars0 <- ABM$init
  lik0 <- ABM$lik(ABM$data$tab, x0, pars0)
  prior0 <- unlist(mapply(do.call, ABM$priors, lapply(c(ABM$nodes[,1], ABM$nodes[,2], ABM$init), list))[1,])
  
  # mcmc parameters
  cat("generation\tposterior\n")
  cat(sprintf("%s\t", ABM$header), "\n", append = FALSE, file = log.file)
  update.freq <- cumsum(ABM$update.freq/sum(ABM$update.freq))
  proposals <- c(0,0,0,0,0,0) # 1st number: update ancestral states; 2nd: update sig2, 3rd: update rho, 4th: update delta
  proposals.accepted <- c(0,0,0,0,0,0)
  
  # posterior level
  post0 <- lik0 + sum(prior0)
  
  ## Start iterations
  for (i in 1:ngen) {
    
    # test whether we update ancestral states (r[1]), sig2 (r[2]), nu (r[3]) or omega (r[4])
    r <- runif(1) <= update.freq
    r[2] <- r[2] & !r[1]
    r[3] <- r[3] & !r[2] & !r[1]
    r[4] <- r[4] & !r[3] & !r[2] & !r[1]
    r[5] <- r[5] & !r[4] & !r[3] & !r[2] & !r[1]
    r[6] <- r[6] & !r[5] & !r[4] & !r[3] & !r[2] & !r[1]
    
    proposals[r] <- proposals[r] + 1
    lik1 <- lik0
    prior1 <- prior0  
    
    if (r[1]) # update ancestral midpoints
    { 
      ind <- sample(ABM$data$nodes, 2, replace = F) # 5 is just for now number of ancestral states updated
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      x1 <- x0
      tmp <- ABM$proposals[[1]](i = x0[ind,1], d = ABM$ws[1], u)
      x0[ind,1] <- tmp$v
      hasting.ratio <- sum(tmp$lnHastingsRatio)
    } else {
      pars1 <- pars0 # ancient parameter values are kept
    }
    
    if (r[2]) # update ancestral ranges
    { 
      ind <- sample(ABM$data$nodes, 2, replace = F) # 5 is just for now number of ancestral states updated
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      x1 <- x0
      tmp <- ABM$proposals[[2]](i = x0[ind,2], d = ABM$ws[2], u)
      x0[ind,2] <- tmp$v
      hasting.ratio <- sum(tmp$lnHastingsRatio)
    } else {
      pars1 <- pars0 # ancient parameter values are kept
    }
    
    if (r[3]) # update sig2_mdp
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[3]](i = pars0[1], d = ABM$ws[3], u) #update with proposal function
      pars0[1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[4]) # update sig2_rge
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[4]](i = pars0[2], d = ABM$ws[4], u) #update with proposal function
      pars0[2] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[5]) # update nu
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[5]](i = pars0[3], d = ABM$ws[5], u) #update with proposal function
      pars0[3] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if (r[6]) # update omega
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[6]](i = pars0[4], d = ABM$ws[6], u) #update with proposal function
      pars0[4] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    # Posterior calculation
    prior0 <- unlist(mapply(do.call, ABM$priors, lapply(c(x0[ABM$data$nodes,1], x0[ABM$data$nodes,2], pars0), list))[1,])
    lik0 <- ABM$lik(ABM$data$tab, x0, pars0)
    post1 <- post0
    post0 <- lik0 + sum(prior0)
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0, prior0)))){
      pr <- -Inf
    } else {
     pr <- post0 - post1 + hasting.ratio
    }
    
    if(psml) f <- 0
    else f <- log(runif(1))
    if (pr >= f)#log(runif(1))) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
    } else # cancel changes
    {
      post0 <- post1
      prior0 <- prior1
      lik0 <- lik1
      
      if (r[1] | r[2]){
        x0 <- x1
      } else {
        pars0 <- pars1
      }
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0) {
      cat(sprintf("%s\t", c(i, post0, x0[ABM$data$nodes,1], x0[ABM$data$nodes,2], pars0, lik0,(sum(proposals.accepted)/i))),
          "\n", append=TRUE, file=log.file) 
    }
    
    # Print to screen
    if (i %% print.freq == 0) {
      cat(i,'\t',post0,'\n') 
    }
    
  } # end of for 
  
  # Calculate acceptance rate
  names(proposals) <- c("Ancestral midpoint","Ancestral range", "sigma^2_mdp","sigma^2_rge","nu","omega")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("Ancestral midpoint","Ancestral range", "sigma^2_mdp","sigma^2_rge","nu","omega")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
}
