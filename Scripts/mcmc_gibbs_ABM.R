require(hypervolume)
require(ape)

## Likelihood under an ABM model
# pars : vector with /[1] sig2nu /[2] sig2omega /[3] sig2theta /[4] nu0 /[5] omega0 /[6:n] thetas (ancestral states)
# x : table of observations with /[,1] volume of the union of descendants' hypervolumes /[,2] asymmetry between descendants' hypervolumes /[,3] overalp between descendants' hypervolumes
# tab : 
lik_ABM <- function(tab, x, pars){
  sig.nu <- pars[1]
  sig.omega <- pars[2]
  sig.theta  <- pars[3]
  nu0   <- pars[4]
  omega0 <- pars[5]
  logthetas <- pars[6:length(pars)]
  
  theta.prob <- dnorm(log(x[tab[,3],1]), logthetas[tab[,3]], sqrt(sig.theta * tab[,4]))
  nu.prob <- dnorm(x[tab[,3],2], exp(logthetas[tab[,3]])*nu0, sqrt(sig.nu * tab[,4]))
  omega.prob <- dnorm(x[tab[,3],3], exp(logthetas[tab[,3]])*omega0, sqrt(sig.omega * tab[,4]))
  
  lik <- sum(log(theta.prob) + log(nu.prob) + log(omega.prob))
  
  return(lik)
} 


build_table <- function (tree){
  
  tree  <- reorder(tree, "cladewise")
  ntips <- length(tree$tip.label)
  root <- ntips+1
  
  # the table with ancestor, the two descendants and the two branch length
  table.tree<-matrix(NA,tree$Nnode,4)
  colnames(table.tree)<-c("anc", "brlen", "node", "nodages")
  
  #contains the root and all the nodes 
  table.tree[1,1] <- NA
  table.tree[1,2] <- NA
  for(i in 2:tree$Nnode){
    table.tree[i,1] <- tree$edge[tree$edge[,2] == i+ntips ,1] - ntips
    table.tree[i,2] <- tree$edge.length[tree$edge[,2] == i+ntips] 
  }
  table.tree[,3] <- 1:tree$Nnode
  table.tree[,4] <- branching.times(tree)
  
  return(table.tree)
}


build_traits <- function(traits, tree){
  
  ## finding bifurcations
  tree  <- reorder(tree, "cladewise")
  ntips <- length(tree$tip.label)
  pp <- prop.part(tree)
  rk <- order(sapply(pp, length))
  
  ## Constructing hypervolumes (1:ntips) -> tips ((ntips+1):(ntips+nodes)) -> nodes. Hypervolumes of nodes are the unions of hypervolumes of descendants tips.
  unions <- list()
  sets <- list()
  sizes <- numeric(ntips + tree$Nnode)
  overlap <- numeric(tree$Nnode)
  asymmetry <- numeric(tree$Nnode)
  
  ## Tips hypervolumes
  for(i in 1:ntips){
    unions[[i]] <- hypervolume(traits[traits[,1] == tree$tip.label[i],-1])
    sizes[i] <- get_volume(unions[[i]])
  }
  
  ## Nodes hypervolumes
  for(i in rk + ntips){
    sets[[i-ntips]] <- hypervolume_set(unions[[tree$edge[tree$edge[,1] == i,2][1]]], unions[[tree$edge[tree$edge[,1] == i,2][2]]], check.memory = F)
    unions[[i]] <- sets[[i-ntips]]@HVList$Union
    sizes[i] <- get_volume(unions[[i]])
    overlap[i-ntips] <- get_volume(sets[[i - ntips]]@HVList$Intersection)
    asymmetry[i-ntips] <- abs(sizes[tree$edge[tree$edge[,1] == i,2][1]]-sizes[tree$edge[tree$edge[,1] == i,2][2]])
  }  
 
  out <- cbind(sizes[-(1:ntips)], asymmetry, overlap)
  dimnames(out) <- list(sprintf("n%s", 1:length(rk)), c("theta", "nu", "omega"))
  return(out)
}

joint_dist <- function (N1, N2){
  mu1 = N1[1]
  sig1 = N1[2]
  mu2 = N2[1]
  sig2 = N2[2]
  sig12 = (sig1*sig2)/(sig1+sig2)
  mu12 = (mu1/sig1 + mu2/sig2) * sig12
  return(c(mu12,sig12))
}

gibbs_sampler <- function(x, tab, pars, priors){
  theta0 <- pars[6:length(pars)]
  sig2 <- pars[3]
  
  ## Loop from root to most recent nodes
  for(i in 2:length(theta0)){
    # Get expected distribution fron descendants and ancestor
    desc.lik <- c(log(x[i,1]), sig2*tab[i,4])
    stem.lik <- c(theta0[tab[i,1]], sig2*tab[i,2])
    # Combine distributions to get estimates
    lik.val <- joint_dist(desc.lik, stem.lik)
    prior.prm <- priors[[5+i]](0)[[2]][[2]]
    # Calculate posterior
    post.prm <- joint_dist(prior.prm, lik.val)
    # Sample from posterior
    theta0[i] <- rnorm(1, post.prm[1], sqrt(post.prm[2]))
  }
  return(theta0)
}


## Make a mapping object to parse to mcmc ABM
make_ABM <- function(phy, traits, scale = F){
  
  ### scale height ###
  if(scale){
    t.len <- max(branching.times(phy))
    phy$edge.length <- phy$edge.length/t.len
  }
  
  ABM <- list()
  
  ### Global variables ###
  ABM$data$traits <- build_traits(traits, tree)
  ABM$data$n      <- length(phy$tip.label)
  ABM$data$tree   <- phy
  ABM$data$nodes  <- (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
  ABM$data$scale  <- scale
  ABM$data$tab <- build_table(phy)
  
  
  ### Likelihood parameters ###
  ABM$lik <- lik_ABM
  ABM$ws <- c(0.5, 0.5, 0.5, 0.1, 0.1, var(ABM$data$traits[,1]))
  ABM$init[1:5] <- c(runif(3, 0.1, 3), mean(ABM$data$traits[,2])/max(ABM$data$traits[,1]), mean(ABM$data$traits[,3])/max(ABM$data$traits[,1]))
  ABM$init[5+(1:ABM$data$tree$Nnode)] <- log(ABM$data$traits[,1])
  ABM$proposals <- list()
  ABM$priors <- list()
  for (i in 5+(1:ABM$data$tree$Nnode)) {
    ABM$priors[[i]] <- jive:::hpfun("Normal", c(ABM$init[i],.5))
  }
  ABM$priors[[1]] <- jive:::hpfun("Gamma", c(1.1, 5))
  ABM$priors[[2]] <- jive:::hpfun("Gamma", c(1.1, 5))
  ABM$priors[[3]] <- jive:::hpfun("Gamma", c(1.1, 5))
  ABM$priors[[4]] <- jive:::hpfun("Uniform", c(0,0.99))
  ABM$priors[[5]] <- jive:::hpfun("Uniform", c(0,1))
  ABM$proposals[[1]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[2]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[3]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[4]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[5]] <- jive:::proposal("slidingWinAbs")
  ABM$proposals[[6]] <- jive:::proposal("multiplierProposal")
  
  ## update frequencies : signu, sigomega, sigtheta, nu0, omega0, theta0, Gibbs
  ABM$update.freq <- c(0.16,0.16,0.16,0.16,0.16,0.16,1-(6*0.16))
  
  #### Prepare headers of log file ####
  
  ABM$header <- c("iter", "posterior", "sig2.nu", "sig2.omega","sig2.theta", "nu0", "omega0", sprintf("anc_state_node_n%s", 1:phy$Nnode), "log.lik", "acc")			
  
  return(ABM)
  
}


## mcmc chain
mcmc_ABM <- function(ABM, log.file = "my_ABM_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                     ngen = 5000000){
  
  # General syntax
  # 0 : used for claculations
  # 1 : not used for calculations
  
  ## checking
  if (length(ABM$update.freq) != 7 && !is.null(ABM$update.freq)) {
    stop("update.freq must contain 7 elements" )
  }
  
  ## initial conditions
  cat("setting initial conditions\n")
  
  # likelihood level
  pars0 <- ABM$init
  lik0 <- ABM$lik(ABM$data$tab, ABM$data$traits, pars0)
  prior0 <- unlist(mapply(do.call, ABM$priors, lapply(ABM$init, list))[1,])
  
  # mcmc parameters
  cat("generation\tposterior\n")
  cat(sprintf("%s\t", ABM$header), "\n", append = FALSE, file = log.file)
  update.freq <- cumsum(ABM$update.freq/sum(ABM$update.freq))
  proposals <- c(0,0,0,0,0,0,0) # 1st number: update ancestral states; 2nd: update sig2, 3rd: update rho, 4th: update delta
  proposals.accepted <- c(0,0,0,0,0,0,0)
  
  # posterior level
  post0 <- lik0 + sum(prior0)
  
  ## Start iterations
  for (i in 1:ngen) {
    
    # test whether we update ancestral states (r[1]), sig2 (r[2]), nu (r[3]) or omega (r[4])
    r <- which(runif(1) <= update.freq)[1]

    proposals[r] <- proposals[r] + 1
    lik1 <- lik0
    prior1 <- prior0  
    pars1 <- pars0
    
    if (r<7) # update
    { 
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[r]](i = pars0[r], d = ABM$ws[r], u)
      pars0[r] <- tmp$v
      hasting.ratio <- sum(tmp$lnHastingsRatio)
      gibbs <- F
    } else { # Gibbs
      pars0[-(1:5)] <- gibbs_sampler(ABM$data$traits, ABM$data$tab, pars0, ABM$priors)
      hasting.ratio <- 0
      gibbs <- T
    }
    
    # Posterior calculation
    prior0 <- unlist(mapply(do.call, ABM$priors, lapply(pars0, list))[1,])
    lik0 <- ABM$lik(ABM$data$tab, ABM$data$traits, pars0)
    post1 <- post0
    post0 <- lik0 + sum(prior0)
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0, prior0)))){
      pr <- -Inf
    } else {
      pr <- post0 - post1 + hasting.ratio
    }
    
    if (pr >= log(runif(1)) | gibbs) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
    } else # cancel changes
    {
      post0 <- post1
      prior0 <- prior1
      lik0 <- lik1
      pars0 <- pars1
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0) {
      cat(sprintf("%s\t", c(i, post0, pars0[1:5], exp(pars0[-(1:5)]), lik0,(sum(proposals.accepted)/i))),
          "\n", append=TRUE, file=log.file) 
    }
    
    # Print to screen
    if (i %% print.freq == 0) {
      cat(i,'\t',post0,'\n') 
    }
    
  } # end of for 
  
  # Calculate acceptance rate
  names(proposals) <- c("sig2.nu", "sig2.omega", "sig2.theta", "nu0", "omega0", "theta0", "thetai>0")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("sig2.nu", "sig2.omega", "sig2.theta", "nu0", "omega0", "theta0", "thetai>0")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
}
