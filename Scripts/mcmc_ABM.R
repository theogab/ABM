## Likelihood under an ABM model
# pars : vector with /[1] sig2_mdp /[2] sig2_rge /[3] nu /[4] omega
lik_ABM <- function(tree, logX, m, pars, SX, Sm, I){
  sig2  <- pars[1:2]
  nu   <- pars[3] * I[1]
  omega <- 1 + (pars[4] - 1) * I[2]
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  
  rg <- dnorm(logX[e2], rnge(omega, nu, logX[e1],SX[e2]), sqrt(sig2[2] * tree$edge.length), log = T)
  mp <- dnorm(m[e2], midpt(omega, nu, logX[e1], m[e1], SX[e2], Sm[e2]), sqrt(sig2[1] * tree$edge.length), log = T)
  
  lik <- sum(c(rg, mp))
  return(c(lik, rg, mp))
} 

## Range inheritance function
# Takes ancestral value and give descendant value for the range
rnge <- function(omega, nu, logX0, SX){
  X0 <- exp(logX0)
  return(log(X0/2*((omega+1)*(1-nu) + (1+SX)*nu)))
}

# Takes descendant and give ancestor value for the range
rev_rnge <- function(omega, nu, logX, SX){
  X <- exp(logX)
  return(log(X*2/((omega+1)*(1-nu) + (1+SX)*nu)))
}

# Expected range value for a node given the root and the sampled Scenarios for its ancestors
rnge_from_root <- function(omega, nu, logX0, SX){
  X0 <- exp(logX0)
  return(log(X0*prod(((omega+1)*(1-nu) + (1+SX)*nu)/2)))
}

# Midpoint inheritance function
midpt <- function(omega, nu, logX0, m0, SX, Sm){
  rg <- rnge(omega, nu, logX0, SX)
  X0 <- exp(logX0)
  return(m0 + Sm/2*(X0 - exp(rg)))
}

# Takes descendant and give ancestor value for the midpoint
rev_mdpt <- function(omega, nu, logX0, logX, m, SX, Sm){
  X0 <- exp(logX0)
  X <- exp(logX)
  return(m - Sm/2*(X0 - X))
}

midpt_from_root <- function(omega, nu, logX0, m0, SX, Sm){
  Dx <- -diff(c(exp(logX0),exp(sapply(1:length(SX), function(i) rnge_from_root(omega, nu, logX0, SX[1:i])))))
  return(m0 + sum(Sm/2*Dx))
}

## Make a mapping object to parse to mcmc ABM
make_ABM <- function(phy, traits, scale = F, switch = F){
  
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
  ABM$data$m <- traits[phy$tip.label,1]
  ABM$data$logX <- log(traits[phy$tip.label,2])
  ABM$data$n      <- length(phy$tip.label)
  ABM$data$tree   <- phy
  ABM$data$nodes  <- (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
  ABM$data$scale  <- scale
  ABM$data$switch  <- switch
  
  ### Likelihood parameters ###
  ABM$lik <- lik_ABM
  # window sizes
  ABM$ws <- c(0.001, 0.001, 0.1, 0.1)
  # proposals
  ABM$proposals <- list()
  ABM$proposals[[1]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[2]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[3]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[4]] <- bite:::proposal("slidingWinAbs")

  # initial conditions
  ABM$init$pars <- c(runif(2, 0.01, 0.99), .99, 0)
  ABM$init$Sm <- numeric(nrow(phy$edge))
  ABM$init$SX <- numeric(nrow(phy$edge))
  d <- nrow(phy$edge)
  i <- 1
  # random starting point for scenarios
  while(d>0){
    if(any(ABM$init$Sm[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] %in% 0)){
      ABM$init$Sm[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] <- sample(c(-1,1), 2)
      ABM$init$SX[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] <- sample(c(-1,1), 2)
      d <- d - 2
    }
    i <- i + 1
  }
  m0 <- c(ABM$data$m, mean(ABM$data$m), rep(NA, ABM$data$tree$Nnode-1))
  logX0 <- c(ABM$data$logX, mean(ABM$data$logX), rep(NA, ABM$data$tree$Nnode-1))
  ABM$init$I <- c(ifelse(switch,0,1), ifelse(switch,0,1))
  
  # default priors
  ABM$priors <- list(m = list(), logX = list(), pars = list(), Sm = list(), SX = list(), I = list())
  for (i in 1:phy$Nnode) {
    ABM$priors$m[[i]] <- bite:::hpfun("Normal", c(m0[ABM$data$n  + 1],10*sd(ABM$data$m)))
    ABM$priors$logX[[i]] <- bite:::hpfun("Normal", c(logX0[ABM$data$n  + 1],10*sd(ABM$data$logX)))
  }
  ABM$priors$pars[[1]] <- bite:::hpfun("Gamma", c(1.1, .1))
  ABM$priors$pars[[2]] <- bite:::hpfun("Gamma", c(1.1, .1))
  ABM$priors$pars[[3]] <- bite:::hpfun("Uniform", c(0,0.99))
  ABM$priors$pars[[4]] <- bite:::hpfun("Uniform", c(0,1))
  pp <- prop.part(ABM$data$tree)
  
  ABM$priors$I <- list(function(I){if(I == 0) list(log(0.95), 0.95) else list(ifelse(switch, log(0.05), log(1)), 0.95)}, function(I){if(I == 0) list(log(0.95), 0.95) else list(ifelse(switch, log(0.05), log(1)), 0.95)})
  
  ## Get realistic ancestral states depending on parameters
  
  st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = ABM$init$Sm, SX = ABM$init$SX, pars = ABM$init$pars, I = ABM$init$I, priors = ABM$priors, sc = T)
  m0 <- st$m
  logX0 <- st$logX
  ABM$init$m <- m0[1:ABM$data$tree$Nnode + ABM$data$n]
  ABM$init$logX <- logX0[1:ABM$data$tree$Nnode + ABM$data$n]
  ABM$init$SX <- st$SX
  ABM$init$Sm <- st$Sm
  
  ABM$update.freq <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4)
  
  #### Prepare headers of log file ####
  
  ABM$header <- c("iter", "posterior", "m_root", "logX_root", "sigma^2_m", "sigma^2_logX", "nu", "omega", "log.lik", sprintf("m_anc_n%s", ABM$data$n + 2:phy$Nnode), sprintf("logX_anc_n%s", ABM$data$n + 2:phy$Nnode), sprintf("S_m%s",  phy$edge[order(phy$edge[,2]),2]), sprintf("S_X%s",  phy$edge[order(phy$edge[,2]),2]), "I_nu", "I_omega", "acc", "update", "temperature")			
  
  return(ABM)
  
}

## join two normal distributions
jnorm <- function (means,vars){
  var = prod(vars)/sum(vars)
  mean = sum(means/vars) * var
  return(c(mean,var))
}

gibbs_move <- function(tree, n, m, logX, Sm, SX, pars, I, priors,sc){
  sig2  <- pars[1:2]
  nu   <- pars[3] * I[1]
  omega <- 1 + (pars[4] - 1) * I[2]
  ## Start with the most recent nodes and go to the root
  pp <- prop.part(tree)
  tip_to_root <- order(sapply(pp, length))
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  for(i in tip_to_root){# Do not sample the root states
    
    ### Range
    # variance accumulated during evolutionary history
    desc.logX.var <- sig2[2] * tree$edge.length[e2 %in% pp[[i]]]
    
    ## Sample SX (approximation not taking omega and nu)
    if(sc){
      P_logX1_lw_logX2 <- dnorm(diff(logX[pp[[i]]]), log((omega-omega*nu+1+nu)/(omega-omega*nu+1-nu)), sqrt(sum(desc.logX.var)))
      P_logX1_gt_logX2 <- dnorm(diff(logX[pp[[i]]]), log((omega-omega*nu+1-nu)/(omega-omega*nu+1+nu)), sqrt(sum(desc.logX.var)))
      if(P_logX1_lw_logX2 == P_logX1_gt_logX2) P_logX1_lw_logX2 <- P_logX1_gt_logX2 <- 0.5
      SX[pp[[i]]] <-sample(c(1,-1), 1, prob = c(P_logX1_gt_logX2, P_logX1_lw_logX2))*c(1,-1)
        # c(1,-1)[which.max(c(P_logX1_gt_logX2, P_logX1_lw_logX2))]*c(1,-1)
    }
    
    ## Conditional distributions on descendants
    # Expected ancestor distribution according to each descendant
    desc.logX <- rev_rnge(omega, nu, logX[pp[[i]]], SX[pp[[i]]])
    # joint ancestor distribution according to descendants
    join.logX <- jnorm(desc.logX, desc.logX.var)
      
    # prior distribution of the range
    prior.logX <- priors$logX[[i]](1)[[2]][[2]]
    # posterior distribution of the range
    post.logX <- jnorm(c(prior.logX[1], join.logX[1]), c(prior.logX[2], join.logX[2]))
    # Sample from posterior
    logX[i + n] <-  post.logX[1]
      #rnorm(1, post.logX[1], sqrt(post.logX[2]))
    
    ## Midpoint
    # variance accumulated during evolutionary history
    desc.m.var <- sig2[1] * tree$edge.length[e2 %in% pp[[i]]]
    ## Sample Sm (approximation not taking omega and nu)
    if(sc){
      P_m1_lw_m2 <- pnorm(diff(m[pp[[i]]]), (exp(logX[i + n])/2)*(1-omega+omega*nu),  sqrt(sum(desc.m.var)))
      P_m1_gt_m2 <- pnorm(diff(m[pp[[i]]]), (-exp(logX[i + n])/2)*(1-omega+omega*nu),  sqrt(sum(desc.m.var)))
      if(P_m1_lw_m2 == P_m1_gt_m2) P_m1_lw_m2 <- P_m1_gt_m2 <- 0.5
      Sm[pp[[i]]] <- sample(c(1,-1), 1, prob = c(P_m1_gt_m2, P_m1_lw_m2))*c(1,-1)
        #c(1,-1)[which.max(c(P_m1_gt_m2, P_m1_lw_m2))]*c(1,-1)
    }
    
    ## Conditional distributions on descendants
    # Expected ancestor distribution according to each descendant
    desc.m <- rev_mdpt(omega, nu, logX[i + n], logX[pp[[i]]], m[pp[[i]]], SX[pp[[i]]], Sm[pp[[i]]])
    # joint ancestor distribution according to descendants
    join.m <- jnorm(desc.m, desc.m.var)
    
    # prior distribution of the midpoint
    prior.m <- priors$m[[i]](1)[[2]][[2]]
    # posterior distribution of the range
    post.m <- jnorm(c(prior.m[1], join.m[1]), c(prior.m[2], join.m[2]))
    # Sample from posterior
    m[i + n] <- post.m[1]
      # rnorm(1, post.m[1], sqrt(post.m[2]))
    ## replace descending tips by j in pp
    pp <- lapply(pp, function(d){
      if(any(pp[[i]] %in% d)){
        d <- d[-which(d %in% pp[[i]])]
        c(d, i + n)    
      } else {
        d
      }
    }) 
  }
  return(list(m0 = m, logX0 = logX, SX0 = SX, Sm0 = Sm))
}

## mcmc chain
mcmc_ABM <- function(ABM, log.file = "my_ABM_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                     ngen = 5000000, psml = F, ncat = 1, beta.param = 0.3, burnin = 0){
  
  # General syntax
  # 0 : used for calculations
  # 1 : not used for calculations
  
  ## define the chain length for each category
  it <- ngen/ncat
  
  ## burnin
  if(burnin < 1) burnin <- burnin*it
  
  ## get the heating parameter for a chain - scaling classes
  if (ncat > 1) {
    beta.class <- bite:::heat_par(ncat, beta.param) 
  } else {
    beta.class <- 1
  }
  
  
  ## initial conditions
  cat("setting initial conditions\n")
  
  # likelihood level
  m0 <- c(ABM$data$m, ABM$init$m)
  logX0 <- c(ABM$data$logX, ABM$init$logX)
  Sm0 <- ABM$init$Sm
  SX0 <- ABM$init$SX
  pars0 <- ABM$init$pars
  I0 <- ABM$init$I
  lik0 <- ABM$lik(ABM$data$tree, logX0, m0, pars0, SX0, Sm0, I0)
  prior.m0 <- unlist(mapply(do.call, ABM$priors$m, lapply(c(ABM$init$m), list))[1,])
  prior.logX0 <- unlist(mapply(do.call, ABM$priors$logX, lapply(c(ABM$init$logX), list))[1,])
  prior.pars0 <- unlist(mapply(do.call, ABM$priors$pars, lapply(c(ABM$init$pars), list))[1,])
  prior.I0 <- unlist(mapply(do.call, ABM$priors$I, lapply(c(ABM$init$I), list))[1,])
  cat(lik0[1], "\n")
  
  # mcmc parameters
  cat("generation", "\t", "posterior", "\n")
  cat(paste(ABM$header, collapse = "\t"), "\n", append = FALSE, file = log.file)
  update.freq <- cumsum(ABM$update.freq/sum(ABM$update.freq))
  proposals <- c(0,0,0,0,0,0,0)
  proposals.accepted <- c(0,0,0,0,0,0,0)
  
  ## Path sampling
  it.beta <- 1
  bet <- beta.class[it.beta]
  if(ncat > 1) cat("beta = ", bet, "\n")
  
  # posterior level
  post0 <- lik0[1] * bet + sum(prior.m0) + sum(prior.logX0) + sum(prior.pars0) + sum(prior.I0)
  
  ## Start iterations
  for (i in 1:(it*ncat)) {
    
    # test whether we update ancestral states (r[1]), sig2 (r[2]), nu (r[3]) or omega (r[4])
    r <- runif(1) <= update.freq
    r[-min(which(r))] <- F
    
    proposals[r] <- proposals[r] + 1
    pars1 <- pars0
    m1 <- m0
    logX1 <- logX0
    SX1 <- SX0
    Sm1 <- Sm0
    I1 <- I0
    prior.m1 <- prior.m0 
    prior.logX1 <- prior.logX0 
    prior.pars1 <- prior.pars0
    prior.I1 <- prior.I0
    lik1 <- lik0
    post1 <- post0

    if (r[1]) # update sig2_m
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[1]](i = pars0[1], d = ABM$ws[1], u) #update with proposal function
      pars0[1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[2]) # update sig2_X
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[2]](i = pars0[2], d = ABM$ws[2], u) #update with proposal function
      pars0[2] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[3]) # update nu
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[3]](i = pars0[3], d = ABM$ws[3], u) #update with proposal function
      if(tmp$v < 1) pars0[3] <- tmp$v
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = SX0, pars = pars0, I = I0, priors = ABM$priors, sc = T)
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      SX0 <- new_st$SX0
      Sm0 <- new_st$Sm0
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if (r[4]) # update omega
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[4]](i = pars0[4], d = ABM$ws[4], u) #update with proposal function
      if(tmp$v <= 1) pars0[4] <- tmp$v
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = SX0, pars = pars0, I = I0, priors = ABM$priors, sc = T)
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      SX0 <- new_st$SX0
      Sm0 <- new_st$Sm0
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if(r[5]) # update Inu
    {
      I0[1] <- abs(I1[1] - 1)
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = SX0, pars = pars0, I = I0, priors = ABM$priors, sc = T)
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      SX0 <- new_st$SX0
      Sm0 <- new_st$Sm0
      hasting.ratio <- 0
    }
    
    if(r[6]) # update Iom
    {
      I0[2] <- abs(I1[2] - 1)
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = SX0, pars = pars0, I = I0, priors = ABM$priors, sc = T)
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      SX0 <- new_st$SX0
      Sm0 <- new_st$Sm0
      hasting.ratio <- 0
    }
    
    if(r[7]) # Gibbs move
    {
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = SX0, pars = pars0, I = I0, priors = ABM$priors, sc = T)
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      SX0 <- new_st$SX0
      Sm0 <- new_st$Sm0
      hasting.ratio <- 0
    }
    
    # Posterior calculation
    prior.m0 <- unlist(mapply(do.call, ABM$priors$m, lapply(m0[ABM$data$nodes], list))[1,])
    prior.logX0 <- unlist(mapply(do.call, ABM$priors$logX, lapply(logX0[ABM$data$nodes], list))[1,])
    prior.pars0 <- unlist(mapply(do.call, ABM$priors$pars, lapply(pars0, list))[1,])
    prior.I0 <- unlist(mapply(do.call, ABM$priors$I, lapply(I0, list))[1,])
    lik0 <- ABM$lik(ABM$data$tree, logX0, m0, pars0, SX0, Sm0, I0)
    post0 <- lik0[1] * bet + sum(prior.m0) + sum(prior.logX0) + sum(prior.pars0) + sum(prior.I0)
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0[1], prior.m0, prior.logX0, prior.pars0)))){
      pr <- -Inf
    } else {
      pr <- post0 - post1 + hasting.ratio
    }
    
    if(psml) f <- 0
    else f <- log(runif(1))
    if (pr >= f | r[7])#log(runif(1))) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
      mr <- rep(NA, length(ABM$data$nodes))
    } else # cancel changes
    {
      post0 <- post1
      prior.m1 <- prior.m0 
      prior.logX1 <- prior.logX0 
      prior.pars1 <- prior.pars0 
      lik0 <- lik1
      pars0 <- pars1
      I0 <- I1
      m0 <- m1
      logX0 <- logX1
      SX0 <- SX1
      Sm0 <- Sm1
     }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0 & i >= burnin) {
      cat(paste(c(i, post0, m0[ABM$data$n + 1], logX0[ABM$data$n + 1], pars0, lik0[1], m0[ABM$data$nodes[-1]], logX0[ABM$data$nodes[-1]], Sm0[Sm0!=0], SX0[SX0!=0], I0, (sum(proposals.accepted)/i), c("sig2_m", "sig2_logX", "nu", "omega", "I_nu", "I_omega")[r], bet), collapse = "\t"),
          "\n", append=TRUE, file=log.file) 
    }
    
    # Print to screen
    if (i %% print.freq == 0) {
      cat(i, post0, '\n') 
    }
    
    
    # change beta value if the length of the category is reach 
    if(i%%it == 0 & i < ngen){
      it.beta = it.beta+1
      bet <- beta.class[it.beta]
      cat("beta = ", bet, "\n")
    }
    
  } # end of for 
  
  # Calculate acceptance rate
  names(proposals) <- c("sigma^2_m","sigma^2_logX","nu","omega", "I_nu", "I_omega", "Gibbs")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("sigma^2_m","sigma^2_logX","nu","omega", "I_nu", "I_omega", "Gibbs")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
}
