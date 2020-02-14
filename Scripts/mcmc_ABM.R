## Likelihood under an ABM model
# pars : vector with /[1] sig2_mdp /[2] sig2_rge /[3] nu /[4] omega
lik_ABM <- function(tree, logX, m, pars, Sx, Sm){
  sig2  <- pars[1:2]
  nu   <- pars[3]
  omega <- pars[4]
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]

  rg <- dnorm(logX[e2], rnge(omega, nu, logX[e1],Sx[e2]), sqrt(sig2[2] * tree$edge.length), log = T)
  mp <- dnorm(m[e2], midpt(omega, nu, logX[e1], m[e1], Sx[e2], Sm[e2]), sqrt(sig2[1] * tree$edge.length), log = T)
  
  #if(any(is.na(rg), is.na(mp))){
  #  cat("pars=", pars, "\n")
  #  cat("X=", X, "\n")
  #  cat("m=", m, "\n")
  #  cat("Sx=", Sx, "\n")
  #  cat("Sm=", Sm, "\n")
  #}
    
  lik <- sum(c(rg, mp))
  return(lik)
} 

# Range inheritance function
rnge <- function(omega, nu, logX0, Sx){
  X0 <- exp(logX0)
  return(log(X0/2*((omega+1)*(1-nu) + (1+Sx)*nu)))
}

anc_rnge<- function(nu,omega, logX1, logX2){
  X1 <- exp(logX1)
  X2 <- exp(logX2)
  return(log((X1 + X2)/((1-nu)*omega+1)))
}

# Midpoint inheritance function
midpt <- function(omega, nu, logX0, m0, Sx, Sm){
  rg <- rnge(omega, nu, logX0, Sx)
  X0 <- exp(logX0)
  return(m0 + Sm/2*(X0 - exp(rg)))
}

anc_midpt <- function(logX1, logX2, logXa, m1, m2){
  X1 <- exp(logX1)
  X2 <- exp(logX2)
  Xa <- exp(logXa)
  return((m1*(Xa-X2)+m2*(Xa-X1))/(2*Xa - X1 - X2))
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
  ABM$data$m <- traits[phy$tip.label,1]
  ABM$data$logX <- log(traits[phy$tip.label,2])
  ABM$data$n      <- length(phy$tip.label)
  ABM$data$tree   <- phy
  ABM$data$nodes  <- (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
  ABM$data$scale  <- scale

  ### Likelihood parameters ###
  ABM$lik <- lik_ABM
  # window sizes
  ABM$ws <- c(0.2, 1.2, 0.1, 0.1, 0.1, 0.1)
  # proposals
  ABM$proposals <- list()
  ABM$proposals[[1]] <- bite:::proposal("slidingWin")
  ABM$proposals[[2]] <- bite:::proposal("slidingWin")
  ABM$proposals[[3]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[4]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[5]] <- bite:::proposal("slidingWinAbs")
  ABM$proposals[[6]] <- bite:::proposal("slidingWinAbs")
  
  
  # initial conditions
  ABM$init$pars <- c(runif(2, 0.01, 0.99), 0.5, 0.5)
  ABM$init$Sm <- numeric(nrow(phy$edge))
  ABM$init$Sx <- numeric(nrow(phy$edge))
  d <- nrow(phy$edge)
  i <- 1
  while(d>0){
    if(any(ABM$init$Sm[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] %in% 0)){
      ABM$init$Sm[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] <- sample(c(-1,1), 2)
      ABM$init$Sx[phy$edge[phy$edge[,1] == phy$edge[i,1],2]] <- sample(c(-1,1), 2)
      d <- d - 2
    }
    i <- i + 1
  }
  
  ## Get realistic ancestral states depending on parameters
  pp <- prop.part(tree)
  I <- order(sapply(pp, length))
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  m0 <- c(ABM$data$m, rep(NA, ABM$data$tree$Nnode))
  logX0 <- c(ABM$data$logX, rep(NA, ABM$data$tree$Nnode))
  
  for(i in I){
    j <- i + ABM$data$n
    # value of descendants
    desc.logX <- logX0[pp[[i]]]
    logX0[j] <- anc_rnge(ABM$init$pars[3], ABM$init$pars[4], desc.logX[1], desc.logX[2])
    # value of descendants
    desc.m <- m0[pp[[i]]]
    m0[j] <- anc_midpt(desc.logX[1], desc.logX[2], logX0[j], desc.m[1], desc.m[2])
    
    ## replace descending tips by j in pp
    pp <- lapply(pp, function(d){
      if(any(pp[[i]] %in% d)){
        d <- d[-which(d %in% pp[[i]])]
        c(d, j)    
      } else {
        d
      }
    }) 
  }
  ABM$init$m <- m0[1:ABM$data$tree$Nnode + ABM$data$n]
  ABM$init$logX <- logX0[1:ABM$data$tree$Nnode + ABM$data$n]
  
  # default priors
  ABM$priors <- list(m = list(), logX = list(), pars = list())
  for (i in 1:phy$Nnode) {
    ABM$priors$m[[i]] <- bite:::hpfun("Normal", c(mean(ABM$data$m),10*var(ABM$data$m)))
    ABM$priors$logX[[i]] <- bite:::hpfun("Normal", c(log(mean(exp(ABM$data$logX))),10*var(ABM$data$logX)))
  }
  ABM$priors$pars[[1]] <- bite:::hpfun("Gamma", c(1.1, .5))
  ABM$priors$pars[[2]] <- bite:::hpfun("Gamma", c(1.1, .5))
  ABM$priors$pars[[3]] <- bite:::hpfun("Uniform", c(0,0.99))
  ABM$priors$pars[[4]] <- bite:::hpfun("Uniform", c(0,1))
  ABM$update.freq <- c(0.9, 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
  
  #### Prepare headers of log file ####
  
  ABM$header <- c("iter", "posterior", sprintf("m_anc_n%s", 1:phy$Nnode), sprintf("logX_anc_n%s", 1:phy$Nnode), sprintf("S_m%s", 1:(nrow(phy$edge))), sprintf("S_X%s", 1:(nrow(phy$edge))), "sigma^2_m", "sigma^2_logX", "nu", "omega", "log.lik", "acc")			
  
  return(ABM)
  
}

## join two normal distributions
jnorm <- function (means,vars){
  var = prod(vars)/sum(vars)
  mean = sum(means/vars) * var
  return(c(mean,var))
}

gibbs_move <- function(tree, n, m, logX, Sm, SX, pars, priors){
  
  sig2  <- pars[1:2]
  nu   <- pars[3]
  omega <- pars[4]
  
  pp <- prop.part(tree)
  I <- order(sapply(pp, length))
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  
  ## Conditional distributions of node's logX based on root state (P(logXi|logX0))
  anc.logX.var <- numeric(length(I))
  anc.logX.lik <- numeric(length(I))
  anc.logX.var[1] <- 0
  anc.logX.lik[1] <- logX[n+1]
  ## Start with the root nodes and go to the most recent nodes
  for(i in rev(I[-length(I)])){# ignore root state (length(I))
    j <- i + n
    anc.logX.var[i] <- anc.logX.var[e1[e2 == j]-n] + sig2[2] * tree$edge.length[e2 == j]
    anc.logX.lik[i] <- rnge(omega, nu, anc.logX.lik[e1[e2 == j]-n], SX[j])
  }
  
  ## Conditional distributions of node's logX based on descendants + combining and sampling
  samp.logX1 <- numeric(length(I))
  samp.logX2 <- numeric(length(I))
  ## Start with the most recent nodes and go to the root
  for(i in I[-length(I)]){# ignore root state (length(I))
    
    j <- i + n  # real node number
    # value of descendants
    desc.logX <- logX[pp[[i]]]
    # variance accumulated during evolutionary history
    desc.logX.var <- sig2[2] * tree$edge.length[e2 %in% pp[[i]]]
    
    # joint distribution of joint descendants anc ancestor giving conditional distribution of i
    join.logX1 <- jnorm(c(anc.logX.lik[i], desc.logX[1]), c(anc.logX.var[i], desc.logX.var[1]))
    join.logX2 <- jnorm(c(anc.logX.lik[i], desc.logX[2]), c(anc.logX.var[i], desc.logX.var[2]))
    
    # prior dstribution of the range
    prior.logX.anc <- priors$logX[[i]](1)[[2]][[2]]
    prior.logX.desc <- rnge(omega, nu, prior.logX.anc[1], SX[e2 %in% pp[[i]]])
    # posterior distribution of the range
    post.logX1 <- jnorm(c(prior.logX.desc[1], join.logX1[1]), c(prior.logX.anc[2], join.logX1[2]))
    post.logX2 <- jnorm(c(prior.logX.desc[2], join.logX2[1]), c(prior.logX.anc[2], join.logX2[2]))
    # sample logX from posterior
    samp.logX1[i] <- rnorm(1, post.logX1[1], sqrt(post.logX1[2]))
    samp.logX2[i] <- rnorm(1, post.logX2[1], sqrt(post.logX2[2]))
    # get logXa from descendants values
    logX[j] <- anc_rnge(nu, omega, samp.logX1[i], samp.logX2[i]) 
    
   ## replace descending tips by j in pp
    pp <- lapply(pp, function(d){
      if(any(pp[[i]] %in% d)){
        d <- d[-which(d %in% pp[[i]])]
        c(d, j)    
      } else {
        d
      }
    }) 
  }
  
  pp <- prop.part(tree)
  ## Conditional distributions of node's m based on root state (P(mi|m0))
  anc.m.var <- numeric(length(I))
  anc.m.lik <- numeric(length(I))
  anc.m.var[1] <- 0
  anc.m.lik[1] <- m[n+1]
  ## Start with the root nodes and go to the most recent nodes
  for(i in rev(I[-length(I)])){# ignore root state (length(I))
    j <- i + n
    anc.m.var[i] <- anc.m.var[e1[e2 == j]-n] + sig2[1] * tree$edge.length[e2 == j]
    anc.m.lik[i] <- midpt(omega, nu, logX[e1[e2 == j]], anc.m.lik[e1[e2 == j]-n], SX[j], Sm[j])
  }
  
  ## Conditional distributions of node's m based on descendants + combining and sampling
  ## Start with the most recent nodes and go to the root
  for(i in I[-length(I)]){# ignore root state (length(I))
    
    # value of descendants
    desc.m <- m[pp[[i]]]
    # variance accumulated during evolutionary history
    desc.m.var <- sig2[1] * tree$edge.length[e2 %in% pp[[i]]]
    
    # value of ancestor
    j <- i + n # real node number
    anc.m <- m[e1[e2 == j]]
    # variance accumulated during evolutionary history
    anc.m.var <- sig2[1] * tree$edge.length[e2 == j]
    
    # effect of asymmetric inheritance for the midpoint
    join.m1 <- jnorm(c(anc.m.lik[i], desc.m[1]), c(anc.m.var, desc.m.var[1]))
    join.m2 <- jnorm(c(anc.m.lik[i], desc.m[2]), c(anc.m.var, desc.m.var[2]))
    #plot(density(rnorm(1000, join.m1[1], sqrt(join.m1[2]))), xlim = c(0,20), col = "blue")
    #lines(density(rnorm(1000, join.m2[1], sqrt(join.m2[2]))), col = "red")
    # prior dstribution of the midpoint
    prior.m.anc <- priors$logX[[i]](1)[[2]][[2]]
    prior.m.desc <- midpt(omega, nu, logX[j], prior.m.anc[1], SX[e2 %in% pp[[i]]], Sm[e2 %in% pp[[i]]])
    # posterior distribution of the midpoint
    post.m1 <- jnorm(c(prior.m.desc[1], join.m1[1]), c(prior.m.anc[2], join.m1[2]))
    post.m2 <- jnorm(c(prior.m.desc[2], join.m2[1]), c(prior.m.anc[2], join.m2[2]))
    # sample midpoint from posterior
    samp.m1 <- rnorm(1, post.m1[1], sqrt(post.m1[2]))
    samp.m2 <- rnorm(1, post.m2[1], sqrt(post.m2[2]))
    
    # get logXa from descendants values
    m[j] <- anc_midpt(samp.logX1[i], samp.logX2[i], logX[j], samp.m1, samp.m2)
    
    ## replace descending tips by j in pp
    pp <- lapply(pp, function(d){
      if(any(pp[[i]] %in% d)){
        d <- d[-which(d %in% pp[[i]])]
        c(d, j)    
      } else {
        d
      }
    }) 
  }
  
  return(list(m0 = m, logX0 = logX))
}

## mcmc chain
mcmc_ABM <- function(ABM, log.file = "my_ABM_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                     ngen = 5000000, psml = F){
  
  # General syntax
  # 0 : used for calculations
  # 1 : not used for calculations
  
  ## initial conditions
  cat("setting initial conditions\n")
  
  # likelihood level
  m0 <- c(ABM$data$m, ABM$init$m)
  logX0 <- c(ABM$data$logX, ABM$init$logX)
  Sm0 <- ABM$init$Sm
  Sx0 <- ABM$init$Sx
  pars0 <- ABM$init$pars
  lik0 <- ABM$lik(ABM$data$tree, logX0, m0, pars0, Sx0, Sm0)
  prior.m0 <- unlist(mapply(do.call, ABM$priors$m, lapply(c(ABM$init$m), list))[1,])
  prior.logX0 <- unlist(mapply(do.call, ABM$priors$logX, lapply(c(ABM$init$logX), list))[1,])
  prior.pars0 <- unlist(mapply(do.call, ABM$priors$pars, lapply(c(ABM$init$pars), list))[1,])
  cat(lik0, "\n")
  
  # mcmc parameters
  cat("generation\tposterior\tlikelihood\n")
  cat(paste(ABM$header, collapse = "\t"), "\n", append = FALSE, file = log.file)
  update.freq <- cumsum(ABM$update.freq/sum(ABM$update.freq))
  proposals <- c(0,0,0,0,0,0,0,0,0)
  proposals.accepted <- c(0,0,0,0,0,0,0,0,0)
  
  # posterior level
  post0 <- lik0 + sum(prior.m0) + sum(prior.logX0) + sum(prior.pars0)
  
  ## Start iterations
  for (i in 1:ngen) {
    
    # test whether we update ancestral states (r[1]), sig2 (r[2]), nu (r[3]) or omega (r[4])
    r <- runif(1) <= update.freq
    r[-min(which(r))] <- F
    
    proposals[r] <- proposals[r] + 1
    lik1 <- lik0
    prior.m1 <- prior.m0 
    prior.logX1 <- prior.logX0 
    prior.pars1 <- prior.pars0 
    gibbs <- F
    
    if(r[1]) # Gibbs sampling
    {
      gibbs <- T
      new_st <- gibbs_move(tree = ABM$data$tree, n = ABM$data$n, m = m0, logX = logX0, Sm = Sm0, SX = Sx0, pars = pars0, priors = ABM$priors)
      m1 <- m0
      logX1 <- logX0
      m0 <- new_st$m0
      logX0 <- new_st$logX0
      hasting.ratio <- 0
    }
    
      
    if(r[2]) #update scenario midpoint
     {
      ## Change five scenarios randomly for midpoints
      anc <- sample(unique(ABM$data$tree$edge[,1]), 2)
      Sm1 <- Sm0
      Sm0[ABM$data$tree$edge[ABM$data$tree$edge[,1] %in% anc,2]] <- Sm0[ABM$data$tree$edge[ABM$data$tree$edge[,1] %in% anc,2]] * (-1)
      hasting.ratio <- 1
     }
    
    if(r[3]) #update scenario range
    {
      ## Change five scenarios randomly for ranges
      anc <- sample(unique(ABM$data$tree$edge[,1]), 2)
      Sx1 <- Sx0
      Sx0[ABM$data$tree$edge[ABM$data$tree$edge[,1] %in% anc,2]] <- Sx0[ABM$data$tree$edge[ABM$data$tree$edge[,1] %in% anc,2]] * (-1)
      hasting.ratio <- 1
    }
    
    if (r[4]) # update root midpoint
    { 
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      m1 <- m0
      tmp <- ABM$proposals[[1]](i = m0[ABM$data$n + 1], d = ABM$ws[1], u)
      m0[ABM$data$n + 1] <- tmp$v
      hasting.ratio <- sum(tmp$lnHastingsRatio)
    }
    
    if (r[5]) # update root range
    { 
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      logX1 <- logX0
      tmp <- ABM$proposals[[2]](i = logX0[ABM$data$n + 1], d = ABM$ws[2], u)
      logX0[ABM$data$n + 1] <- tmp$v
      hasting.ratio <- sum(tmp$lnHastingsRatio)
    }
    
    if (r[6]) # update sig2_m
    {
      pars1 <- pars0
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[3]](i = pars0[1], d = ABM$ws[3], u) #update with proposal function
      pars0[1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[7]) # update sig2_X
    {
      pars1 <- pars0
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[4]](i = pars0[2], d = ABM$ws[4], u) #update with proposal function
      pars0[2] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[8]) # update nu
    {
      pars1 <- pars0
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[5]](i = pars0[3], d = ABM$ws[5], u) #update with proposal function
      pars0[3] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if (r[9]) # update omega
    {
      pars1 <- pars0
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ABM$proposals[[6]](i = pars0[4], d = ABM$ws[6], u) #update with proposal function
      pars0[4] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    # Posterior calculation
    prior.m0 <- unlist(mapply(do.call, ABM$priors$m, lapply(m0[ABM$data$nodes], list))[1,])
    prior.logX0 <- unlist(mapply(do.call, ABM$priors$logX, lapply(logX0[ABM$data$nodes], list))[1,])
    prior.pars0 <- unlist(mapply(do.call, ABM$priors$pars, lapply(pars0, list))[1,])
    lik0 <- ABM$lik(ABM$data$tree, logX0, m0, pars0, Sx0, Sm0)
    post1 <- post0
    post0 <- lik0 + sum(prior.m0) + sum(prior.logX0) + sum(prior.pars0)
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0, prior.m0, prior.logX0, prior.pars0)))){
      pr <- -Inf
    } else {
     pr <- post0 - post1 + hasting.ratio
    }
    
    if(psml) f <- 0
    else f <- log(runif(1))
    if (pr >= f | gibbs)#log(runif(1))) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
    } else # cancel changes
    {
      post0 <- post1
      prior.m1 <- prior.m0 
      prior.logX1 <- prior.logX0 
      prior.pars1 <- prior.pars0 
      lik0 <- lik1
      
      if (r[2]){
        Sm0 <- Sm1
      } else if(r[3]){
        Sx0 <- Sx1
      } else if(r[4]){
        m0 <- m1
      } else if(r[5]){
        logX0 <- logX1
      } else {
        pars0 <- pars1
      }
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0) {
      cat(paste(c(i, post0, m0[ABM$data$nodes], logX0[ABM$data$nodes], Sm0[Sm0!=0], Sx0[Sx0!=0], pars0, lik0,(sum(proposals.accepted)/i)), collapse = "\t"),
          "\n", append=TRUE, file=log.file) 
    }
    
    # Print to screen
    if (i %% print.freq == 0) {
      #plot_states(ABM$data$tree, X0)
      #plot_states(ABM$data$tree, m0, "#5b7553")
      cat(i,'\t',post0, lik0, '\n') 
    }
    
  } # end of for 
  
  # Calculate acceptance rate
  names(proposals) <- c("Gibbs", "Sm", "Sx", "Ancestral midpoint","Ancestral range", "sigma^2_m","sigma^2_logX","nu","omega")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("Gibbs", "Sm", "Sx", "Ancestral midpoint","Ancestral range", "sigma^2_m","sigma^2_logX","nu","omega")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
}
