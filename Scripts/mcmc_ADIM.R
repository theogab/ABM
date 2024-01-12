### ADIM Asymmetric and Displaced Inheritance Model
## Likelihood under an ADIM model
# pars : vector with /[1] sig2_mdp /[2] sig2_rge /[3] nu /[4] omega
lik_ADIM <- function(tree, logv, m, pars, Sn, So, I, maps){
  sig2 <- pars[1:2]
  nu   <- pars[3] * I[1]
  omega <- pars[4] * I[2]
  n <- length(tree$tip.label)
  a <- sqrt(2)*erfinv(0.95)
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  Sp <- (maps$desc.map[e2,]*maps$se.map)%*%rep(1,n-1)
  ## Likelihood of logv evolution
  Elogv <- logv[n + 1] + maps$anc.map[1:n,]%*%log(0.5*(2 + pmin(0, Sp* Sn[e1 - n]*nu - Sp *So[e1 - n]*omega*(2-nu)) + pmin(0, Sp*Sn[e1 - n]*nu + Sp* So[e1 - n]*omega*(2-nu))))
  Elogv_nodes <- logv[n + 1] - (maps$node.map%*%rep(1,n-1)-1)*log(2) + maps$anc.map[n+1:tree$Nnode,]%*%log(2 + pmin(0, Sp* Sn[e1 - n]*nu - Sp *So[e1 - n]*omega*(2-nu)) + pmin(0, Sp*Sn[e1 - n]*nu + Sp* So[e1 - n]*omega*(2-nu)))
  V <- vcv.phylo(tree)
  Vlogv <- V*sig2[2]
  detlogv <- determinant(Vlogv)$modulus[1]
  invlogv <- solve(Vlogv)
  logliklogv <- (-n/2 * log(2 * pi) - 1/2 * detlogv - 1/2 * (crossprod(logv[1:n] - Elogv,invlogv)%*%(logv[1:n] - Elogv)))
  ## Likelihood of m evolution
  Em <- m[n + 1] + a/2*maps$anc.map[1:n,]%*%(exp(Elogv_nodes[e1-n])*(pmin(0, Sp * Sn[e1 - n]*nu + Sp*So[e1 - n]*omega*(2-nu)) - pmin(0, Sp*Sn[e1 - n]*nu - Sp*So[e1 - n]*omega*(2-nu))))
  Em_nodes <- m[n + 1] + a/2*maps$anc.map[n+1:tree$Nnode,]%*%(exp(Elogv_nodes[e1-n])*(pmin(0, Sp * Sn[e1 - n]*nu + Sp*So[e1 - n]*omega*(2-nu)) - pmin(0, Sp*Sn[e1 - n]*nu - Sp*So[e1 - n]*omega*(2-nu))))
  Vm <- V*sig2[1]
  detm <- determinant(Vm)$modulus[1]
  invm <- solve(Vm)
  loglikm <- (-n/2 * log(2 * pi) - 1/2 * detm - 1/2 * (crossprod(m[1:n] - Em,invm)%*%(m[1:n] - Em)))
  return(list(loglik = loglikm + logliklogv, nodelik = c(loglikm , logliklogv)))
} 
## Variance inheritance function
# Takes ancestral value and give descendant value for the variance
# sp determines wether the descending value is given for the left descendant sp = -1 or the right descendant sp = 1
sigd <- function(n, o, siga, sn, so, sp){
  N <- sp*sn*n
  O <- sp*so*o*(2 - n)
  return(siga/2*(2 + pmin(0,N + O)  + pmin(0,N - O)))
}

# Mean inheritance function
# sp determines wether the descending value is given for the left descendant sp = -1 or the right descendant sp = 1
mud <- function(n, o, siga, mua, sn, so, sp){
  N <- sp*sn*n
  O <- sp*so*o*(2 - n)
  a <- sqrt(2)*erfinv(0.95)
  return(mua + a*siga/2*(pmin(0,N + O) - pmin(0, N - O)))
}


## Inverse error function
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)

## Make a mapping object to parse to mcmc ADIM
make_ADIM <- function(phy, traits, scale = F, switch = F){
  
  ### validity test ###
  if (geiger::name.check(phy, traits) != "OK") {
    stop("Species do not match in tree and traits")
  }
  
  
  ### scale height ###
  if(scale){
    t.len <- max(branching.times(phy))
    phy$edge.length <- phy$edge.length/t.len
  }
  
  ADIM <- list()
  
  ### Global variables ###
  ADIM$data$m <- traits[phy$tip.label,1]
  ADIM$data$logv <- log(traits[phy$tip.label,2])
  ADIM$data$n      <- length(phy$tip.label)
  ADIM$data$tree   <- phy
  ADIM$data$nodes  <- (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
  ADIM$data$scale  <- scale
  ADIM$data$switch  <- switch
  
  ### Likelihood parameters ###
  ADIM$lik <- lik_ADIM
  # window sizes
  ADIM$ws <- c(1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1)
  # proposals
  ADIM$proposals <- list()
  ADIM$proposals[[1]] <- bite:::proposal("slidingWinAbs")
  ADIM$proposals[[2]] <- bite:::proposal("slidingWinAbs")
  ADIM$proposals[[3]] <- bite:::proposal("slidingWinAbs")
  ADIM$proposals[[4]] <- bite:::proposal("slidingWinAbs")
  ADIM$proposals[[5]] <- bite:::proposal("slidingWin")
  ADIM$proposals[[6]] <- bite:::proposal("slidingWin")
  
  # initial conditions
  ADIM$init$pars <- c(runif(2, 1e-6, 1-1e-6), 0, 0)
  # random starting point for scenarios
  ADIM$init$So <- sample(c(-1,1), phy$Nnode, rep = T)
  ADIM$init$Sn <- sample(c(-1,1), phy$Nnode, rep = T)
  m0 <- c(ADIM$data$m, mean(ADIM$data$m), rep(NA, ADIM$data$tree$Nnode-1))
  logv0 <- c(ADIM$data$logv, mean(ADIM$data$logv), rep(NA, ADIM$data$tree$Nnode-1))
  ADIM$init$I <- c(ifelse(switch,0,1), ifelse(switch,0,1))
  
  # default priors
  ADIM$priors <- list(m = list(), logv = list(), pars = list(), So = list(), Sn = list(), I = list())
  pp <- prop.part(phy)
  for (i in 1:phy$Nnode) {
    ADIM$priors$m[[i]] <- bite:::hpfun("Normal", c(mean(ADIM$data$m[pp[[i]]]),2*sd(ADIM$data$m[pp[[i]]])))
    ADIM$priors$logv[[i]] <- bite:::hpfun("Normal", c(mean(ADIM$data$logv[pp[[i]]]),length(pp[[i]])*2*sd(ADIM$data$logv[pp[[i]]])))
  }
  ADIM$priors$pars[[1]] <- bite:::hpfun("Gamma", c(1.1, .1))
  ADIM$priors$pars[[2]] <- bite:::hpfun("Gamma", c(1.1, .1))
  ADIM$priors$pars[[3]] <- bite:::hpfun("Uniform", c(0,0.99))
  ADIM$priors$pars[[4]] <- bite:::hpfun("Uniform", c(0,0.99))
  pp <- prop.part(ADIM$data$tree)
  
  ADIM$priors$I <- list(function(I){if(I == 0) list(log(0.95), 0.95) else list(ifelse(switch, log(0.05), log(1)), 0.95)}, function(I){if(I == 0) list(log(0.95), 0.95) else list(ifelse(switch, log(0.05), log(1)), 0.95)})
  
  ### map for the tree Matrices 
  # map of tips descending from each node tips x nodes
  tip.map <- matrix(0, ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$n, dimnames = list(ADIM$data$tree$tip.label, ADIM$data$n + 1:ADIM$data$tree$Nnode))
  # map of tips descending from each node with indicators tips x nodes
  s.map <- matrix(0, ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$n, dimnames = list(ADIM$data$tree$tip.label, ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of edges between the root and each node + tip edge x node + tip
  anc.map <- matrix(0,nrow = ADIM$data$n + ADIM$data$tree$Nnode, ncol = length(ADIM$data$tree$edge.length), dimnames = list(c(ADIM$data$tree$tip.label, ADIM$data$n + 1:ADIM$data$tree$Nnode), 1:length(ADIM$data$tree$edge.length)))
  # map of edges descending from each node edges x nodes
  edge.map <- matrix(0,ncol = ADIM$data$tree$Nnode, nrow = length(ADIM$data$tree$edge.length), dimnames = list(1:length(ADIM$data$tree$edge.length), ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of edges descending from each node with indicators edges x nodes
  se.map <- matrix(0, ncol = ADIM$data$tree$Nnode, nrow = length(ADIM$data$tree$edge.length), dimnames = list(1:length(ADIM$data$tree$edge.length), ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of direct descendants from each node (tips + nodes) x nodes
  desc.map <- matrix(0, ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$tree$Nnode + ADIM$data$n, dimnames = list(c(ADIM$data$tree$tip.label, ADIM$data$n + 1:ADIM$data$tree$Nnode), ADIM$data$n + 1:ADIM$data$tree$Nnode))
  # map of nodes descending from each nodes nodes x nodes
  node.map <- matrix(0,ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$tree$Nnode, dimnames = list(ADIM$data$n +1:ADIM$data$tree$Nnode, ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of the number of internal nodes between each tip and each node tips x nodes
  ni.tip <- matrix(0,ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$n, dimnames = list(ADIM$data$tree$tip.label, ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of the number of internal nodes between each edge and each node edges x nodes
  ni.edge <- matrix(0,ncol = ADIM$data$tree$Nnode, nrow = length(ADIM$data$tree$edge.length), dimnames = list(1:length(ADIM$data$tree$edge.length), ADIM$data$n +1:ADIM$data$tree$Nnode))
  # map of the number of internal nodes between each node nodes x nodes
  ni.node <- matrix(0,ncol = ADIM$data$tree$Nnode, nrow = ADIM$data$tree$Nnode, dimnames = list(ADIM$data$n +1:ADIM$data$tree$Nnode, ADIM$data$n +1:ADIM$data$tree$Nnode))
  
  e1 <- ADIM$data$tree$edge[,1]
  e2 <- ADIM$data$tree$edge[,2]
  
  for(i in 1:ADIM$data$tree$Nnode){
    tip.map[pp[[i]],i] <- 1
    desc.map[e2[e1 == i + ADIM$data$n],i] <- 1
  }
  
  for(i in 1:length(ADIM$data$tree$edge.length)){
    anc.map[e2[i],] <- anc.map[e1[i],]
    anc.map[e2[i], i] <- 1
  }
  
  node.rank <- colSums(tip.map)
  for(i in order(node.rank)){
    # edges that descend from the descendants of node i
    desc.i <- e2[e1 == i + ADIM$data$n]
    if(any(desc.i > ADIM$data$n)){
      if(all(desc.i > ADIM$data$n)){
        s.map[,i] <- -1*abs(s.map[,desc.i[1] - ADIM$data$n]) + 1*abs(s.map[,desc.i[2] - ADIM$data$n])
        se.map[,i] <- -1*abs(se.map[,desc.i[1] - ADIM$data$n]) + 1*abs(se.map[,desc.i[2] - ADIM$data$n])
        s.map[is.na(s.map[,desc.i[1] - ADIM$data$n]) & is.na(s.map[,desc.i[2] - ADIM$data$n]),i] <- NA
        edge.map[,i] <- rowSums(edge.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n])
        node.map[,i] <- rowSums(node.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n])
      } else {
        s.map[,i] <- c(-1,1)[which(desc.i > ADIM$data$n)]*abs(s.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n])
        se.map[,i] <- c(-1,1)[which(desc.i > ADIM$data$n)]*abs(se.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n])
        s.map[desc.i[desc.i <= ADIM$data$n],i] <- c(-1,1)[which(desc.i <= ADIM$data$n)]
        edge.map[,i] <- edge.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n]
        node.map[,i] <- node.map[,desc.i[desc.i > ADIM$data$n] - ADIM$data$n]
      }
    } else {
      s.map[desc.i,i] <- c(-1,1)
    }
    node.map[c(i,desc.i[desc.i > ADIM$data$n] - ADIM$data$n),i] <- 1
    # edges that connect node i and its direct descendants
    se.map[which(e1 == i + ADIM$data$n), i] <- c(-1,1)
    edge.map[which(e1 == i + ADIM$data$n), i] <- 1
  }
  
  for(i in 1:ADIM$data$tree$Nnode){
    ni.tip[pp[[i]],i] <- rowSums(tip.map[pp[[i]], node.rank < node.rank[i]])  + 1
    ni.edge[edge.map[,i] == 1,i] <- rowSums(edge.map[edge.map[,i] == 1, node.rank < node.rank[i]])  + 1
    if(sum(node.map[,i]) > 0){
      if(sum(node.map[,i]) > 1){
        ni.node[node.map[,i] == 1,i] <- rowSums(node.map[node.map[,i] == 1, node.rank < node.rank[i]])  + 1
      } else {
        ni.node[node.map[,i] == 1,i] <- sum(node.map[node.map[,i] == 1, node.rank < node.rank[i]])  + 1
      }
    }
  }
  ADIM$data$maps$tip.map <- tip.map
  ADIM$data$maps$s.map <- s.map
  ADIM$data$maps$edge.map <- edge.map
  ADIM$data$maps$anc.map <- anc.map
  ADIM$data$maps$se.map <- se.map
  ADIM$data$maps$node.map <- node.map
  ADIM$data$maps$desc.map <- desc.map
  ADIM$data$maps$ni.tip <- ni.tip
  ADIM$data$maps$ni.edge <- ni.edge
  ADIM$data$maps$ni.node <- ni.node
  
  ## Get realistic ancestral states depending on parameters
  st <- gibbs_move(tree = ADIM$data$tree, n = ADIM$data$n, m = m0, logv = logv0, So = ADIM$init$So, Sn = ADIM$init$Sn, pars = ADIM$init$pars, I = ADIM$init$I, priors = ADIM$priors, maps = ADIM$data$maps)
  ADIM$init$m <- st$expec$m
  ADIM$init$logv <- st$expec$logv
  ADIM$init$Sn <- st$expec$Sn
  ADIM$init$So <- st$expec$So
  
  ADIM$update.freq <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1,0.1,0.2)
  
  #### Prepare headers of log file ####
  
  ADIM$header <- c("iter", "posterior", "m_root", "logv_root", "sigma^2_m", "sigma^2_logv", "nu", "omega", "log.lik", sprintf("m_anc_n%s", ADIM$data$n + 2:phy$Nnode), sprintf("logv_anc_n%s", ADIM$data$n + 2:phy$Nnode), sprintf("So%s",  ADIM$data$n + 1:phy$Nnode), sprintf("Sn%s",  ADIM$data$n + 1:phy$Nnode), "In", "Io", "acc", "update")			
  
  return(ADIM)
  
}

## join two normal distributions
jnorm <- function (means,vars){
  var = (vars[1,] * vars[2,])/(vars[1,] + vars[2,])
  mean = (means[1,]/vars[1,] + means[2,]/vars[2,]) * var
  return(rbind(mean,var))
}

## weight function for mu conditional distributions
muwei <- function(Sn, nu, omega, a){
  if(omega==0|nu==0) out <- rep(1/2, length(Sn))
  else out <- (pmin(0,-Sn*nu + omega*(2-nu)) - pmin(0,-Sn*nu - omega*(2-nu)))/(pmin(0,-Sn*nu + omega*(2-nu)) - pmin(0,-Sn*nu - omega*(2-nu)) + pmin(0, Sn*nu+omega*(2-nu))- pmin(0, Sn*nu - omega*(2-nu)))
  out[out==0] <- 1e-6
  out[out==1] <- 1 - 1e-6
  return(out)
}

gibbs_move <- function(tree, n, m, logv, So, Sn, pars, I, priors, maps){
  sig2  <- pars[1:2]
  nu   <- pars[3] * I[1]
  omega <- pars[4] * I[2]
  a <- sqrt(2)*erfinv(0.95)
  ## Start with the most recent nodes and go to the root
  pp <- prop.part(tree)
  e1 <- tree$edge[,1]
  e2 <- tree$edge[,2]
  Sp <- (maps$desc.map[e2,]*maps$se.map)%*%rep(1,n-1)
  #C <- log(2) -1/2*(log(2 + pmin(0, nu-omega*(2-nu)) - pmax(0,omega*(2-nu)+nu) + omega*(2-nu) + nu)  + log(2 + pmin(0, nu-omega*(2-nu)) - pmax(0,omega*(2-nu)+nu) + omega*(2-nu) - nu))
  C <- -2*log(2) + (log(2 + pmin(0, nu + omega*(2-nu)) + pmin(0, nu - omega*(2-nu))) + log(2 + pmin(0, -nu + omega*(2-nu)) + pmin(0, -nu - omega*(2-nu))))
  #D <- c(pmin(0,  nu - omega*(2-nu)) + pmax(0, nu + omega*(2-nu)) - nu)
  #K <- (omega*(2-nu) - min(0, nu - omega*(2-nu)) - max(0, nu + omega*(2-nu)) + nu)/(omega*(2-nu) + min(0, nu - omega*(2-nu)) + max(0, nu + omega*(2-nu)) - nu)
  #if(is.na(K)) K <- 1
  
  # variance accumulated during evolutionary history
  desc.logv.var <- sig2[2] * tree$edge.length
  desc.m.var <- sig2[1] * tree$edge.length
  ## Sample log(variance) parameters
  # E[logv] and V[logv] for every node
  anc.logv.mean <- rep(1,n)%*%((maps$tip.map*logv[1:n] - maps$ni.tip*C/2)*1/2^maps$ni.tip)
  anc.logv.var <- rep(1,length(tree$edge.length))%*%(maps$edge.map*desc.logv.var*1/4^(maps$ni.edge))
  # Priors for logv on every node
  logv.prior <- sapply(priors$logv, function(x) x(1)[[2]][[2]])
  # Posteriors for logv on every node
  logv.post <- jnorm(rbind(anc.logv.mean, logv.prior[1,]), rbind(anc.logv.var, logv.prior[2,]^2))
  # Gibbs sampling
  samp.logva <- rnorm(tree$Nnode, logv.post[1,], sqrt(logv.post[2,]))
  # E[Sn] and V[Sn] for every node
  logv.mean <- (c(logv[1:n], anc.logv.mean)[e2]%*%(maps$desc.map[e2,]*maps$se.map))[1,]
  logv.var <- (c(desc.logv.var[e2 <= n], anc.logv.var)%*%maps$desc.map)[1,]
  P.sn <- pnorm(0, logv.mean, sqrt(logv.var), lower.tail = F)
  # Gibbs sampling
  samp.sn <- c(1,-1)[rbinom(tree$Nnode, 1, 1 - P.sn) + 1]
  E.sn <- c(1,-1)[apply(cbind(P.sn,rep(0.5,n-1)), 1, which.max)]
  
  ## Sample mean parameters given logv and Sn samples
  # E[m] and V[m] for every node
  we.edge <- muwei(Sp*samp.sn[e1-n], nu, omega, a)
  we.map <- maps$anc.map[1:n,]%*%diag(we.edge)
  we.map[we.map>0] <- log(we.map[we.map>0])
  # Products of branch weights per descending node
  we.prod <- exp(we.map%*%maps$edge.map)
  anc.m.mean <- rep(1,n)%*%((maps$tip.map*m[1:n])*we.prod)
  #anc.m.mean <-  rep(1,n)%*%(maps$tip.map*m[1:n]*K^((maps$ni.tip+maps$s.map%*%(maps$node.map*samp.sn))/2)/(1+K)^maps$ni.tip)
  anc.m.var <- rep(1,length(tree$edge.length))%*%(maps$edge.map*desc.m.var*we.edge^2)
  #anc.m.var <- rep(1,length(tree$edge.length))%*%((maps$edge.map*desc.m.var*K^(maps$ni.edge+maps$se.map%*%(maps$node.map*samp.sn)))/(1+K)^(2*maps$ni.edge))
  m.mean <- (c(m[1:n], anc.m.mean)[e2]%*%(maps$desc.map[e2,]*maps$se.map))[1,]
  m.var <- (c(desc.m.var[e2 <= n], anc.m.var)%*%maps$desc.map)[1,]
  # Priors for m on every node
  m.prior <- sapply(priors$m, function(x) x(1)[[2]][[2]])
  # Posteriors for m on every node
  m.post <- jnorm(rbind(anc.m.mean, m.prior[1,]), rbind(anc.m.var, m.prior[2,]^2))
  # Gibbs sampling
  samp.mua <- rnorm(tree$Nnode, m.post[1,], sqrt(m.post[2,]))
  # E[So] and V[So] for every node
  P.so <- pnorm(0, m.mean, sqrt(m.var), lower.tail = F)
  # Gibbs sampling
  E.so <- c(1,-1)[apply(cbind(P.so,rep(0.5,n-1)), 1, which.max)]
  samp.so <- c(1,-1)[rbinom(tree$Nnode, 1, 1 - P.so) + 1]
  
  return(list(gibbs = list(m = samp.mua, logv = samp.logva, Sn = samp.sn, So = samp.so),
              expec = list(m = anc.m.mean, logv = anc.logv.mean, Sn = E.sn, So = E.so)))
}

## mcmc chain
mcmc_ADIM <- function(ADIM, log.file = "my_ADIM_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                     ngen = 5000000, ncat = 1, beta.param = 0.3, burnin = 0){
  
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
  m0 <- c(ADIM$data$m, ADIM$init$m)
  logv0 <- c(ADIM$data$logv, ADIM$init$logv)
  So0 <- ADIM$init$So
  Sn0 <- ADIM$init$Sn
  pars0 <- ADIM$init$pars
  I0 <- ADIM$init$I
  lik0 <- ADIM$lik(tree = ADIM$data$tree, logv = logv0, m =  m0, pars =  pars0, Sn = Sn0, So = So0, I = I0, maps = ADIM$data$maps)
  prior.m0 <- unlist(mapply(do.call, ADIM$priors$m, lapply(c(ADIM$init$m), list))[1,])
  prior.logv0 <- unlist(mapply(do.call, ADIM$priors$logv, lapply(c(ADIM$init$logv), list))[1,])
  prior.pars0 <- unlist(mapply(do.call, ADIM$priors$pars, lapply(c(ADIM$init$pars), list))[1,])
  prior.I0 <- unlist(mapply(do.call, ADIM$priors$I, lapply(c(ADIM$init$I), list))[1,])
  cat(lik0[[1]], "\n")
  
  # mcmc parameters
  cat("generation", "\t", "posterior", "\n")
  cat(paste(ADIM$header, collapse = "\t"), "\n", append = FALSE, file = log.file)
  update.freq <- cumsum(ADIM$update.freq/sum(ADIM$update.freq))
  proposals <- c(0,0,0,0,0,0,0,0,0)
  proposals.accepted <- c(0,0,0,0,0,0,0,0,0)
  
  ## Path sampling
  it.beta <- 1
  bet <- beta.class[it.beta]
  if(ncat > 1) cat("beta = ", bet, "\n")
  
  # posterior level
  post0 <- lik0[[1]] * bet + prior.m0[1] + prior.logv0[1] + sum(prior.pars0) + sum(prior.I0)
  
  ## Start iterations
  for (i in 1:(it*ncat)) {
    
    # test whether we update ancestral states (r[1]), sig (r[2]), nu (r[3]) or omega (r[4])
    r <- runif(1) <= update.freq
    r[-min(which(r))] <- F
    
    proposals[r] <- proposals[r] + 1
    pars1 <- pars0
    m1 <- m0
    logv1 <- logv0
    Sn1 <- Sn0
    So1 <- So0
    I1 <- I0
    prior.m1 <- prior.m0 
    prior.logv1 <- prior.logv0 
    prior.pars1 <- prior.pars0
    prior.I1 <- prior.I0
    lik1 <- lik0
    post1 <- post0
    
    if (r[1]) # update sig2_m
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[1]](i = pars0[1], d = ADIM$ws[1], u) #update with proposal function
      pars0[1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[2]) # update sig2_logv
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[2]](i = pars0[2], d = ADIM$ws[2], u) #update with proposal function
      pars0[2] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[3]) # update nu
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[3]](i = pars0[3], d = ADIM$ws[3], u) #update with proposal function
      if(tmp$v < 1) pars0[3] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if (r[4]) # update omega
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[4]](i = pars0[4], d = ADIM$ws[4], u) #update with proposal function
      if(tmp$v < 1) pars0[4] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    if (r[5]) # update rootm
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[5]](i = m0[ADIM$data$n+1], d = ADIM$ws[5], u) #update with proposal function
      m0[ADIM$data$n+1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[6]) # update rootlogv
    {
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      tmp <- ADIM$proposals[[6]](i = logv0[ADIM$data$n +1], d = ADIM$ws[6], u) #update with proposal function
      logv0[ADIM$data$n+1] <- tmp$v
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    
    if(r[7]) # update Inu
    {
      I0[1] <- abs(I1[1] - 1)
      hasting.ratio <- 0
    }
    
    if(r[8]) # update Iom
    {
      I0[2] <- abs(I1[2] - 1)
      hasting.ratio <- 0
    }
    
    new_st <- gibbs_move(tree = ADIM$data$tree, n = ADIM$data$n, m = m0, logv = logv0, So = So0, Sn = Sn0, pars = pars0, I = I0, priors = ADIM$priors, maps = ADIM$data$maps)
    
    if(r[9]){ # Gibbs
      m0 <- c(m1[1:(ADIM$data$n+1)], new_st$gibbs$m[-1])
      logv0 <- c(logv1[1:(ADIM$data$n+1)], new_st$gibbs$logv[-1])
      Sn0 <- new_st$gibbs$Sn
      So0 <- new_st$gibbs$So
      hasting.ratio <- 0
    } else { # expectation
      logv0 <- c(logv0[1:(ADIM$data$n+1)], new_st$expec$logv[-1])
      m0 <- c(m0[1:(ADIM$data$n+1)], new_st$expec$m[-1])
      Sn0 <- new_st$expec$Sn
      So0 <- new_st$expec$So
    }
      
      
    
    # Posterior calculation
    prior.m0 <- unlist(mapply(do.call, ADIM$priors$m, lapply(m0[ADIM$data$nodes], list))[1,])
    prior.logv0 <- unlist(mapply(do.call, ADIM$priors$logv, lapply(logv0[ADIM$data$nodes], list))[1,])
    prior.pars0 <- unlist(mapply(do.call, ADIM$priors$pars, lapply(pars0, list))[1,])
    prior.I0 <- unlist(mapply(do.call, ADIM$priors$I, lapply(I0, list))[1,])
    lik0 <- ADIM$lik(ADIM$data$tree, logv0, m0, pars0, Sn0, So0, I0, ADIM$data$maps)
    post0 <- lik0[[1]] * bet + prior.m0[1] + prior.logv0[1] + sum(prior.pars0) + sum(prior.I0)
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0[[1]], prior.m0, prior.logv0, prior.pars0)))){
      pr <- -Inf
    } else {
      pr <- post0 - post1 + hasting.ratio
    }
    
    f <- log(runif(1))
    if (pr >= f | r[9])#log(runif(1))) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
      mr <- rep(NA, length(ADIM$data$nodes))
      if(any(r[c(1:4,7,8)])){
        m0 <- m1
        logv0 <- logv1
        Sn0 <- Sn1
        So0 <- So1 
      }
    } else # cancel changes
    {
      post0 <- post1
      prior.m1 <- prior.m0 
      prior.logv1 <- prior.logv0 
      prior.pars1 <- prior.pars0 
      lik0 <- lik1
      pars0 <- pars1
      I0 <- I1
      m0 <- m1
      logv0 <- logv1
      Sn0 <- Sn1
      So0 <- So1
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0 & i >= burnin) {
      cat(paste(c(i, post0, m0[ADIM$data$n + 1], logv0[ADIM$data$n + 1], pars0, lik0[[1]], m0[ADIM$data$nodes[-1]], logv0[ADIM$data$nodes[-1]], So0, Sn0, I0, (sum(proposals.accepted)/i), c("sig2_m", "sig2_logv", "nu", "omega", "I_nu", "I_omega")[r]), collapse = "\t"),
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
  names(proposals) <- c("sigma^2_m","sigma^2_logv","nu","omega","rootm", "rootlogv", "I_nu", "I_omega", "Gibbs")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("sigma^2_m","sigma^2_logv","nu","omega","rootm", "rootlogv", "I_nu", "I_omega", "Gibbs")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
}

