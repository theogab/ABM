### Code that will be added as example code
lib <- c("ape", "phytools", "OUwie")
sapply(lib, library, character.only = T)

# BM and VOU
#simulate a tree
N <- 50
tree <- pbtree(n = N, scale = 100)

#simulate species mean
the0.mbm = 1
sig.mbm = 0.5
m <- abs(fastBM(tree, a = the0.mbm, sig2 = sig.mbm))

#simulate species variance
alp.vou = c(0.5,0.5)
sig.vou = c(1, 1)
the.vou = c(1.0,0.5,0.5)
reg <- sample(1:2, 50, replace = T)
names(reg) <- tree$tip.label
tree_map <- make.simmap(tree, reg)
v <- OUwie.sim(tree_map, simmap.tree = T, alpha = alp.vou, sigma.sq = sig.vou, theta0 = the.vou[1], theta = the.vou[-1])
v[,2] <- abs(v[,2])

#simulate individual traits by species
spec.obs <- rpois(N, 20)
traits <- do.call(rbind, lapply(1:length(spec.obs), function(i){
  obs <- rep(NA, max(spec.obs))
  obs[sample(1:max(spec.obs), spec.obs[i])] <- abs(rnorm(spec.obs[i], m[i], v[i,2]))
  return(obs)
}))
rownames(traits) <- tree_map$tip.label
map <- matrix(tree_map$edge.length)
rownames(map) <- apply(tree_map$edge, 1, paste, collapse = ",")
tree_map$mapped.edge <- map

#treeOU1 <- tree_map
#traitsOU1 <- traits
#regimesOU1 <- tree_map$mapped.edge
#save(treeOU1, file = "data/treeOU1.rda")
#save(regimesOU1, file = "data/regimesOU1.rda")
#save(traitsOU1, file = "data/traitsOU1.rda")



## MBM and VOUM
#simulate a tree
N <- 50
tree <- pbtree(n = N, scale = 100)

#simulate species mean
the0.mbm = 1
sig.mbm = 0.5
m <- abs(fastBM(tree, a = the0.mbm, sig2 = sig.mbm))

#simulate species variance
alp.vou = c(0.5,0.5)
sig.vou = c(0.5, 1)
the.vou = c(1.0,0.5,2.0)
reg <- sample(1:2, 50, replace = T)
names(reg) <- tree$tip.label
tree_map <- make.simmap(tree, reg)
v <- OUwie.sim(tree_map, simmap.tree = T, alpha = alp.vou, sigma.sq = sig.vou, theta0 = the.vou[1], theta = the.vou[-1])
v[,2] <- abs(v[,2])

#simulate individual traits by species
spec.obs <- rpois(N, 20)
traits <- do.call(rbind, lapply(1:length(spec.obs), function(i){
  obs <- rep(NA, max(spec.obs))
  obs[sample(1:max(spec.obs), spec.obs[i])] <- abs(rnorm(spec.obs[i], m[i], v[i,2]))
  return(obs)
}))
rownames(traits) <- tree_map$tip.label

#treeOU2 <- tree_map
#traitsOU2 <- traits
#regimesOU2 <- tree_map$mapped.edge
#save(treeOU2, file = "data/treeOU2.rda")
#save(regimesOU2, file = "data/regimesOU2.rda")
#save(traitsOU2, file = "data/traitsOU2.rda")


#see what likOU does
pars <- c(alp.vou[1],sig.vou[1],the.vou[-1])
x <- apply(traits, 1, sd, na.rm = T)
jive:::lik_ou(pars, x, tree_map, tree_map$mapped.edge, root.station = T)


#see what lik_bm does
#BMS
pars <- c(sig.vou,the.vou[1])
x <- apply(traits, 1, sd, na.rm = T)
jive:::lik_bm(pars, x, tree_map, tree_map$mapped.edge, root.station = T)

#BM
pars <- c(sig.vou[1],the.vou[1])
x <- apply(traits, 1, sd, na.rm = T)
jive:::lik_bm(pars, x, tree_map, tree_map$mapped.edge)

res <- microbenchmark(v = v.reg(tree_map, tree_map$mapped.edge, 50, 1), w = w.reg.for(tree_map, tree_map$mapped.edge, 50, 1, .5), times =1000)
t1 <- w.reg(tree_map, tree_map$mapped.edge, 50, 1, .5)
t2 <- w.reg.for(tree_map, tree_map$mapped.edge, 50, 1, .5)


## testing the make_jive and control functions
jive <- make_jive(tree_map, traits, model.var = "OUM", model.mean = "BMM")


data(traitsOU2)
data(treeOU2)
data(regimesOU2)

## Run a simple MCMC chain
my.jive <- make_jive(treeOU2, traitsOU2,  model.var="OUM", model.mean="BMM", root.station = T)
mcmc_jive(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 



### Profiling
data(traitsOU2)
data(treeOU2)
data(regimesOU2)

## Run a simple MCMC chain
my.jive <- make_jive(treeOU2, traitsOU2,  model.var="OUM", model.mean="BMM", root.station = T)

setwd("/Users/admin/Dropbox/Documents/JIVE/jive-testing")
Rprof("mcmc_jive")
mcmc_jive(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 
Rprof()

summaryRprof("mcmc_jive")



## MBM and VOUM
#simulate a tree
N <- 500
tree <- pbtree(n = N, scale = 100)

#simulate species mean
the0.mbm = 1
sig.mbm = 0.5
m <- abs(fastBM(tree, a = the0.mbm, sig2 = sig.mbm))

#simulate species variance
alp.vou = c(0.5,0.5)
sig.vou = c(0.5, 1)
the.vou = c(1.0,0.5,2.0)
reg <- sample(1:2, N, replace = T)
names(reg) <- tree$tip.label
tree_map <- make.simmap(tree, reg)
tree <- tree_map
v <- OUwie.sim(tree_map, simmap.tree = T, alpha = alp.vou, sigma.sq = sig.vou, theta0 = the.vou[1], theta = the.vou[-1])
v[,2] <- abs(v[,2])

#simulate individual traits by species
spec.obs <- rpois(N, 20)
traits <- do.call(rbind, lapply(1:length(spec.obs), function(i){
  obs <- rep(NA, max(spec.obs))
  obs[sample(1:max(spec.obs), spec.obs[i])] <- abs(rnorm(spec.obs[i], m[i], v[i,2]))
  return(obs)
}))
rownames(traits) <- tree_map$tip.label

#save(tree, file = "tree.rda")
#save(traits, file = "traits.rda")

load("tree.rda")
load("traits.rda")
## Run a simple MCMC chain
my.jive <- make_jive(tree, traits,  model.var="OUM", model.mean="BMM", root.station = F)

setwd("/Users/admin/Dropbox/Documents/JIVE/jive-testing")
Rprof("mcmc_jive3")
mcmc_jive(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 
Rprof()

summaryRprof("mcmc_jive2")


