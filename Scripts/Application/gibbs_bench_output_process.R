## ABM simulations Output processing

lib <- c("bite", "ape", "coda","phytools", "sns")
sapply(lib, library, character.only = T)

res_BF <- matrix(0, nrow = 3*3*5*100, ncol = 38)
colnames(res_BF) <- c("n", "siglogv", "sigm", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                         "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                         "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                         "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","post_ESS", "rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                         "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")

get_res <- function(file, n, siglogv, sigm, nu, omega,sim,tree){
  out <- c(n, siglogv, sigm, nu, omega,sim)
  res <- read.table(file.path("/Volumes/Extreme\ SSD/ABM/Results/old_tests/", file), header = T, sep = "\t")
  res <- as.mcmc(res[-(1:500),!grepl("update", colnames(res))])
  res[,"nu"] <- res[,"nu"] * res[,"In"]
  res[,"omega"] <- res[,"omega"] * res[,"Io"]
  est <- colMeans(res)
  hpd <- HPDinterval(res)
  ess_col <- ess(res)
  out <- c(out, est[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega", "In", "Io")], sum(res[,"In"] == 0 & res[,"Io"] == 0)/nrow(res))
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),1])
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),2])
  out <- c(out, ess_col[c("post","logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega")])
  Sosim <- dat[[sim]]$So
  Snsim <- dat[[sim]]$Sn
  Soest <- est[grepl("So", names(est))]
  Snest <- est[grepl("Sn", names(est))]
  out <- c(out, sum(abs(Sosim - Soest)<0.1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<0.1)/length(Snest))
  out <- c(out, sum(abs(Sosim - Soest)<1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<1)/length(Snest))
  return(out)
}
sce <- c("SCI", "ADI", "ACI", "INT", "SDI")
i <- 1
for(fol in sce){
  files_bm <- list.files(file.path("/Volumes/Extreme\ SSD/ABM/Results/old_tests", fol))
  nsim_bm <- length(files_bm)
  cat(nsim_bm, "\n")
  if(nsim_bm > 0){
    old_n <- 0
    old_pars <- rep(0, 5)
    for(j in 1:nsim_bm){
      file_bm <- files_bm[j]
      sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][7])))
      n <- as.numeric(gsub("n", "", strsplit(file_bm, "_")[[1]][2]))
      nu <- as.numeric(gsub("nu", "", strsplit(file_bm, "_")[[1]][5]))
      omega <- as.numeric(gsub("omega-", "", strsplit(file_bm, "_")[[1]][6]))
      sigm <- siglogv <- as.numeric(gsub("sigm", "", strsplit(file_bm, "_")[[1]][3]))
      pars <- c(n,sigm,siglogv,nu,omega)
      if(n != old_n){
        tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
      }
      if(!all(old_pars == pars)){
        dat <- readRDS(sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_nu%s_omega%s.rds", n, sigm, siglogv, nu, omega))
      }
      old_n <- n
      old_pars <- pars
      cat("nu = ", nu, "omega = ", omega, "siglogv = ", sigm, "sigm = ", siglogv, "n = ", n, "sim = ", sim, "\n")
    if(all(res_BF[i,]==0)){
        try(res_BF[i,] <- get_res(file = file.path(fol,file_bm), n= n, siglogv = sigm, sigm = siglogv, nu = nu, omega =omega, sim = sim, tree = tree[[sim]]))
    }
    i = i + 1
    }
  }
}


saveRDS(res_BF, "Results/res_BF_v2.rds")

## With heterogeneous processes 
res_BF <- matrix(0, nrow = 1200, ncol = 35)
colnames(res_BF) <- c("prop", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")

get_res <- function(file, prop, nu, omega,sim,tree,dat){
  out <- c(prop, nu, omega,sim)
  res <- read.table(file.path("/Volumes/Extreme\ SSD/ABM/Results/hetero/", file), header = T, sep = "\t")
  res <- as.mcmc(res[-(1:500),!grepl("update", colnames(res))])
  res[,"nu"] <- res[,"nu"] * res[,"In"]
  res[,"omega"] <- res[,"omega"] * res[,"Io"]
  est <- colMeans(res)
  hpd <- HPDinterval(res)
  ess_col <- ess(res)
  out <- c(out, est[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega", "In", "Io")], sum(res[,"In"] == 0 & res[,"Io"] == 0)/nrow(res))
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),1])
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),2])
  out <- c(out, ess_col[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega")])
  obs_tips <- which(round(diag(vcv(tree)),9)==1)
  obs_nodes <- unique(tree$edge[tree$edge[,1] %in% mrca(tree)[obs_tips,obs_tips],1]  - length(tree$tip.label))
  Sosim <- dat$So[obs_nodes]
  Snsim <- dat$Sn[obs_nodes]
  Soest <- est[grepl("So", names(est))]
  Snest <- est[grepl("Sn", names(est))]
  out <- c(out, sum(abs(Sosim - Soest)<0.1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<0.1)/length(Snest))
  out <- c(out, sum(abs(Sosim - Soest)<1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<1)/length(Snest))
  return(out)
}

i <- 1
list_sc <- list(SDI = c(0,0.5), ADI = c(0.5, 0.5), ACI = c(0.5, 0), INT = c(0.2,0.2))
for(sc in 1:4){
  nu <- list_sc[[sc]][1]
  omega <- list_sc[[sc]][2]
  for(prop in c(0.2,0.5,0.8)){
    sig <- 0.1 
    n <- 20
    files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/hetero/", pattern = sprintf("^adim_h_prop%s_nu%s_omega-%s_", prop, nu, omega))
    nsim_bm <- length(files_bm)
    cat(nsim_bm, "\n")
    if(nsim_bm > 0){
      tree <- readRDS(sprintf("Results/Simulations/tree_prop%s_ran.rds", prop))
      dat <- readRDS(file = sprintf("Results/Simulations/sim_nu%s_omega%s_prop%s_ran.rds", nu, omega, prop))
      for(j in 1:nsim_bm){
        file_bm <- files_bm[j]
        sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][6])))
        cat("prop = ", prop, "nu = ", nu, "omega = ", omega, "sim = ", sim, "\n")
        if(all(res_BF[i,]==0)){
          try(res_BF[i,] <- get_res(file = file_bm, prop = prop, nu = nu, omega = omega, sim = sim, tree = tree[[sim]], dat = dat[[sim]]))
        }
        
        i = i + 1
      }
    }
  }
}

res_BF <- res_BF[res_BF[,1]!=0,]
saveRDS(res_BF, "Results/res_BF_hetero_sim.rds")

##### Output BF with indicators BM_vs_ABM resampled
res_BF <- matrix(0, nrow = 3*3*5*100, ncol = 37)
colnames(res_BF) <- c("n", "siglogv", "sigm", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")
i <- 1

for(nu in c(0, 0.2, 0.5)){
  for(omega in c(0, 0.2, 0.5)){
    #pdf(sprintf("Figures/sim_logX_nu%s_omega%s.pdf", nu, omega),height =10, width =4)
    #par(mfrow = c(3,2))
    for(sig in c(0.1, 0.5, 1)){
      for(n in c(20,50,100)){
        files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/", pattern = sprintf("^adim_resamp_n%s_sigm%s_siglogv%s_nu%s_omega-%s_", n, sig, sig, nu, omega))
        nsim_bm <- length(files_bm)
        cat(nsim_bm, "\n")
        tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n, sig, sig, nu, omega))
        dat <- readRDS(sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_nu%s_omega%s.rds", n, sig, sig, nu, omega))
        #x <- exp(dat[[1]]$logX)
        #names(x) <- c(tree[[1]]$tip.label, (n+1):(2*n-1))
        #phenogram(tree[[1]], x, main = sprintf("siglogX = %s", siglogX)) 
        #phenogram(tree[[1]], x[1:n], main = sprintf("siglogX = %s", siglogX)) 
        if(nsim_bm > 0){
          for(j in 1:nsim_bm){
            file_bm <- files_bm[j]
            sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][8])))
            cat("nu = ", nu, "omega = ", omega, "siglogv = ", sig, "sigm = ", sig, "n = ", n, "sim = ", sim, "\n")
            if(all(res_BF[i,]==0)){
              try(res_BF[i,] <- get_res(file_bm, n, sig, sig, nu, omega,sim, tree[[sim]]))
            }
            
            i = i + 1
          }
        }
      }
    }
    #dev.off()
  }
}

saveRDS(res_BF, "Results/res_BF_resamp_sim.rds")




##### Output BF with extinction
res_BF <- matrix(0, nrow = 1200, ncol = 35)
colnames(res_BF) <- c("mu", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")

get_res <- function(file,mu, nu, omega,sim,tree,dat){
  out <- c(mu, nu, omega,sim)
  res <- read.table(file.path("/Volumes/Extreme\ SSD/ABM/Results/ext/", file), header = T, sep = "\t")
  res <- as.mcmc(res[-(1:500),!grepl("update", colnames(res))])
  res[,"nu"] <- res[,"nu"] * res[,"In"]
  res[,"omega"] <- res[,"omega"] * res[,"Io"]
  est <- colMeans(res)
  hpd <- HPDinterval(res)
  ess_col <- ess(res)
  out <- c(out, est[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega", "In", "Io")], sum(res[,"In"] == 0 & res[,"Io"] == 0)/nrow(res))
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),1])
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),2])
  out <- c(out, ess_col[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega")])
  obs_tips <- which(round(diag(vcv(tree)),9)==1)
  obs_nodes <- unique(tree$edge[tree$edge[,1] %in% mrca(tree)[obs_tips,obs_tips],1]  - length(tree$tip.label))
  Sosim <- dat$So[obs_nodes]
  Snsim <- dat$Sn[obs_nodes]
  Soest <- est[grepl("So", names(est))]
  Snest <- est[grepl("Sn", names(est))]
  out <- c(out, sum(abs(Sosim - Soest)<0.1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<0.1)/length(Snest))
  out <- c(out, sum(abs(Sosim - Soest)<1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<1)/length(Snest))
  return(out)
}

i <- 1
list_sc <- list(SDI = c(0,0.5), ADI = c(0.5, 0.5), ACI = c(0.5, 0), INT = c(0.2,0.2))
for(sc in 1:4){
  nu <- list_sc[[sc]][1]
  omega <- list_sc[[sc]][2]
  for(mu in c(0.2,0.5,0.8)){
      sig <- 0.1 
      n <- 20
      files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/ext/", pattern = sprintf("^adim_ext_mu%s_nu%s_omega-%s_", mu, nu, omega))
      nsim_bm <- length(files_bm)
      cat(nsim_bm, "\n")
      if(nsim_bm > 0){
        tree <- read.tree(sprintf("Results/Simulations/tree_mu%s.tre", mu))
        dat <- readRDS(file = sprintf("Results/Simulations/sim_mu%s_nu%s_omega%s.rds", mu, nu, omega))
        for(j in 1:nsim_bm){
          file_bm <- files_bm[j]
          sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][6])))
          cat("mu = ", mu, "nu = ", nu, "omega = ", omega, "sim = ", sim, "\n")
          if(all(res_BF[i,]==0)){
            try(res_BF[i,] <- get_res(file = file_bm, mu = mu, nu = nu, omega = omega, sim = sim, tree = tree[[sim]], dat = dat[[sim]]))
          }
          
          i = i + 1
        }
      }
    }
}

res_BF <- res_BF[res_BF[,1]!=0,]
saveRDS(res_BF, "Results/res_BF_ext_sim.rds")


##### Output BF with incomplete taxon sampling
res_BF <- matrix(0, nrow = 1200, ncol = 35)
colnames(res_BF) <- c("n", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")

get_res <- function(file,n, nu, omega,sim,tree,dat){
  out <- c(n, nu, omega,sim)
  res <- read.table(file.path("/Volumes/Extreme\ SSD/ABM/Results/its/", file), header = T, sep = "\t")
  res <- as.mcmc(res[-(1:500),!grepl("update", colnames(res))])
  res[,"nu"] <- res[,"nu"] * res[,"In"]
  res[,"omega"] <- res[,"omega"] * res[,"Io"]
  est <- colMeans(res)
  hpd <- HPDinterval(res)
  ess_col <- ess(res)
  out <- c(out, est[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega", "In", "Io")], sum(res[,"In"] == 0 & res[,"Io"] == 0)/nrow(res))
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),1])
  out <- c(out, hpd[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega"),2])
  out <- c(out, ess_col[c("logv_root", "m_root", "sigma.2_logv", "sigma.2_m", "nu", "omega")])
  obs_tips <- which(round(diag(vcv(tree)),9)==1)
  obs_nodes <- unique(tree$edge[tree$edge[,1] %in% mrca(tree)[obs_tips,obs_tips],1]  - length(tree$tip.label))
  Sosim <- dat$So[obs_nodes]
  Snsim <- dat$Sn[obs_nodes]
  Soest <- est[grepl("So", names(est))]
  Snest <- est[grepl("Sn", names(est))]
  out <- c(out, sum(abs(Sosim - Soest)<0.1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<0.1)/length(Snest))
  out <- c(out, sum(abs(Sosim - Soest)<1)/length(Soest))
  out <- c(out, sum(abs(Snsim - Snest)<1)/length(Snest))
  return(out)
}

i <- 1
list_sc <- list(SDI = c(0,0.5), ADI = c(0.5, 0.5), ACI = c(0.5, 0), INT = c(0.2,0.2))
for(sc in 1:4){
  nu <- list_sc[[sc]][1]
  omega <- list_sc[[sc]][2]
    for(n in c(25,30,40)){
      files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/its/", pattern = sprintf("^adim_its_n%s_nu%s_omega-%s_", n, nu, omega))
      nsim_bm <- length(files_bm)
      cat(nsim_bm, "\n")
      tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n))
      dat <- readRDS(file = sprintf("Results/Simulations/sim_n%s_nu%s_omega%s.rds", n, nu, omega))
      if(nsim_bm > 0){
        for(j in 1:nsim_bm){
          file_bm <- files_bm[j]
          sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][6])))
          cat("n = ", n, "nu = ", nu, "omega = ", omega, "sim = ", sim, "\n")
          if(all(res_BF[i,]==0)){
            try(res_BF[i,] <- get_res(file = file_bm, n = n, nu = nu, omega = omega, sim = sim, tree = tree[[sim]], dat = dat[[sim]]))
          }
          
          i = i + 1
      }
    }
  }
}

res_BF <- res_BF[res_BF[,1]!=0,]
saveRDS(res_BF, "Results/res_BF_its_sim.rds")


### ##### Output BF with indicators BM_vs_ABM OU
res_BF <- matrix(0, nrow = 3*3*100, ncol = 37)
colnames(res_BF) <- c("n", "siglogv", "sigm", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")
i <- 1

for(sig in c(0.1, 0.5, 1)){
  for(n in c(20,50,100)){
    files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/", pattern = sprintf("^adim_OU_n%s_sigm%s_siglogv%s_", n, sig, sig))
    nsim_bm <- length(files_bm)
    cat(nsim_bm, "\n")
    tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n, sig, sig, 0, 1))
    dat <- readRDS(sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_OU.rds", n, sig, sig))
    if(nsim_bm > 0){
      for(j in 1:nsim_bm){
        file_bm <- files_bm[j]
        sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][6])))
        cat("nu = ", 0, "omega = ", 0, "siglogv = ", sig, "sigm = ", sig, "n = ", n, "sim = ", sim, "\n")
        if(all(res_BF[i,]==0)){
          try(res_BF[i,] <- get_res(file_bm, n, sig, sig, 0, 0,sim, tree[[sim]]))
        }
        
        i = i + 1
      }
    }
  }
}

saveRDS(res_BF, "Results/res_BF_OU.rds")



### ##### Output BF with indicators resampled OU
res_BF <- matrix(0, nrow = 3*3*100, ncol = 37)
colnames(res_BF) <- c("n", "siglogv", "sigm", "nu", "omega", "sim", "rootlogv_mean","rootm_mean", "siglogv_mean",
                      "sigm_mean", "nu_mean", "omega_mean",  "In", "Io", "Ino", "rootlogv_HPDmin", "rootm_HPDmin", "siglogv_HPDmin",
                      "sigm_HPDmin", "nu_HPDmin", "omega_HPDmin", "rootlogv_HPDmax", "rootm_HPDmax",
                      "siglogv_HPDmax", "sigm_HPDmax", "nu_HPDmax", "omega_HPDmax","rootlogv_ESS","rootm_ESS", "siglogv_ESS",
                      "sigm_ESS", "nu_ESS", "omega_ESS", "So_90","Sn_90", "So", "Sn")
i <- 1

for(sig in c(0.1, 0.5, 1)){
  for(n in c(20,50,100)){
    files_bm <- list.files("/Volumes/Extreme\ SSD/ABM/Results/", pattern = sprintf("^adim_OU_resamp_n%s_sigm%s_siglogv%s_", n, sig, sig))
    nsim_bm <- length(files_bm)
    cat(nsim_bm, "\n")
    tree <- read.tree(sprintf("Results/Simulations/tree_n%s.tre", n, sig, sig, 0, 1))
    dat <- readRDS(sprintf("Results/Simulations/sim_n%s_sigm%s_siglogv%s_OU.rds", n, sig, sig))
    if(nsim_bm > 0){
      for(j in 1:nsim_bm){
        file_bm <- files_bm[j]
        sim <- as.numeric(gsub(".log", "", gsub("sim", "", strsplit(file_bm, "_")[[1]][7])))
        cat("nu = ", 0, "omega = ", 0, "siglogv = ", sig, "sigm = ", sig, "n = ", n, "sim = ", sim, "\n")
        if(all(res_BF[i,]==0)){
          try(res_BF[i,] <- get_res(file_bm, n, sig, sig, 0, 0,sim, tree[[sim]]))
        }
        
        i = i + 1
      }
    }
  }
}

saveRDS(res_BF, "Results/res_BF_resamp_OU_sim.rds")

