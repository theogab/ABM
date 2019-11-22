## ABM simulations Output processing

lib <- c("jive", "ape", "coda")
sapply(lib, library, character.only = T)

NU <- c(0.1, 0.5, 0.9)
OM <- c(0.1, 0.5, 0.9)
SIG_RGE <- c(0.01, 0.1, 0.5)
SIG_MDP <- c(0.01, 0.1, 0.5)
TR <- 1:3
N <- c(20,100,500)
SIM <- 1:10

## Prediction of nu
jpeg("Figures/sim_out_nu_rge.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(nu), xlab = expression({sigma^2}[rge]))
    mtext(side = 1, at = c(2,6,10), text = SIG_RGE, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig_rge in SIG_RGE){
      lines(x1+c(-0.5,2.5), y = rep(nu,2), lty = 2, col = "#626066")
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          res <- do.call(rbind,lapply(SIG_MDP, function(sig_mdp){
            cat(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log\n", nu, om, sig_rge, sig_mdp, tr, sim))
            if(file.exists(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim))){
              read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim), header = T)
            }
          }))
          mean(res$nu, na.rm = T)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)

dev.off()


## Prediction of nu
jpeg("Figures/sim_out_nu_mdp.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(nu), xlab = expression({sigma^2}[mdp]))
    mtext(side = 1, at = c(2,6,10), text = SIG_MDP, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig_mdp in SIG_MDP){
      lines(x1+c(-0.5,2.5), y = rep(nu,2), lty = 2, col = "#626066")
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          res <- do.call(rbind,lapply(SIG_RGE, function(sig_rge){
            cat(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log\n", nu, om, sig_rge, sig_mdp, tr, sim))
            if(file.exists(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim))){
              read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim), header = T)
            }
          }))
          mean(res$nu, na.rm = T)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)

dev.off()

## Prediction of omega
jpeg("Figures/sim_out_omega_rge.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(omega), xlab = expression({sigma^2}[rge]))
    mtext(side = 1, at = c(2,6,10), text = SIG_RGE, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig_rge in SIG_RGE){
      lines(x1+c(-0.5,2.5), y = rep(om,2), lty = 2, col = "#626066")
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          res <- do.call(rbind,lapply(SIG_MDP, function(sig_mdp){
            cat(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log\n", nu, om, sig_rge, sig_mdp, tr, sim))
            if(file.exists(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim))){
              read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim), header = T)
            }
          }))
          mean(res$omega, na.rm = T)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)
dev.off()


## Prediction of omega
jpeg("Figures/sim_out_omega_mdp.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(omega), xlab = expression({sigma^2}[mdp]))
    mtext(side = 1, at = c(2,6,10), text = SIG_MDP, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig_mdp in SIG_MDP){
      lines(x1+c(-0.5,2.5), y = rep(om,2), lty = 2, col = "#626066")
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          res <- do.call(rbind,lapply(SIG_RGE, function(sig_rge){
            cat(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log\n", nu, om, sig_rge, sig_mdp, tr, sim))
            if(file.exists(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim))){
              read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig2_rge-%s_sig2_mdp-%s_tree-%s_sim-%s.log", nu, om, sig_rge, sig_mdp, tr, sim), header = T)
            }
          }))
          mean(res$omega, na.rm = T)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)
dev.off()


## Prediction of sigma^2
jpeg("Figures/sim_out_sig2.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(sigma^2), xlab = expression({sigma^2}[sim]))
    mtext(side = 1, at = c(2,6,10), text = SIG, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig in SIG){
      lines(x1+c(-0.5,2.5), y = rep(sig,2), lty = 2, col = "#626066")
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          res <- read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig-%s_tree-%s_sim-%s.log", nu, om, sig, tr, sim), header = T)
          mean(res$sigma.2)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)

dev.off()

## Predictions error ancestral states
jpeg("Figures/sim_out_anc_st.jpg", height = 20, width = 20, unit = "cm", res = 500)
par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
cols = c("#8cb369", "#f4a259", "#e71d36")
plots <- T
for(om in OM){
  plotn <- T
  for(nu in NU){
    plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = c(0,1), xaxt = "n", ylab = expression(R^2), xlab = expression({sigma^2}[sim]))
    mtext(side = 1, at = c(2,6,10), text = SIG, line = 1, cex = .6)
    if(plotn) mtext(bquote(omega == .(om)), 2, line = 5, cex = .8)
    if(plots) mtext(bquote(nu == .(nu)), 3, line = 3, cex = .8)
    x1 <- 1
    plotn <- F
    for(sig in SIG){
      x2 <- 0
      for(tr in TR){
        out <- sapply(SIM, function(sim){
          load(sprintf("Results/sim_tree_nu-%s_omega-%s_sig-%s_sim-%s.rda", nu, om, sig, sim))
          val <- sim_tree[[tr]][-(1:N[tr])]
          res <- read.table(sprintf("Results/ABM_mcmc_nu-%s_omega-%s_sig-%s_tree-%s_sim-%s.log", nu, om, sig, tr, sim), header = T)
          st <- colMeans(res[,grepl("anc_state", colnames(res))])
          1-sum((val-st)^2)/sum((val-mean(val))^2)
        })
        lines(x = rep(x1+x2, 2), y = quantile(out, c(.1, .9)), col = "#0f7173")
        points(x1+x2, mean(out), col = cols[x2+1], pch = 16)
        x2 <- x2 + 1
      }
      x1 <- x1 + 4
    }
  }
  plots <- F
}

mtext("N = ", side = 1, cex = 1, line = 5, at = -22)
legend(-20, -0.55, legend = N, col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)

dev.off()
