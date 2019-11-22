## Test fct° mcmc_ABM
lib <- c("ape", "phytools", "bite", "nloptr")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ABM.R")


ex_lik_cor <- function(vari, range_i, varj, range_j){
  
  rangei <- seq(range_i[1],range_i[2], length.out = 101)
  rangej <- seq(range_j[1],range_j[2], length.out = 101)
  out <- matrix(NA, nrow = length(rangei), ncol = length(rangej))
  for(i in 1:length(rangei)){
    for(j in 1:length(rangej)){
      x0 <- x
      pars0 <- pars
      
      if(vari == "node_mdp"){
        x0[5,1] <- x[5,1] + rangei[i]
        tvi <- x[5,1]
      }
      if(vari == "ancestor_mdp"){
        x0[4,1] <- x[4,1] + rangei[i]
        tvi <- x[4,1]
      }
      if(vari == "node_rge"){
        x0[5,2] <- x[5,2] + rangei[i]
        tvi <- x[5,2]
      }
      if(vari == "ancestor_rge"){
        x0[4,2] <- x[4,2] + rangei[i]
        tvi <- x[4,2]
      }
      if(vari == "sig_mdp"){
        pars0[1] <- pars[1] + rangei[i]
        tvi <- pars[1]
      }
      if(vari == "sig_rge"){
        pars0[2] <- pars[2] + rangei[i]
        tvi <- pars[2]
      }
      if(vari == "nu"){
        pars0[3] <- pars[3] + rangei[i]
        tvi <- pars[3]
      }
      if(vari == "omega"){
        pars0[4] <- pars[4] + rangei[i]
        tvi <- pars[4]
      }
      if(varj == "node_mdp"){
        x0[5,1] <- x[5,1] + rangej[j]
        tvj <- x[5,1]
      }
      if(varj == "ancestor_mdp"){
        x0[4,1] <- x[4,1] + rangej[j]
        tvj <- x[4,1]
      }
      if(varj == "node_rge"){
        x0[5,2] <- x[5,2] + rangej[j]
        tvj <- x[5,2]
      }
      if(varj == "ancestor_rge"){
        x0[4,2] <- x[4,2] + rangej[j]
        tvj <- x[4,2]
      }
      if(varj == "sig_mdp"){
        pars0[1] <- pars[1] + rangej[j]
        tvj <- pars[1]
      }
      if(varj == "sig_rge"){
        pars0[2] <- pars[2] + rangej[j]
        tvj <- pars[2]
      }
      if(varj == "nu"){
        pars0[3] <- pars[3] + rangej[j]
        tvj <- pars[3]
      }
      if(varj == "omega"){
        pars0[4] <- pars[4] + rangej[j]
        tvj <- pars[4]
      }
      
      out[i,j] <- lik_ABM(tab, x0, pars0)
    }
  }
  
  out[out <= 7] <- NA
  out2 <- out - min(out, na.rm = T)
  
  colVec <- rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4",
                  "#88ddaa", "#66c2a5", "#37b69b", "#3288bd", "#71acbc")) #c("#E83509","#F9ED11","#B4D92A")
  jet.colors <- colorRampPalette(colVec)
  colour <- jet.colors(101)
  plot(1:length(rangei), col = "white", xlab = vari, ylab = varj, xaxt = "n", yaxt = "n")
  axis(1, at = seq(1, length(rangei), length.out = 5), labels = round(tvi + seq(range_i[1],range_i[2],length.out = 5), 2))
  axis(2, at = seq(1, length(rangei), length.out = 5), labels = round(tvj + seq(range_j[1],range_j[2],length.out = 5), 2))
  for(i in 1:length(rangei)){
    for(j in 1:length(rangej)){
      points(i,j, pch = 15, col = colour[round(out2[i,j]*100/max(out2, na.rm = T),0)+1])
    }
  }
  points(1:length(rangei), rep(length(rangej)*1.1, length(rangei)) , pch = 15, xpd = T, col = colour)
  axis(3, at = seq(1, length(rangei), length.out = 4), labels = round(seq(min(out, na.rm = T), max(out, na.rm = T), length.out = 4), 2), line = 1)
}


# lik fct°
for(n in c(0.1, 0.5, 0.9)){
  for(o in c(0.1, 0.5, 0.9)){
    pars <- c(0.001, 0.001, n, o)
    tree <- pbtree(n = 3)
    tab <- build_table(tree)
    sig2 <- pars[1:2]
    sig2_mdp <- pars[1]
    sig2_rge  <- pars[2]
    nu <- pars[3]
    omega <- pars[4]
    
    sim1 <- sim_ABM(tree, nu = nu, omega = omega, pars_rge = c(sig2_rge, 2, 0, Inf), pars_mdp = c(sig2_mdp, 10, -Inf, Inf), nsim = 1)
    sim <- sim1[[1]][1:3,]
    x <- sim1[[1]]
    rownames(sim) <- tree$tip.label
    
    mdp <- c("ancestor_mdp", "node_mdp", "sig_mdp", "nu", "omega")
    range_mdp <- list(c(-0.5,0.5), c(-0.5,0.5), c(0, 0.2), c(-0.05,0.05), c(-0.05,0.05))
    rge <- c("ancestor_rge", "node_rge", "sig_rge", "nu", "omega")
    range_rge <- list(c(-0.5,0.5), c(-0.5,0.5), c(0, 0.2), c(-0.05,0.05), c(-0.05,0.05))
    
    jpeg(sprintf("Figures/Likelihood_space_om-%s_nu-%s_zoom.jpg", o, n), height = 30, width = 30, unit = "cm", res = 400)
    par(mfrow = c(5,5), mar = c(4,4,4,1))
    for(k in 1:5){
      cat("k =", k, "\n")
      for(l in 1:5){
        cat(l, "\n")
        if(k < l){
          if(mdp[k] == mdp[l]){
            plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
          } else {
            ex_lik_cor(mdp[k], range_mdp[[k]], mdp[l], range_mdp[[l]])
          }
        } else {
          if(rge[k] == rge[l]){
            plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
          } else {
            ex_lik_cor(rge[k], range_rge[[k]], rge[l], range_rge[[l]])
          }
        }
      }
    }
    dev.off()
    
  }
}


### Likelihood search
fn <- function(x){
  x0 <- val
  x0[4,] <- x[1:2]
  x0[5,] <- x[3:4]
  pars0 <- x[-(1:4)]
  cat("pars =", pars0, "\n")
  out <- lik_ABM(tab, x0, pars0)
  cat("lik =", out, "\n")
  if(is.infinite(out)) return(1e6)
  else return(-out)
}

sim_est <- matrix(0, nrow = 10*3^4, ncol = 18)
colnames(sim_est) <- paste(rep(c("anc_mdp", "node_mdp", "anc_rge", "node_rge", "sig_mdp", "sig_rge", "nu", "omega", "loglik"), times = 2), rep(c("_sim", "_est"), each = 9), sep = "")
i <- 1
for(n in c(0.1, 0.5, 0.9)){
  for(o in c(0.1, 0.5, 0.9)){
    for(m in c(0.001, 0.1, 1)){
      for(r in c(0.001, 0.1, 1)){
        for(it in 1:10){
          pars <- c(m, r, n, o)
          tree <- pbtree(n = 3)
          tab <- build_table(tree)
          sig2 <- pars[1:2]
          sig2_mdp <- pars[1]
          sig2_rge  <- pars[2]
          nu <- pars[3]
          omega <- pars[4]
          sim1 <- sim_ABM(tree, nu = nu, omega = omega, pars_rge = c(sig2_rge, 2, 0, Inf), pars_mdp = c(sig2_mdp, 10, -Inf, Inf), nsim = 1)
          sim <- sim1[[1]][1:3,]
          rownames(sim) <- tree$tip.label
          val <- sim1[[1]]
          opt <- lbfgs(c(mean(sim[,1]), mean(sim[,2]), mean(sim[,1]), mean(sim[,2]), 0.1, 0.1, 0.5, 0.5), fn, lower = rep(0, 8), upper = c(12, 8, 12, 8, 1,1,1,1), control= list(maxeval = 10000))
          sim_est[i,] <- c(val[4,],val[5,], pars, -fn(c(val[4,],val[5,], pars)), opt$par, -opt$value)
          i <- i + 1
        }
      }
    }
  }
}


for(v in c("anc_mdp", "node_mdp", "anc_rge", "node_rge", "sig_mdp", "sig_rge", "nu", "omega", "loglik")){
  jpeg(sprintf("Figures/sim_lik_%s.jpg", v), height = 20, width = 20, unit = "cm", res = 500)
  par(mar = c(5,5,3,1), las = 1, mfrow = c(3,3), oma = c(4,4,4,0))
  cols = c("#8cb369", "#f4a259", "#e71d36")
  plots <- T
  for(sig_mdp in unique(sim_est[,"sig_mdp_sim"])){
    plotn <- T
    for(sig_rge in unique(sim_est[,"sig_mdp_sim"])){
      plot(0, type = "n", bty = "n", xlim = c(0,12), ylim = range(sim_est[,sprintf("%s_est", v)]), xaxt = "n", ylab = v, xlab = expression(nu))
      mtext(side = 1, at = c(2,6,10), text = unique(sim_est[,"nu_sim"]), line = 1, cex = .6)
      if(plotn) mtext(bquote({sigma^2}[mdp] == .(sig_mdp)), 2, line = 5, cex = .5)
      if(plots) mtext(bquote({sigma^2}[rge] == .(sig_rge)), 3, line = 3, cex = .5)
      x1 <- 1
      plotn <- F
      for(nu in unique(sim_est[,"nu_sim"])){
        x2 <- 0
        for(om in unique(sim_est[,"omega_sim"])){
          cond <- sim_est[,"sig_mdp_sim"] == sig_mdp & sim_est[,"sig_rge_sim"] == sig_rge & sim_est[,"nu_sim"] == nu & sim_est[,"omega_sim"] == om
          cat("trait", mean(sim_est[cond,sprintf("%s_sim", v)]), "error", quantile(sim_est[cond,sprintf("%s_est", v)], c(.1, .9)), "\n")
          lines(x1+x2+c(-0.1,0.1), y = rep(mean(sim_est[cond,sprintf("%s_sim", v)]),2), lty = 2, col = "#626066")
          lines(x = rep(x1+x2, 2), y = quantile(sim_est[cond,sprintf("%s_est", v)], c(.1, .9)), col = "#0f7173")
          points(x1+x2, mean(sim_est[cond,sprintf("%s_est", v)]), col = cols[x2+1], pch = 16)
          x2 <- x2 + 1
        }
        x1 <- x1 + 4
      }
    }
    plots <- F
  }
  
  mtext(expression(omega), side = 1, cex = 1, line = 5, at = -22)
  legend(-20, -0.55, legend = unique(sim_est[,"omega_sim"]), col = cols, pch = 16, xpd = NA, bty = "n", horiz = T , cex = 1.2)
  
  dev.off()
}
