## ABM simulations

# asymetric inheritance of x during speciation
# poisson process for speciation
sim_ABM <- function(theta0 = 0, sigma2 = 1, rho = 1, delta = 1, lambda = 5, nspmax = 100, ngen = 100){
  x <- matrix(NA, nrow = nspmax, ncol = ngen+1)
  x[1,1] <- theta0 
  t.sp <- rep(NA,nspmax)
  t.sp[1] <- rpois(1, lambda)
  i.sp <- rep(NA, nspmax)
  i.sp[1] <- 0
  spe <- rep(F, nspmax)
  for (i in 2:(ngen+1)){
    nsp <- sum(!is.na(i.sp))
    # test for speciation
    spe[!is.na(i.sp)] <- i.sp[!is.na(i.sp)] == t.sp[!is.na(i.sp)]
    i.sp <- i.sp + 1
    if (any(spe) & nsp + sum(spe) <= nspmax){
      # asymetric inheritance of trait
      x[,i] <- x[,i-1]
      x[spe,i] <- rho * x[spe, i-1]
      x[nsp+1:sum(spe),i] <- (rho + delta + 1) * x[spe, i-1]
      # trait evolution under brownian motion
      x[,i] <- rnorm(nrow(x), x[,i], sigma2)
      # fill vector for speciation process
      i.sp[spe] <- 0
      t.sp[spe] <- rpois(sum(spe),lambda)
      i.sp[nsp+1:sum(spe)] <- rep(0, sum(spe))
      t.sp[nsp+1:sum(spe)] <- rpois(sum(spe),lambda)
    } else {
      # trait evolution under brownian motion
      x[,i] <- rnorm(nrow(x), x[,i-1], sigma2)
    }
  }
  return(x)
}

plot_ABM_sim <- function(sim, col, new = T, xlim, ylim, reverse = F){
  if(new){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    axis(ifelse(reverse, 4,2))
    mtext("trait variance", side = ifelse(reverse, 4,2), line = 2)
    axis(1)
    mtext("generations", side = 1, line = 2)
  } 
  if(reverse){
    apply(sim, 1, function(x){
      lines(rev(x), col = col)
    })
  } else {
    apply(sim, 1, function(x){
      lines(x, col = col)
    })
  }
}


pars <- cbind(delta = c(0.1, 0.9), rho = c(0.1, 0.9), sigma = c(0.01, 0.5))
comb <- rbind(c(1,1,1), c(1,2,1), c(1,2,2), c(1,1,2),c(2,1,1),c(2,2,1), c(2,2,2),c(2,1,2))
sim <- lapply(1:nrow(comb), function(i){
  x <- comb[i,]
  cat(x, "\n")
  set.seed(300)
  sim_ABM(lambda = 4000, delta = pars[x[1],1], rho = pars[x[2],2], ngen = 1e4, nspmax = 3, sigma = pars[x[3],3], theta0 = 2)
})


jpeg("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Figures/rho_effect.jpg", width = 30, height = 20, units = "cm", res = 600)
layout(cbind(c(1,2), c(3,4)))
# delta = 0.1; sigma = 0.01
par(mar = c(2,6,2,-1)+1)
plot_ABM_sim(sim[[1]], "#518eb6", T, xlim = c(0, ncol(sim[[1]])), ylim = c(0, max(sim[[1]], na.rm = T)))
plot_ABM_sim(sim[[2]], "#ac5e63", F, xlim = c(0, ncol(sim[[1]])), ylim = c(0, max(sim[[1]], na.rm = T)))
mtext(expression(sigma^2 == 0.01), 3, line = 1)
mtext(expression(delta == 0.1), 2, line = 3.5, las = 1)
legend(100, 1, legend = c(expression(rho == 0.1), expression(rho == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# delta = 0.9; sigma = 0.01
par(mar = c(3,6,1,-1)+1)
plot_ABM_sim(sim[[5]], "#518eb6", T, xlim = c(0, ncol(sim[[5]])), ylim = c(0, max(sim[[5]], na.rm = T)))
plot_ABM_sim(sim[[6]], "#ac5e63", F, xlim = c(0, ncol(sim[[5]])), ylim = c(0, max(sim[[5]], na.rm = T)))
mtext(expression(delta == 0.9), 2, line = 3.5, las = 1)
legend(100, 1, legend = c(expression(rho == 0.1), expression(rho == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# delta = 0.1; sigma = 0.5
par(mar = c(2,3,2,2)+1)
plot_ABM_sim(sim[[3]], "#518eb6", T, xlim = c(0, ncol(sim[[3]])), ylim = c(0, max(sim[[3]], na.rm = T)))
plot_ABM_sim(sim[[4]], "#ac5e63", F, xlim = c(0, ncol(sim[[3]])), ylim = c(0, max(sim[[3]], na.rm = T)))
mtext(expression(sigma^2 == 0.5), 3, line = 1)
legend(100, 80, legend = c(expression(rho == 0.1), expression(rho == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# delta = 0.9; sigma = 0.5
par(mar = c(3,3,1,2)+1)
plot_ABM_sim(sim[[7]], "#518eb6", T, xlim = c(0, ncol(sim[[7]])), ylim = c(0, max(sim[[7]], na.rm = T)))
plot_ABM_sim(sim[[8]], "#ac5e63", F, xlim = c(0, ncol(sim[[7]])), ylim = c(0, max(sim[[7]], na.rm = T)))
legend(100, 90, legend = c(expression(rho == 0.1), expression(rho == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
dev.off()

jpeg("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Figures/delta_effect.jpg", width = 30, height = 20, units = "cm", res = 600)
layout(cbind(c(1,2), c(3,4)))
# rho = 0.05; sigma = 0.01
par(mar = c(2,6,2,-1)+1)
plot_ABM_sim(sim[[1]], "#518eb6", T, xlim = c(0, ncol(sim[[1]])), ylim = c(0, max(sim[[1]], na.rm = T)))
plot_ABM_sim(sim[[5]], "#ac5e63", F, xlim = c(0, ncol(sim[[1]])), ylim = c(0, max(sim[[1]], na.rm = T)))
mtext(expression(sigma^2 == 0.01), 3, line = 1)
mtext(expression(rho == 0.1), 2, line = 3.5, las = 1)
legend(100, 1, legend = c(expression(delta == 0.1), expression(delta == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# rho = 0.5; sigma = 0.01
par(mar = c(3,6,1,-1)+1)
plot_ABM_sim(sim[[2]], "#518eb6", T, xlim = c(0, ncol(sim[[2]])), ylim = c(0, max(sim[[2]], na.rm = T)))
plot_ABM_sim(sim[[6]], "#ac5e63", F, xlim = c(0, ncol(sim[[2]])), ylim = c(0, max(sim[[2]], na.rm = T)))
mtext(expression(rho == 0.9), 2, line = 3.5, las = 1)
legend(100, 1, legend = c(expression(delta == 0.1), expression(delta == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# rho = 0.05; sigma = 0.5
par(mar = c(2,3,2,2)+1)
plot_ABM_sim(sim[[4]], "#518eb6", T, xlim = c(0, ncol(sim[[8]])), ylim = c(0, max(sim[[8]], na.rm = T)))
plot_ABM_sim(sim[[8]], "#ac5e63", F, xlim = c(0, ncol(sim[[8]])), ylim = c(0, max(sim[[8]], na.rm = T)))
mtext(expression(sigma^2 == 0.5), 3, line = 1)
legend(100, 90, legend = c(expression(delta == 0.1), expression(delta == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
# rho = 0.05; sigma = 0.5
par(mar = c(3,3,1,2)+1)
plot_ABM_sim(sim[[3]], "#518eb6", T, xlim = c(0, ncol(sim[[7]])), ylim = c(0, max(sim[[7]], na.rm = T)))
plot_ABM_sim(sim[[7]], "#ac5e63", F, xlim = c(0, ncol(sim[[7]])), ylim = c(0, max(sim[[7]], na.rm = T)))
legend(100, 90, legend = c(expression(delta == 0.1), expression(delta == 0.9)), lty = 1, col = c("#518eb6", "#ac5e63"), bty = "n")
dev.off()



v1 <- function(rho, Va){
  return(Va*ifelse(rho<.5, rho, 2*rho-1))
}
v2 <- function(rho, Va){
  return(Va*ifelse(rho<.5, 1-rho, 1))
}

colvec <- colorRampPalette(rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4",
                              "#88ddaa", "#66c2a5", "#37b69b", "#3288bd", "#71acbc")))(101)
rhos <- seq(1e-3, 1, length = 100)
Va <- 10
out_v1 <- v1(rhos, Va)
names(out_v1) <- rhos
out_v2 <- v2(rhos, Va)
names(out_v2) <- rhos

layout(t(1:2))
plot(rep(rhos, times = 100), rep(deltas, each = 100), pch = 15, col = colvec[round(out_v1/max(c(out_v1, out_v2))*100)+1],
     main = expression(paste("Effect of ", rho, " and ", delta, " on ", v[1])), xlab = expression(rho), ylab = expression(delta)) 
plot(rep(rhos, times = 100), rep(deltas, each = 100), pch = 15, col = colvec[round(out_v2/max(c(out_v1, out_v2))*100)+1],
     main = expression(paste("Effect of ", rho, " and ", delta, " on ", v[2])), xlab = expression(rho), ylab = expression(delta))

hist(out_v1)
hist(out_v2)
