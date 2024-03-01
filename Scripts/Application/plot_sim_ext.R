## ABM simulations Output processing -> extgeneous processes

lib <- c("bite", "ape", "coda","vioplot")
sapply(lib, library, character.only = T)

cols = c("#8cb369", "#f4a259", "#e71d36")
symb <- list(n = "N", nu = expression(nu), omega = expression(omega), rootlogv = expression(theta[zeta]), rootm = expression(theta[mu]), siglogv = expression({sigma^2}[zeta]), sigm = expression({sigma^2}[mu]))
res_BF_h <- readRDS("Results/res_BF_ext_sim.rds")
res_BF <- readRDS("Results/res_BF.rds")
res_BF <- rbind(res_BF_h, cbind(rep(0, 500), res_BF[res_BF[,"n"] == 20 & res_BF[,"sigm"] == 0.1 & res_BF[,"siglogv"] == 0.1,-(1:3)]))
sub_ACI <- res_BF[,"nu"] == 0.5 & res_BF[,"omega"] == 0
sub_ADI <- res_BF[,"nu"] == 0.5 & res_BF[,"omega"] == 0.5
sub_SDI <- res_BF[,"nu"] == 0 & res_BF[,"omega"] == 0.5
sub_INT <- res_BF[,"nu"] == 0.2 & res_BF[,"omega"] == 0.2
list_res <- list('Sym/Dis' = res_BF[sub_SDI,], 'Asym/Cons' = res_BF[sub_ACI,], 'Asym/Dis(-)' = res_BF[sub_INT, ], 'Asym/Dis(+)' = res_BF[sub_ADI,])
list_res_ess <- lapply(list_res, function(x){
  x[x[,"siglogv_ESS"] > 200 & x[,"sigm_ESS"] > 200,]
})

## % of rejected BM / nu = 0 / omega = 0
BF <- lapply(list_res, function(x){ 
  sigs <- lapply(c(0, 0.2, 0.5, 0.8), function(prop){
    x <- x[x[,"mu"] == prop,] 
    post <- (x[,"In"])*(x[,"Io"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05^2)/(0.95^2)))
  })
  names(sigs) <- c(0, 0.2, 0.5, 0.8)
  return(sigs)
})



BF_nu <- lapply(list_res, function(x){ 
  sigs <- lapply(c(0, 0.2, 0.5, 0.8), function(prop){
    x <- x[x[,"mu"] == prop,] 
    post <- (x[,"In"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05)/(0.95)))
  })
  names(sigs) <- c(0, 0.2, 0.5, 0.8)
  return(sigs)
})

BF_omega <- lapply(list_res, function(x){ 
  sigs <- lapply(c(0, 0.2, 0.5, 0.8), function(prop){
    x <- x[x[,"mu"] == prop,] 
    post <- (x[,"Io"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05)/(0.95)))
  })
  names(sigs) <- c(0, 0.2, 0.5, 0.8)
  return(sigs)
})


prop <- lapply(BF, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
prop_nu <- lapply(BF_nu, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
prop_omega <- lapply(BF_omega, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
colors <- c("#61988E", "#0A2472", "#E84855", "#AFB3F7")
cex.ax = 0.5
cex.la = 0.8

png("Figures/prop_BF_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(prop[[i]], beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[nu], " = 0 &  ", I[omega], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0,0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext(expression(mu), 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()

png("Figures/prop_BF_nu_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(prop_nu[[i]], beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[nu], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0,0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext(expression(mu), 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()

png("Figures/prop_BF_omega_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(prop_omega[[i]], beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[omega], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0, 0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext(expression(mu), 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()

est <- lapply(list_res, function(x){
  sigs <- lapply(c(0, 0.2, 0.5, 0.8), function(sig){
    x <- x[x[,"mu"] == sig, ] 
    x[, grepl("mean", colnames(x))]
  })
  names(sigs) <- c(0.2, 0.5, 0.8, 1)
  return(sigs)
})

exp <- t(sapply(list_res, function(x){
  unique(x[, c("nu", "omega")])
}))



pdf("Figures/BF_per_scenario_ext.pdf", height = 5, width = 5)
par(bty = "n", las = 1, mgp = c(1.2,0.4,0),cex.axis = cex.ax, cex.lab = cex.la, mar = c(0,0,0,0), tck = -0.03)
layout(rbind(c(1, 1, 1,1,1,3,3,3),c(1, 1, 1,1,1,3,3,3),c(1, 1, 1,1,1,4,4,4), c(2,2,2,2,2,4,4,4), c(2,2,2,2,2,5,5,5), c(2,2,2,2,2,5,5,5)))

par(mar = c(4,3,2,2))
numod <- c(3,3,3,2)
boxplot(rep(-20,4)~as.factor(c("Sym/Dis", "Asym/Cons", "Asym/Dis(-)", "Asym/Dis(+)")), ylim = c(-10,25), xlab = "", ylab = expression(paste(2 * "log" * ("BF"[nu]))), main =  "", xaxt = "n")
mtext("a. asymmetry", font = 2, cex = 0.6, at = .5)
axis(1, 1:4, labels = rep("",4))
x <- 0.7
for(i in c(2,4,3,1)){
  for(k in 1:4){
    boxplot(unlist(do.call(cbind,BF_nu[[i]])[,k]), at = x, col = adjustcolor(colors[k], exp(1)/exp(numod[i])), add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(unlist(do.call(cbind,BF_nu[[i]])[,k])),2), col = "#FFFFFF")
    x <- x + 0.175
  }
  x <- x + 0.3
}
lines(c(0.3, 4.3), y = rep(6,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2, xpd = T)
text(c(4.3, 4.3), c(4.5,7.5), c( expression(nu == 0), expression(nu > 0)), cex = 0.6, pos=4, xpd = T)

par(mar = c(4,3,2,2))
numod <- c(2,3,3,3)
boxplot(rep(-20,4)~as.factor(c("Sym/Dis", "Asym/Cons", "Asym/Dis(-)", "Asym/Dis(+)")), ylim = c(-10,25), xlab = "", ylab = expression(paste(2 * "log" * ("BF"[omega]))), main =  "", xaxt = "n")
mtext("b. displacement", font = 2, cex = 0.6, at = .5)
axis(1, 1:4, labels = rep("", 4))
x <- 0.7
for(i in c(2,4,3,1)){
  for(k in 1:4){
    boxplot(unlist(do.call(cbind,BF_omega[[i]])[,k]), at = x, col = adjustcolor(colors[k], exp(1)/exp(numod[i])), add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(unlist(do.call(cbind,BF_omega[[i]])[,k])),2), col = "#FFFFFF")
    x <- x + 0.175
  }
  x <- x + 0.3
}
lines(c(0.3, 4.3), y = rep(6,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2, xpd = T)
text(c(4.3, 4.3), c(4.5,7.5), c( expression(omega == 0), expression(omega > 0)), cex = 0.6, pos=4, xpd = T)

par(mar= c(1,2.5,2.5,1.5))
boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = expression(mu), ylab = expression(bar(nu)), main =  "")
mtext("c.", 3, at = 0, font = 2, cex = 0.6)
for(k in 1:4){
  boxplot(est[[4]][[k]][,"nu_mean"], at = k, col = adjustcolor(colors[k], exp(1)/exp(numod[i])), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
  lines(k + c(-0.1, 0.1), y = rep(median(est[[4]][[k]][,"nu_mean"]),2), col = "#FFFFFF")
  lines(k - c(-0.4, 0.4), y = rep(exp[4,1],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1)
}
par(mar= c(2,2.5,1.5,1.5))
boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = expression(mu), ylab = expression(bar(omega)), main =  "")
mtext("d.", 3, at = 0, font = 2, cex = 0.6)
for(k in 1:4){
  boxplot(est[[4]][[k]][,"omega_mean"], at = k, col = adjustcolor(colors[k], exp(1)/exp(numod[i])), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
  lines(k + c(-0.1, 0.1), y = rep(median(est[[4]][[k]][,"omega_mean"]),2), col = "#FFFFFF")
  lines(k - c(-0.4, 0.4), y = rep(exp[4,2],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1)
}

par(mar= c(2,2.5,1.5,1.5))
boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,2), xlab = expression(mu), ylab = expression(bar(sigma[zeta]^2)), main =  "")
mtext("e.", 3,at = 0, font = 2, cex = 0.6)
for(k in 1:4){
  boxplot(est[[4]][[k]][,"siglogv_mean"], at = k, col = adjustcolor(colors[k], exp(1)/exp(numod[i])), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
  lines(k + c(-0.1, 0.1), y = rep(median(est[[4]][[k]][,"siglogv_mean"]),2), col = "#FFFFFF")
  lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1)
}

dev.off()

#### Parameter estimation

png("Figures/est_par_nu_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = expression(mu), ylab = expression(bar(nu)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"nu_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"nu_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(exp[i,1],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()



png("Figures/est_par_omega_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = expression(mu), ylab = expression(bar(omega)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"omega_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"omega_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(exp[i,2],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_sigm_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1.5), xlab = expression(mu), ylab = expression(bar(sigma[mu]^2)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"sigm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"sigm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_siglogv_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1.5), xlab = expression(mu), ylab = expression(bar(sigma[zeta]^2)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"siglogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"siglogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_rootm_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(8,12), xlab = expression(mu), ylab = expression(bar(theta[mu])), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"rootm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"rootm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(10,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_rootlogv_ext.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(3,4)~c(0,0.2,0.5,0.8), ylim = c(-15,2), xlab = expression(mu), ylab = expression(bar(theta[zeta])), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"rootlogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"rootlogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

