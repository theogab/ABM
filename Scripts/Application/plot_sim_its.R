## ABM simulations Output processing -> incomplete taxon sampling

lib <- c("bite", "ape", "coda","vioplot")
sapply(lib, library, character.only = T)

cols = c("#8cb369", "#f4a259", "#e71d36")
symb <- list(n = "N", nu = expression(nu), omega = expression(omega), rootlogv = expression(theta[eta]), rootm = expression(theta[mu]), siglogv = expression({sigma^2}[eta]), sigm = expression({sigma^2}[mu]))
res_BF_h <- readRDS("Results/res_BF_its_sim.rds")
res_BF <- readRDS("Results/res_BF.rds")
res_BF <- rbind(res_BF_h, res_BF[res_BF[,"n"] == 20 & res_BF[,"sigm"] == 0.1 & res_BF[,"siglogv"] == 0.1,-(2:3)])
sub_ACI <- res_BF[,"nu"] == 0.5 & res_BF[,"omega"] == 0
sub_ADI <- res_BF[,"nu"] == 0.5 & res_BF[,"omega"] == 0.5
sub_SDI <- res_BF[,"nu"] == 0 & res_BF[,"omega"] == 0.5
sub_INT <- res_BF[,"nu"] == 0.2 & res_BF[,"omega"] == 0.2
list_res <- list('Sym/Dis' = res_BF[sub_SDI,], 'Asym/Cons' = res_BF[sub_ACI,], 'Asym/Dis(-)' = res_BF[sub_INT, ], 'Asym/Dis(+)' = res_BF[sub_ADI,])
list_res_ess <- lapply(list_res, function(x){
  x[x[,"siglogv_ESS"] > 200 & x[,"sigm_ESS"] > 200,]
})

## % of rejected BM / nu = 0 / omega = 0
BF_nu <- lapply(list_res, function(x){ 
  sigs <- lapply(c(20, 25, 30, 40), function(prop_nu){
    x <- x[x[,"n"] == prop_nu,] 
    post <- (x[,"In"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05)/(0.95)))
  })
  names(sigs) <- c(0, 0.2, 0.33, 0.5)
  return(sigs)
})

BF_omega <- lapply(list_res, function(x){ 
  sigs <- lapply(c(20, 25, 30, 40), function(prop_nu){
    x <- x[x[,"n"] == prop_nu,] 
    post <- (x[,"Io"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05)/(0.95)))
  })
  names(sigs) <- c(0, 0.2, 0.33, 0.5)
  return(sigs)
})

prop_nu <- lapply(BF_nu, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
prop_omega <- lapply(BF_omega, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
colors <- c("#61988E", "#0A2472", "#E84855", "#AFB3F7")
cex.ax = 0.6
cex.la = 0.8


pdf("Figures/BF_per_scenario_its.pdf", height = 5, width = 4)
par(bty = "n", las = 1, mgp = c(2,0.8,0),cex.axis = cex.ax, cex.lab = cex.la, mar = c(0,0,0,0))
layout(rbind(c(1,1, 2, 2, 2,2,2), c(1,1,3,3,3,3,3)))

plot(0:1, bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n", type = "n")
leg.col <- character(9)
i=1
text(1.1,.994, "a. simulated values", font = 2, cex = 0.9, pos = 4, xpd = T)
text(rep(1.5, 6), seq(0.95, 0.35, len = 6), c("Asymmetric and Conserved\n (Asym/Cons)", "Asymmetric and Displaced\n (Asym/Dis(+))", "Asymmetric and Displaced\n (Asym/Dis(-))", "Symmetric and Conserved\n (Sym/Cons)", "Symmetric and Conserved, OU\n(OU)", "Symmetric and Displaced\n (Sym/Dis)"),
     cex = 0.5, pos = 1, xpd = T)
text(rep(1.8, 6), seq(0.9, 0.3, len = 6), c(expression(paste(nu, " = 0.5")),expression(paste(nu, " = 0.5")),expression(paste(nu, " = 0.2")),expression(paste(nu, " = 0")),expression(paste(nu, " = 0")),expression(paste(nu, " = 0"))),
     cex = 0.5, pos = 1, xpd = T)
text(rep(1.8, 6), seq(0.88, 0.28, len = 6), c(expression(paste(omega, " = 0")),expression(paste(omega, " = 0.5")),expression(paste(omega, " = 0.2")),expression(paste(omega, " = 0")),expression(paste(omega, " = 0")),expression(paste(omega, " = 0.5"))),
     cex = 0.5, pos = 1, xpd = T)
legend(1.1, 0.2, bty = "n", legend = c(1, 0.80, 0.66, 0.50), fill = colors, ncol = 2, x.intersp = 0.5, y.intersp = 1.2, title = "Taxon Sampling")

par(mar = c(4,3,2,2))
numod <- c(3,3,3,2)
boxplot(rep(-20,4)~as.factor(c("Sym/Dis", "Asym/Cons", "Asym/Dis(-)", "Asym/Dis(+)")), ylim = c(-10,25), xlab = "", ylab = expression(paste(2 * "log" * ("BF"[nu]))), main =  "", xaxt = "n")
mtext("b. asymmetry", font = 2, cex = 0.6, at = .5)
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
mtext("c. displacement", font = 2, cex = 0.6, at = .5)
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

dev.off()

#### Parameter estimation
est <- lapply(list_res, function(x){
  sigs <- lapply(c(20, 25, 30, 40), function(sig){
    x <- x[x[,"n"] == sig, ] 
    x[, grepl("mean", colnames(x))]
  })
  names(sigs) <- c(0, 0.2, 0.33, 0.5)
  return(sigs)
})

exp <- t(sapply(list_res, function(x){
  unique(x[, c("nu", "omega")])
}))

png("Figures/est_par_nu_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.33,0.5), ylim = c(0,1), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(nu)), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"nu_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"nu_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(exp[i,1],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()



png("Figures/est_par_omega_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.33,0.5), ylim = c(0,1), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(omega)), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"omega_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"omega_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(exp[i,2],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_sigm_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.33,0.5), ylim = c(0,1.5), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(sigma[mu]^2)), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"sigm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"sigm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_siglogv_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.33,0.5), ylim = c(0,1.5), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(sigma[zeta]^2)), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"siglogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"siglogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_rootm_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.33,0.5), ylim = c(8,12), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(theta[mu])), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"rootm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"rootm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(10,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_rootlogv_its.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-3,4)~c(0,0.2,0.33,0.5), ylim = c(-2,2), xlab = "Incomplete Taxon Sampling", ylab = expression(bar(theta[zeta])), main =  paste0("Simulated model : ", names(prop_nu)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[k]][,"rootlogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[2])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[k]][,"rootlogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

