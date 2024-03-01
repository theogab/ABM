## ABM simulations Output processing -> Heterogeneous processes

lib <- c("bite", "ape", "coda","vioplot")
sapply(lib, library, character.only = T)

cols = c("#8cb369", "#f4a259", "#e71d36")
symb <- list(n = "N", nu = expression(nu), omega = expression(omega), rootlogv = expression(theta[eta]), rootm = expression(theta[mu]), siglogv = expression({sigma^2}[eta]), sigm = expression({sigma^2}[mu]))
res_BF_h <- readRDS("Results/res_BF_hetero_sim.rds")
res_BF <- readRDS("Results/res_BF.rds")
res_BF <- rbind(res_BF_h, cbind(rep(1, 500), res_BF[res_BF[,"n"] == 20 & res_BF[,"sigm"] == 0.1 & res_BF[,"siglogv"] == 0.1,-(1:3)]))
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
  sigs <- lapply(c(0.2, 0.5, 0.8,1), function(prop){
    x <- x[x[,"prop"] == prop,] 
    post <- (x[,"In"])*(x[,"Io"])
    post[post == 1] <- 0.999
    post[post == 0] <- 0.001
    2 * log((post/(1 - post))/((0.05^2)/(0.95^2)))
  })
  names(sigs) <- c(0.2, 0.5, 0.8, 1)
  return(sigs)
})



BF_nu <- lapply(list_res, function(x){ 
  sigs <- lapply(c(0.2, 0.5, 0.8, 1), function(prop){
    x <- x[x[,"prop"] == prop,] 
      post <- (x[,"In"])
      post[post == 1] <- 0.999
      post[post == 0] <- 0.001
      2 * log((post/(1 - post))/((0.05)/(0.95)))
    })
    names(sigs) <- c(0.2, 0.5, 0.8, 1)
    return(sigs)
  })

BF_omega <- lapply(list_res, function(x){ 
  sigs <- lapply(c(0.2, 0.5, 0.8,1), function(prop){
    x <- x[x[,"prop"] == prop,] 
      post <- (x[,"Io"])
      post[post == 1] <- 0.999
      post[post == 0] <- 0.001
      2 * log((post/(1 - post))/((0.05)/(0.95)))
    })
    names(sigs) <- c(0.2, 0.5, 0.8,1)
    return(sigs)
  })


prop <- lapply(BF, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
prop_nu <- lapply(BF_nu, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
prop_omega <- lapply(BF_omega, function(x) sapply(x, function(y) sum(y > 6)/length(y)))
colors <- c("#61988E", "#0A2472", "#E84855", "#AFB3F7")
cex.ax = 0.6
cex.la = 0.8

png("Figures/prop_BF_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(rev(prop[[i]]), beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[nu], " = 0 &  ", I[omega], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0,0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext("Probability of SCI at nodes", 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()

png("Figures/prop_BF_nu_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(rev(prop_nu[[i]]), beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[nu], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0,0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext("Probability of SCI at nodes", 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()

png("Figures/prop_BF_omega_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(mar = c(4,5,2,1), las = 1, mgp = c(3,1,0.5), cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  barplot(rev(prop_omega[[i]]), beside = T, col = sapply(seq(0.2,0.8, len = 4), function(tr) adjustcolor(colors[1],tr)), ylab = expression(paste(I[omega], " = 0 (BF > 6)")), main = paste0("Simulated model : ", names(prop)[i]), ylim = c(0,1), xlab = "", xaxt = "n")
  mtext(c(0, 0.2,0.5,0.8), 1, line = 0.5, cex = cex.ax/1.5, at = c(.7,1.9,3.1,4.3))
  mtext("Probability of SCI at nodes", 1, line = 1.5, cex = cex.la/1.5)
}
dev.off()


pdf("Figures/BF_per_scenario_hetero.pdf", height = 5, width = 4)
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
legend(1.1, 0.2, bty = "n", legend = c(0, 0.2, 0.5, 0.8), fill = colors, ncol = 2, x.intersp = 0.5, y.intersp = 1.2, title = "P(Sym/Cons)")

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
    sigs <- lapply(c(0.2, 0.5, 0.8, 1), function(sig){
      x <- x[x[,"prop"] == sig, ] 
      x[, grepl("mean", colnames(x))]
    })
    names(sigs) <- c(0.2, 0.5, 0.8, 1)
    return(sigs)
})

exp <- t(sapply(list_res, function(x){
  unique(x[, c("nu", "omega")])
}))

png("Figures/est_par_nu_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(nu)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
      boxplot(est[[i]][[5-k]][,"nu_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
      lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"nu_mean"]),2), col = "#FFFFFF")
      lines(k - c(-0.4, 0.4), y = rep(exp[i,1],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
    }
}
dev.off()



png("Figures/est_par_omega_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(omega)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[5-k]][,"omega_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"omega_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(exp[i,2],2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_sigm_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1.5), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(sigma[mu]^2)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[5-k]][,"sigm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"sigm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_siglogv_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(0,1.5), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(sigma[zeta]^2)), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[5-k]][,"siglogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"siglogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0.1,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


png("Figures/est_par_rootm_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-1,4)~c(0,0.2,0.5,0.8), ylim = c(8,12), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(theta[mu])), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[5-k]][,"rootm_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"rootm_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(10,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()

png("Figures/est_par_rootlogv_hetero.png", height = 15, width = 20, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2.6,0.8,0),cex.axis = cex.ax, cex.lab = cex.la)
layout(rbind(1:2, 3:4))
for(i in 1:4){
  boxplot(rep(-3,4)~c(0,0.2,0.5,0.8), ylim = c(-2,2), xlab = "Probability of Sym/Cons at nodes", ylab = expression(bar(theta[zeta])), main =  paste0("Simulated model : ", names(prop)[i]))
  for(k in 1:4){
    boxplot(est[[i]][[5-k]][,"rootlogv_mean"], at = k, col = adjustcolor(colors[2], exp(1)/exp(5-k)), add = T, boxwex = .8, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(k + c(-0.2, 0.2), y = rep(median(est[[i]][[5-k]][,"rootlogv_mean"]),2), col = "#FFFFFF")
    lines(k - c(-0.4, 0.4), y = rep(0,2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  }
}
dev.off()


### Convergence ####
png("Figures/convergence.png", height = 13, width = 18, unit = "cm", res = 500)
par(bty = "n", las = 1, mgp = c(2,0.8,0),cex.axis = cex.ax, cex.lab = cex.la, mar = c(4,4,1,1))
layout(rbind(1:3, 1:3, 1:3, 1:3, 1:3, 4:6, 4:6, 4:6, 4:6, 4:6, 7))
sn <- c(20, 50, 100)
ss <- c(0.1, 0.5, 1)
#### OMEGA #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[omega])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "omega_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "omega_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}

#### Root mean #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[theta[mu]])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "rootm_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "rootm_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}
#### Evolutionary rate mean #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[sigma[mu]^2])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "sigm_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "sigm_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}

#### NU #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[nu])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "nu_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "nu_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}

#### Root log sd #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[theta[log(sigma)]])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "rootlogv_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "rootlogv_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}

#### Evolutionary rate log sd #####
boxplot(rep(-20,3)~as.factor(c(20, 50, 100)), ylim = c(0,10), xlab = "Number of tips", ylab = expression(log(ESS[sigma[log(sigma)]^2])))
x <- 0.8
for(i in 1:3){
  for(k in 1:3){
    boxplot(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "siglogv_ESS"]), at = x, col = colors[k], add = T, boxwex = .2, xaxt = "n", yaxt = "n",outline = F ,outwex = 0, border = colors[k])
    lines(x + c(-0.02, 0.02), y = rep(median(log(res_BF[res_BF[,"n"] == sn[i] & res_BF[,"sigm"] == ss[k], "siglogv_ESS"])),2), col = "#FFFFFF")
    x <- x + 0.2
  }
  lines(x + c(-0.7, 0.1), y = rep(log(200),2), col = adjustcolor("#000000", .5), lty = 2, lwd = 1.2)
  x <- x + 0.4
}

par(mar = c(0,0,0,0))
plot(0:1, bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n", type = "n")
legend(1.4, 1, bty = "n", legend = c(expression(paste(sigma^2, "= 0.1")) , expression(paste(sigma^2, "= 0.5")), expression(paste(sigma^2, "= 1"))), fill = colors, ncol = 3)
dev.off()


