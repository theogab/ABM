## Alpha skew distribution
dasnorm <- function(x, alpha = 0, mean = 0, sd = 1, log = FALSE){
  X <- (x-mean)/sd
  phi <- dnorm(X, log = F)
  (1/sd)^2*((1+alpha*X^2)/(1+alpha))*phi
}

x <- seq(-5,5,length.out = 10000)
plot(dasnorm(x, alpha = 1)~x, type = "l", xlab = "Desnsity")
lines(x, dasnorm(x, alpha = 1, mean = 1))

plot(phy, show.tip.label = F)
md <- get('last_plot.phylo', envir = .PlotPhyloEnv)
points(md$xx, md$yy, col = colvec[id], pch = 16)
points(seq(20, 80, length.out = 101), rep(-1,101), pch = 15, col = colvec, xpd = T)
axis(1, at = seq(19.5, 80.5, length.out = 6), labels = seq(0, 1, length.out = 6), cex.axis = 0.5, line = 1.3, mgp = c(1,0.5, 0))

tra <- sim_tree[[1]]
cl <- extract.clade(phy, 31)
X <- sim_tree[[1]][c(31,9,10),]
Xis(0.5,0.5,X[1])
Xis(0.5,0.5,X[1],T)
dnorm(X[3],Xis(0.5,0.5,X[1],T),phy$edge.length[20]*0.5)


## Simulations truncated normal
sample.size <- c(1e2,1e3,1e5)
mn <- 10
std <- 1
split <- 9
init.sample <- rnorm(sample.size[3], mn, std)
trunc.sample <- init.sample[init.sample < split]
jpeg("Figures/sim_sampling.jpg", height = 50, width = 8, units = "cm", res = 600)
layout(1:11)
par(mar = c(2,2,2,1))
hist(trunc.sample, col = "#2e86ab", main = "Speciation Time")
lines(rep(mean(c(split,min(init.sample))),2),c(0,1e4), lty = 2, col = "red")
for(i in 1:1e4){
  trunc.sample <- trunc.sample + rnorm(length(trunc.sample),0,0.05)
  if(i %in% seq(1e3, 1e4, length.out = 10)){
    hist(trunc.sample, col = "#2e86ab", main = sprintf("ngen = %s",i))
    lines(rep(mean(c(split,min(init.sample))),2),c(0,1e4), lty = 2, col = "red")
  }
}
dev.off()

jpeg("Figures/Range_sampling.jpg", height = 16, width = 16, units = "cm", res = 600)
plot_hp(hpfun("Normal", c(mn,std)), xlab = "Max T°C")
lines(rep(split,2), c(0,0.4), lty = 2, col = "red")
points(init.sample, rep(-0.01, sample.size[1]), pch = 3)
lines(rep(min(init.sample),2), c(0,0.1), col = "red")
arrows(min(init.sample),0.1,split,0.1,code = 3, col = "red", length= .1)
text(mean(c(split,min(init.sample))),0.13,labels = expression(X[min](T[s])), col = "red")
lines(rep(max(init.sample),2), c(0,0.1), col = "red")
arrows(max(init.sample),0.1,split,0.1,code = 3, col = "red", length= .1)
text(mean(c(split,max(init.sample))),0.13,labels = expression(X[max](T[s])), col = "red")
dev.off()

diff(range(init.sample))
diff(c(split,min(init.sample)))
