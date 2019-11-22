##Figures ABM search


## First example
jpeg("Figures/bifurcation.jpg", height = 8, width = 8, units = "cm", res = 600)
par(mar = c(3,1,0,1))
plot.new()
plot.window(xlim = c(1, 9), ylim = c(1, 9))
segments(c(2,5.1,5.1), c(5,5.2,4.8), c(4.9, 8, 8), c(5, 8.1, 1.9))
points(c(4.9,5.1,5.1),c(5,5.2,4.8),pch = 16,cex = .6)
text(c(4.5,5.1,5.1),c(4.8,5.2,4.8),labels = c(expression(X[0[s]]), expression(X[1[s]]), expression(X[2[s]])),pos = c(3,3,1), cex = .8)
axis(1, at = c(2,5,8), labels = c(expression(T[r]), expression(T[s]), 0))
dev.off()



##Parametrization
v1 <- function(omega, nu, X0){
  return(X0/2*(omega*(1-nu^2)/(1+nu)+1+nu))
}

v2 <- function(omega, nu, X0){
  return(X0/2*(omega*(1-nu^2)/(1+nu)+1-nu))
}
colvec <- colorRampPalette(rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4",
                                 "#88ddaa", "#66c2a5", "#37b69b", "#3288bd", "#71acbc")))(101)

nu <- seq(0,1, length.out = 1e4)
omega <- seq(0,1, length.out = 20)

jpeg("Figures/descendants.jpg", height = 7, width = 8, units = "cm", res = 600)
par(mar = c(2,2,1,5), las = 1, mgp = c(0.8,0.3,0), tck = -.01)
plot(v1(1,nu,1)~nu, type = "n", ylim = c(0,1), ylab = expression(X[i[s]]/X[0[s]]), xlab = expression(nu),
     cex.lab = .6, cex.axis = .5, xaxt = "n")
axis(1, cex.axis = .5, mgp = c(0.8,0,0))
for(i in omega){
  lines(nu, v1(i,nu,1), col = adjustcolor("#d95d39",i))
  lines(nu, v2(i,nu,1), col = adjustcolor("#3943b7",i))
  lines(c(1.1,1.2), rep(i,2), xpd = T, col = adjustcolor("#d95d39",i))
  lines(c(1.3,1.4), rep(i,2), xpd = T, col = adjustcolor("#3943b7",i))
  text(1.5,i, labels = round(i,2), cex = 0.5, xpd = T)
}
 text(c(1.15, 1.35, 1.5), rep(1.1, 3), xpd = T, labels = c(expression(X[1[s]]),expression(X[2[s]]), expression(omega)), cex = .5)
dev.off()

nu <- seq(0,1, length.out = 6)
omega <- seq(0,1, length.out = 6)

jpeg("Figures/Parameter_space.jpg", height = 7, width = 8, units = "cm", res = 600)
par(oma = c(2,2,1,1), mfrow = c(6,6), mar =c(1.5,1,0,.5))
for(o in rev(omega)){
  for(n in nu){
    barplot(c(v1(o,n,1),v2(o,n,1)), col = c("#d95d39", "#3943b7"), axes = F, ylim = c(0,1), space = .8,
            border = NA)
    axis(1, cex.axis = .5, at = c(1.5,3.4), labels = c(expression(X[1[s]]), expression(X[2[s]])),
         mgp = c(1,0,0), tck = 0, lwd = 0.5)
    if(n == nu[1]) lines(c(0,0), c(-2,2), xpd = T)
    if(o == omega[1]){
      lines(c(-1,6), c(-1,-1), xpd = T)
      mtext(n, outer = T, side = 1, line = 0, at = n*.85+.08, cex = .3, las = 1)
    } 
  }
  mtext(o, outer = T, side = 2, line = 0, at = o*.8+.12, cex = .3, las = 1)
}
  mtext(expression(omega), side = 2, line = 0.8, cex = .5, outer = T)
  mtext(expression(nu), side = 1, line = 0.8, cex = .5, outer = T)
dev.off()
