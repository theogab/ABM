## Plot simulations ABM

lib <- c("ape")
sapply(lib, library, character.only = T)

main.plot <- function(lik){
  at.gr <- rep(seq(5, 15, length.out = 3), each = 5) + rep(seq(1,5,length.out = 5), 3)
  lik[lik %in% c(-Inf, Inf)] <- NA
  boxplot(t(lik[,-c(1,2)]), at = at.gr, col = rev(rep(c("#d53e4f", "#66c2a5","#3288bd"), each = 5)),
          outcol = rev(rep(c("#d53e4f", "#66c2a5","#3288bd"), each = 5)), names = rep(seq(0.1,0.5, length.out = 5), 3),
          xlab = "rho", ylab = "log(likelihood)", cex.axis = .8, las = 1, ylim =range(lik[,-c(1,2)], na.rm = T))
}

for(n in 1:3){
  jpeg(sprintf("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Figures/lik_abm_tree-%s.jpg", c(10,50,200)[n]), height = 20, width = 40, units = "cm", res = 600)
  cat("tree-", c(10,50,200)[n], "\n")
  par(mfrow = c(2,4), oma = c(0, 0, 8, 0))
  for (rho in c(0.05,0.5)){
    for(h in c(0.1, 0.9)){
      for(sig2 in c(0.01, 0.5)){
        if(sig2 == 0.5) par(mar = c(5, 5, 1, 6))
        else par(mar = c(5, 10, 1, 1))
        cat("rho = ", rho, "sig2 = ", sig2, "h = ", h, "\n")
        load(sprintf("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Results/lik_rho-%s_h-%s_sig-%s.Rdata", rho, h, sig2))
        main.plot(lik[[n]][[2]])
        if(rho == 0.05){
          mtext(bquote(sigma^2 == .(sig2)), side = 3, line = 1)
          if(sig2 == 0.01) mtext(bquote(h == .(h)), side = 3, line = 5, at = 22)
        }
        if(sig2 == 0.01) mtext(bquote(rho == .(rho)), side = 2, line = 5, las = 1)
      }
    }
  }
  dev.off()
}


## plot trees
set.seed(300)
tree1 <- pbtree(n = 10, scale = 100)
tree2 <- pbtree(n = 50, scale = 100)
tree3 <- pbtree(n = 200, scale = 100)
tree <- list(tree1, tree2, tree3)

for(n in 1:3){
  jpeg(sprintf("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Figures/tree_abm_sim-%s.jpg", c(10,50,200)[n]), height = 40, width = 20, units = "cm", res = 600)
  cat("tree-", c(10,50,200)[n], "\n")
  plot(tree[[n]], show.tip.label = F, x.lim = c(0, 180))
  tips.coords <- get('last_plot.phylo', envir = .PlotPhyloEnv)$yy[1:length(tree[[n]]$tip.label)]
  i <- 0
  for(sig2 in c(0.01, 0.5)){
    for (rho in c(0.05,0.5)){
      for(h in c(0.1, 0.9)){
        cat("rho = ", rho, "sig2 = ", sig2, "h = ", h, "\n")
        load(sprintf("/Users/tgaboria/Dropbox/Documents/JIVE/ABM/Results/lik_rho-%s_h-%s_sig-%s.Rdata", rho, h, sig2))
        val <- lik[[n]][[1]][1:length(tree[[n]]$tip.label),1]
        val <- val/max(val)*10
        points(rep(100+i, length(tips.coords)), tips.coords, cex = val, pch = 16, col = ifelse(rho == 0.05, "#d53e4f", "#3288bd"))
        
        i <- i + 10
      }
    }
  }
  dev.off()
}

