### LASEC PLOT FUNCTION ###

plot.lasec <- function(lasec.result) {
  iter <- dim(lasec.result$fit)[1]
  n.lm <- dim(lasec.result$fit)[2]
  plot(x=NA, xlim=c(0, n.lm), ylim=c(0,1), xlab="No. landmarks", ylab="Fit")
  for(i in 1:iter) {  
    par(new=T)
    plot(lasec.result$fit[i,], xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="grey", xlab="", ylab="", axes=F)   # 1-pss to match AWTY
  }
  par(new = T)
  plot(lasec.result$median.fit, xlim = c(0, n.lm), ylim = c(0, 1), type = "l", 
       col = "black", lwd = 3, xlab = "", ylab = "", axes = F)
}
