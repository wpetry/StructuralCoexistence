library(shiny)
library(rgl)

#################################################
## Structural coexistence calculation & plotting code
## from Saavedra et al. 2017 Ecol Mono
## http://dx.doi.org/10.1002/ecm.1263
#################################################
## excerpt from: toolbox_coexistence.R
#################################################
require(mvtnorm)

#input parameters:
#alpha = competition strength matrix 
#r = vector of intrinsic growth rates

#structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

#vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}

#structural fitness difference (in degree)
theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(round(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))),15))*180/pi  # add rounding step to overcome floating point errors
  return(out)
}

#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}

#test which pairs in a system (alpha and r) are feasible (output 1 = feasible, 0 = not feasible)
test_feasibility_pairs <- function(alpha,r){
  n <- length(r)
  c <- combn(n,2)
  nc <- dim(c)[2]
  f <- rep(NA,nc)
  for (i in 1:nc){
    f[i] <- prod(solve(alpha[c[,i],c[,i]],r[c[,i]])>0)
  }
  out <- list(pairs = c, feasibility = f)
  return(out)
}

#compute the feasiblity domain, the feasibility domain of all pairs, and their overlap (Nrand = number of randomization)
compute_overlap <- function(alpha,Nrand){
  
  n <- dim(alpha)[1]
  
  counter_f <- 0
  counter_overlap <- 0
  counter_all <- 0
  
  for (i in 1:Nrand){
    
    r_rand <- abs(rnorm(n))  
    r_rand <- r_rand/sqrt(sum(r_rand^2))
    
    f1 <- test_feasibility(alpha,r_rand)  
    f2 <- test_feasibility_pairs(alpha,r_rand)$feasibility  
    
    counter_f <- counter_f + f1
    counter_all <- counter_all + prod(f2)
    counter_overlap <- counter_overlap + f1*prod(f2)
    
  }
  
  Omega <- counter_f/Nrand
  Omega_all <- counter_all/Nrand
  overlap <- counter_overlap/Nrand
  
  out <- list(Omega = Omega, Omega_all = Omega_all, overlap = overlap)
  return(out)
  
}
#################################################
## exerpt from: toolbox_figure.R
#################################################
require('scatterplot3d')
require('pracma')

######################
#R-code to draw the feasibility domain and its projection on the simplex.
#inputs: alpha = interaction matrix, file_name = output name of the figure (.pdf),
#sp_names = the names of the species.

# Modified by WK Petry to plot centroid of the feasibility domain (r_c) and the
# current vector of intrinsic growth rates (r)

#######################
cone_3sp <- function(alpha){
  par(mar = c(0,0,0,0))
  
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  v3 <- alpha_n[,3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  
  lambda = c(0,1.2)
  
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  s3d <- scatterplot3d(X, -Y, Z, xlim=c(0,1.4), ylim=c(0,-1.4), zlim=c(0,1.4),type = 'l',
                       box = F, angle = 30, axis = F, grid = F, asp = 1)
  
  s3d$points3d(c(0,0),c(0,0),c(0,1.4), type = 'l', col = 'black', lwd = 2) 
  s3d$points3d(c(0,0),c(0,-1.4),c(0,0), type = 'l', col = 'black', lwd = 2) 
  s3d$points3d(c(0,1.4),c(0,0),c(0,0), type = 'l', col = 'black', lwd = 2) 
  
  pp <- s3d$xyz.convert(0.62,0,0)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp1",
       adj = c(0, 1.5), cex = 1)
  pp <- s3d$xyz.convert(0, -1.4, 0)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp2",
       adj = c(0, 1.5), cex = 1, srt = 31)
  pp <- s3d$xyz.convert(0, 0, 0.75)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp3",
       adj = c(0, 1.5), cex = 1, srt = 90)
  
  lambda = c(0,1.2)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  s3d$points3d(X,-Y,Z,type='l',col='mediumseagreen',lwd=4) 
  
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  s3d$points3d(X2,-Y2,Z2,type='l',col='mediumseagreen',lwd=4) 
  
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  s3d$points3d(X3,-Y3,Z3,type='l',col='mediumseagreen',lwd=4) 
  
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  s3d$points3d(X4,-Y4,Z4,type='l',col='blue4',lwd=4) 
  
  
  lambda = c(1.2,10)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  s3d$points3d(X,-Y,Z,type='l',col='mediumseagreen',lwd=4,lty=2) 
  
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  s3d$points3d(X2,-Y2,Z2,type='l',col='mediumseagreen',lwd=4,lty=2) 
  
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  s3d$points3d(X3,-Y3,Z3,type='l',col='mediumseagreen',lwd=4,lty=2) 
  
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  s3d$points3d(X4,-Y4,Z4,type='l',col='blue4',lwd=4,lty=2) 
  
  a <- seq(0,1,by=0.01)
  b <- sqrt(1-a^2)
  c <- rep(0,length(a))
  s3d$points3d(a*1.2,-b*1.2,c*1.2,type='l',col='grey50',lwd=2) 
  s3d$points3d(c*1.2,-a*1.2,b*1.2,type='l',col='grey50',lwd=2) 
  s3d$points3d(b*1.2,-c*1.2,a*1.2,type='l',col='grey50',lwd=2) 
  
  
  mu <- seq(0,1,by=0.01)
  w1 <- t(t(v1)) %*% t(mu) + t(t(v2)) %*% t(1-mu)
  w1 <- w1 %*% diag(1/sqrt(colSums(w1^2)))
  
  w2 <- t(t(v2)) %*% t(mu) + t(t(v3)) %*% t(1-mu)
  w2 <- w2 %*% diag(1/sqrt(colSums(w2^2)))
  
  w3 <- t(t(v3)) %*% t(mu) + t(t(v1)) %*% t(1-mu)
  w3 <- w3 %*% diag(1/sqrt(colSums(w3^2)))
  
  wp1 <- s3d$xyz.convert(w1[1,]*1.2,-w1[2,]*1.2,w1[3,]*1.2)
  wp2 <- s3d$xyz.convert(w2[1,]*1.2,-w2[2,]*1.2,w2[3,]*1.2)
  wp3 <- s3d$xyz.convert(w3[1,]*1.2,-w3[2,]*1.2,w3[3,]*1.2)
  
  XXX <- c(wp1$x, wp2$x, wp3$x)
  YYY <- c(wp1$y, wp2$y, wp3$y)
  
  color <- col2rgb("mediumseagreen")
  polygon(XXX,YYY,col= rgb(color[1,1],color[2,1],color[3,1],90,maxColorValue=255), border=FALSE)
  
  s3d$points3d(w1[1,]*1.2,-w1[2,]*1.2,w1[3,]*1.2,type='l',col='mediumseagreen',lwd=4) 
  s3d$points3d(w2[1,]*1.2,-w2[2,]*1.2,w2[3,]*1.2,type='l',col='mediumseagreen',lwd=4) 
  s3d$points3d(w3[1,]*1.2,-w3[2,]*1.2,w3[3,]*1.2,type='l',col='mediumseagreen',lwd=4) 
  
}

#######################

projection_3sp_with_pairwise <- function(alpha,r){
  par(mar = c(0,0,0,0))
  
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  v3 <- alpha_n[,3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  
  v1 <- v1/sum(v1)
  v2 <- v2/sum(v2)
  v3 <- v3/sum(v3)
  vc <- vc/sum(vc)
  
  Xf <- sqrt(2);
  Yf <- sqrt(6)/2;
  
  XX <- c(-Xf/2,Xf/2,0,-Xf/2)
  YY <- c(0,0,Yf,0)
  
  plot(-XX,YY,axes=F,xlab='',ylab='',xlim=c(-Xf/2-0.05,Xf/2+0.05),ylim=c(0-0.05,Yf+0.05),col='grey50',type='l',lwd=2, asp=1)
  
  #####
  
  v1P <- v1
  v1P[3] <- 0
  v1P <- v1P / sum(v1P)
  
  v2P <- v2
  v2P[3] <- 0
  v2P <- v2P / sum(v2P)
  
  vcP <- v1P/sqrt(sum(v1P^2)) + v2P/sqrt(sum(v2P^2))
  vcP[3] <- 0
  vcP <- vcP / sum(vcP)
  
  v1C <- c((0.5-0.5*v1P[3]-v1P[1])*Xf,v1P[3]*Yf)
  v2C <- c((0.5-0.5*v2P[3]-v2P[1])*Xf,v2P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf,vcP[3]*Yf)
  
  lines(-c(v1C[1],v2C[1]),c(v1C[2],v2C[2]),col='mediumseagreen',lwd=2)
  lines(-c(v1C[1],XX[3]),c(v1C[2],YY[3]),col='mediumseagreen',lty=2,lwd=2)
  lines(-c(v2C[1],XX[3]),c(v2C[2],YY[3]),col='mediumseagreen',lty=2,lwd=2)
  color <- col2rgb("mediumseagreen")
  polygon(-c(v1C[1],v2C[1],XX[3],v1C[1]),c(v1C[2],v2C[2],YY[3],v1C[2]),col=rgb(color[1,1],color[2,1],color[3,1],30,maxColorValue=255) ,border = F  )
  points(-c(v1C[1],v2C[1]),c(v1C[2],v2C[2]),col='dodgerblue',pch=16,cex=1.5)
  #points(-vcC[1],vcC[2],col='blue4',pch=16,cex=1.5)
  
  #####
  
  v1P <- v1
  v1P[2] <- 0
  v1P <- v1P / sum(v1P)
  
  v3P <- v3
  v3P[2] <- 0
  v3P <- v3P / sum(v3P)
  
  vcP <- v1P/sqrt(sum(v1P^2)) + v3P/sqrt(sum(v3P^2))
  vcP[2] <- 0
  vcP <- vcP / sum(vcP)
  
  v1C <- c((0.5-0.5*v1P[3]-v1P[1])*Xf,v1P[3]*Yf)
  v3C <- c((0.5-0.5*v3P[3]-v3P[1])*Xf,v3P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf,vcP[3]*Yf)
  lines(-c(v1C[1],v3C[1]),c(v1C[2],v3C[2]),col='mediumseagreen',lwd=2)
  lines(-c(v1C[1],XX[2]),c(v1C[2],YY[2]),col='mediumseagreen',lty=2,lwd=2)
  lines(-c(v3C[1],XX[2]),c(v3C[2],YY[2]),col='mediumseagreen',lty=2,lwd=2)
  polygon(-c(v1C[1],v3C[1],XX[2],v1C[1]),c(v1C[2],v3C[2],YY[2],v1C[2]),col=rgb(color[1,1],color[2,1],color[3,1],30,maxColorValue=255) ,border = F  )
  points(-c(v1C[1],v3C[1]),c(v1C[2],v3C[2]),col='dodgerblue',pch=16,cex=1.5)
  #points(-vcC[1],vcC[2],col='blue4',pch=16,cex=1.5)
  
  
  #####
  
  v2P <- v2
  v2P[1] <- 0
  v2P <- v2P / sum(v2P)
  
  v3P <- v3
  v3P[1] <- 0
  v3P <- v3P / sum(v3P)
  
  vcP <- v2P/sqrt(sum(v2P^2)) + v3P/sqrt(sum(v3P^2))
  vcP[1] <- 0
  vcP <- vcP / sum(vcP)
  
  
  v2C <- c((0.5-0.5*v2P[3]-v2P[1])*Xf,v2P[3]*Yf)
  v3C <- c((0.5-0.5*v3P[3]-v3P[1])*Xf,v3P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf,vcP[3]*Yf)
  
  lines(-c(v2C[1],v3C[1]),c(v2C[2],v3C[2]),col='mediumseagreen',lwd=2)
  lines(-c(v2C[1],XX[1]),c(v2C[2],YY[1]),col='mediumseagreen',lty=2,lwd=2)
  lines(-c(v3C[1],XX[1]),c(v3C[2],YY[1]),col='mediumseagreen',lty=2,lwd=2)
  polygon(-c(v2C[1],v3C[1],XX[1],v2C[1]),c(v2C[2],v3C[2],YY[1],v2C[2]),col=rgb(color[1,1],color[2,1],color[3,1],30,maxColorValue=255) ,border = F  )
  points(-c(v2C[1],v3C[1]),c(v2C[2],v3C[2]),col='dodgerblue',pch=16,cex=1.5)
  
  
  
  v1C <- c((0.5-0.5*v1[3]-v1[1])*Xf,v1[3]*Yf)
  v2C <- c((0.5-0.5*v2[3]-v2[1])*Xf,v2[3]*Yf)
  v3C <- c((0.5-0.5*v3[3]-v3[1])*Xf,v3[3]*Yf)
  vcC <- c((0.5-0.5*vc[3]-vc[1])*Xf,vc[3]*Yf)
  
  color <- col2rgb("green4")
  
  polygon(-c(v1C[1],v2C[1],v3C[1]),c(v1C[2],v2C[2],v3C[2]),col= rgb(color[1,1],color[2,1],color[3,1],90,maxColorValue=255), border=FALSE)
  
  points(-c(v1C[1],v2C[1],v3C[1],v1C[1]) , c(v1C[2],v2C[2],v3C[2],v1C[2]), col= 'green4', type='l',cex=1.5,lwd=2)
  points(-c(v1C[1],v2C[1],v3C[1]),c(v1C[2],v2C[2],v3C[2]),col='blue4',pch=16,cex=1.5)
  
  points(-vcC[1], vcC[2], col="blue4", pch=1, cex=1.5) # plot centroid of D_F
  
  rX <- Xf*(0.5-0.5*((2*r[2]+r[3])/(r[1]+r[2]+r[3])))
  rY <- Yf*(r[3]/(r[1]+r[2]+r[3]))
  
  points(rX, rY, col= "black", pch=16,cex=1.5)
  
  text(-XX,YY,labels = c("sp1","sp2","sp3"),cex=1.7,pos = c(1,1,3))
  
}

# ported to rgl for interactive 3d plotting

require(rgl)


plot_cone_3D <- function(alpha, 
                         r = c(0,0,0),                    # if no vector of growth rates provided, null vector 
                         sp_name = c('sp1','sp2','sp3')){ # Species names can be customized
  
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  v3 <- alpha_n[,3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  
  lambda = c(0,1.2)
  
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  
  lambda = c(0,1.2)
  Sp1 <- v1[1] * lambda 
  Sp2 <- v1[2] * lambda 
  Sp3 <- v1[3] * lambda 
  
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  
  # feasibility domain
  
  plot3d(Sp1, -Sp2, Sp3, col = 'mediumseagreen',
         xlab = "", ylab = "", zlab = "",
         type = 'l', lwd = 2.5, box = FALSE, axes = FALSE)
  lines3d(X2, -Y2, Z2, col = 'mediumseagreen', lwd = 2.5)
  lines3d(X3, -Y3, Z3, col = 'mediumseagreen', lwd = 2.5)
  lines3d(X4, -Y4, Z4, col = 'orange', lwd = 2.5)
  
  
  lambda = c(1.2,2)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  
  lines3d(X,-Y,Z, col = 'mediumseagreen', lwd = 1)
  lines3d(X2,-Y2,Z2, col = 'mediumseagreen', lwd = 1)
  lines3d(X3,-Y3,Z3, col = 'mediumseagreen', lwd = 1)
  lines3d(X4,-Y4,Z4, col = 'orange', lwd = 1)
  
  
  #axes
  
  lines3d(c(0,1.4),c(0,0),c(0,0), col = 'black', lwd = 1)
  lines3d(c(0,0),c(0,-1.4),c(0,0), col = 'black', lwd = 1)
  lines3d(c(0,0),c(0,0),c(0,1.4), col = 'black', lwd = 1)
  
  #arcs
  
  a <- seq(0,1,by=0.01)
  b <- sqrt(1-a^2)
  c <- rep(0,length(a))
  lines3d(a*1.2,-b*1.2,c*1.2,col='grey', lwd = 2.5) 
  lines3d(c*1.2,-a*1.2,b*1.2,col='grey', lwd = 2.5) 
  lines3d(b*1.2,-c*1.2,a*1.2,col='grey', lwd = 2.5) 
  
  mu <- seq(0,1,by=0.01)
  w1 <- t(t(v1)) %*% t(mu) + t(t(v2)) %*% t(1-mu)
  w1 <- w1 %*% diag(1/sqrt(colSums(w1^2)))
  
  w2 <- t(t(v2)) %*% t(mu) + t(t(v3)) %*% t(1-mu)
  w2 <- w2 %*% diag(1/sqrt(colSums(w2^2)))
  
  w3 <- t(t(v3)) %*% t(mu) + t(t(v1)) %*% t(1-mu)
  w3 <- w3 %*% diag(1/sqrt(colSums(w3^2)))
  
  lines3d(w1[1,]*1.2,-w1[2,]*1.2,w1[3,]*1.2,col='mediumseagreen',lwd=2) 
  lines3d(w2[1,]*1.2,-w2[2,]*1.2,w2[3,]*1.2,col='mediumseagreen',lwd=2) 
  lines3d(w3[1,]*1.2,-w3[2,]*1.2,w3[3,]*1.2,col='mediumseagreen',lwd=2) 
  
  ## Surface of the conical hull (yet to implement in rgl)
  
  # wp1 <- s3d$xyz.convert(w1[1,]*1.2,-w1[2,]*1.2,w1[3,]*1.2) 
  # wp2 <- s3d$xyz.convert(w2[1,]*1.2,-w2[2,]*1.2,w2[3,]*1.2)
  # wp3 <- s3d$xyz.convert(w3[1,]*1.2,-w3[2,]*1.2,w3[3,]*1.2)
  # 
  # XXX <- c(wp1$x, wp2$x, wp3$x)
  # YYY <- c(wp1$y, wp2$y, wp3$y)
  
  # color <- col2rgb("mediumseagreen")
  # polygon(XXX,YYY,col= rgb(color[1,1],color[2,1],color[3,1],90,maxColorValue=255), border=FALSE)
  
  # vector of growth rates
  
  rs <- 1.2*r/sqrt(r[1]^2+r[2]^2+r[3]^2)
  lines3d(c(0,rs[1]),c(0,-rs[2]),c(0,rs[3]), col = 'red', lwd = 2)
  rs2 <- 1.8*r/sqrt(r[1]^2+r[2]^2+r[3]^2)
  lines3d(c(0,rs2[1]),c(0,-rs2[2]),c(0,rs2[3]), col = 'red', lwd = 0.8)
  
  text3d(1.5,0,0, text = sp_name[1], col = 'black', lwd = 0.8)
  text3d(0,-1.5,0, text = sp_name[2], col = 'black', lwd = 0.8)
  text3d(0,0,1.5, text = sp_name[3], col = 'black', lwd = 0.8)
  
  aspect3d("iso")
}




#################################################
# Shiny server logic
#################################################
# This section takes inputs from the ui and displays the desired
# information using the code from Saavedra et al. above

shinyServer(function(input, output) {
  output$stats <- renderTable({
    r <- c(input$r1, input$r2, input$r3)
    feas_pairs <- test_feasibility_pairs(input$alphamat,r)
    feas_txt <- ifelse(sum(feas_pairs$feasibility)==0L,
                       "none",
                       ifelse(sum(feas_pairs$feasibility)==1L,
                              paste(feas_pairs$pairs[,which(feas_pairs$feasibility==1L)],
                                    collapse=" & "),
                              paste(apply(feas_pairs$pairs[,which(feas_pairs$feasibility==1L)],
                                          MARGIN=2,FUN=paste,collapse=" & "),collapse=", ")))
    data.frame(Metric=c("Structural niche difference (Ω)",
                        "Structural fitness difference (θ)",
                        HTML("Centroid of feasibility domain (r_c)"),
                        "Feasible triplet?",
                        "Feasible species pairs"),
               Value=c(paste0(round(exp(Omega(alpha=input$alphamat)),3),"㏛"),
                       paste0(round(theta(input$alphamat,r),3),"°"),
                       paste(round(r_centroid(input$alphamat),3),collapse=", "),
                       ifelse(test_feasibility(input$alphamat,r)==1L,"yes","no"),
                       feas_txt))
  })
  output$cone <- renderPlot({
    cone_3sp(input$alphamat)
  })
  output$proj <- renderPlot({
    projection_3sp_with_pairwise(input$alphamat,r=c(input$r1, input$r2, input$r3))
  })
  output$cone3d <- renderRglwidget({
    open3d(useNULL = T)
    plot_cone_3D(input$alphamat, r = c(input$r1, input$r2, input$r3))
    rglwidget()
  })
  
})
