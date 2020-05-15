#################################################-
## SERVER
#################################################-
## Preliminaries ----
#################################################-
library(shiny)
library(rgl)
library(scatterplot3d)
library(pracma)
library(mvtnorm)
library(igraph)
#################################################-
## Functions: calculate metrics ----
## Structural coexistence calculation & plotting code
## from Saavedra et al. 2017 Ecol Mono
## http://dx.doi.org/10.1002/ecm.1263
#################################################-
# Omega ====
# structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0, n), upper = rep(Inf, n),
               mean = rep(0, n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

# r_centroid ====
# vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- unname(t(t(r_c)))
  return(r_c)
}

# theta ====
# structural fitness difference (in degree)
theta <- function(alpha, r){
  r_c <- r_centroid(alpha)
  out <- acos(round(sum(r_c*r) / (sqrt(sum(r^2)) * sqrt(sum(r_c^2))), 15)) * 180 / pi  # add rounding step to overcome floating point errors
  return(out)
}

# test_feasibility ====
# test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha, r){
  out <- prod(solve(alpha, r) > 0)
  return(out)
}

# test_feasibility_pairs ====
# test which pairs in a system (alpha and r) are feasible (output 1 = feasible, 0 = not feasible)
test_feasibility_pairs <- function(alpha, r){
  n <- length(r)
  c <- combn(n, 2)
  nc <- dim(c)[2]
  f <- rep(NA, nc)
  for (i in 1:nc){
    f[i] <- prod(solve(alpha[c[, i], c[, i]], r[c[, i]]) > 0)
  }
  out <- list(pairs = c, feasibility = f)
  return(out)
}

# compute_overlap ====
# compute the feasiblity domain, the feasibility domain of all pairs, and their overlap (Nrand = number of randomization)
compute_overlap <- function(alpha, Nrand){
  n <- dim(alpha)[1]
  counter_f <- 0
  counter_overlap <- 0
  counter_all <- 0
  for (i in 1:Nrand){
    r_rand <- abs(rnorm(n))  
    r_rand <- r_rand / sqrt(sum(r_rand^2))
    f1 <- test_feasibility(alpha, r_rand)  
    f2 <- test_feasibility_pairs(alpha, r_rand)$feasibility
    counter_f <- counter_f + f1
    counter_all <- counter_all + prod(f2)
    counter_overlap <- counter_overlap + f1*prod(f2)
  }
  Omega <- counter_f / Nrand
  Omega_all <- counter_all / Nrand
  overlap <- counter_overlap / Nrand
  out <- list(Omega = Omega, Omega_all = Omega_all,
              overlap = overlap)
  return(out)
}

#################################################-
## Functions: make plots ----
#################################################-
# mk_graph_3sp ====
# make interaction network plot
mk_graph_3sp <- function(alphamat, rs){
  alphamat <- unname(alphamat)
  g <- igraph::graph_from_adjacency_matrix(alphamat > 0)
  E(g)$weight <- as.numeric(alphamat)
  widths <- E(g)$weight * 5
  #widths[widths > 1] <- sqrt(widths)
  plot(g,
       main = "Species interaction network",
       margin = c(0, -0.15, -0.3, -0.15),
       xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
       vertex.label.cex = 2.5,
       vertex.label.color = "black",
       vertex.size = 50 * rs,
       vertex.color = "grey80",
       vertex.frame.color = "transparent",
       edge.curved = TRUE,
       edge.width = widths,
       edge.arrow.size = 3,
       edge.arrow.mode = c(0, 2, 2,
                           2, 0, 2,
                           2, 2, 0),
       edge.color = "black",
       edge.loop.angle = 0.75,
       layout = matrix(c(4, 0, 0, 0, 2, sqrt(3)/2), ncol = 2,
                       byrow = TRUE))
}

# mk_graph_4sp ====
# make interaction network plot
mk_graph_4sp <- function(alphamat, rs){
  alphamat <- unname(alphamat)
  g <- igraph::graph_from_adjacency_matrix(alphamat > 0)
  E(g)$weight <- as.numeric(alphamat)
  widths <- E(g)$weight * 5
  #widths[widths > 1] <- sqrt(widths)
  plot(g,
       main = "Species interaction network",
       margin = c(0, -0.15, -0.3, -0.15),
       xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
       vertex.label.cex = 2.5,
       vertex.label.color = "black",
       vertex.size = 50 * rs,
       vertex.color = "grey80",
       vertex.frame.color = "transparent",
       edge.curved = TRUE,
       edge.width = widths,
       edge.arrow.size = 3,
       edge.arrow.mode = c(0, 2, 2, 2,
                           2, 0, 2, 2,
                           2, 2, 0, 2,
                           2, 2, 2, 0),
       edge.color = "black",
       edge.loop.angle = 0.5,
       layout = matrix(c(4, 0,
                         0, 0,
                         0, 4,
                         4, 4), ncol = 2,
                       byrow = TRUE))
}

# cone_3sp ====
# (helper) make static feasibility cone plot for 3 species
cone_3sp <- function(alpha){
  par(mar = c(0, 0, 0, 0))
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[, 1]
  v2 <- alpha_n[, 2]
  v3 <- alpha_n[, 3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  lambda = c(0, 1.2)
  X <- v1[1] * lambda
  Y <- v1[2] * lambda
  Z <- v1[3] * lambda
  s3d <- scatterplot3d(X, -Y, Z, xlim = c(0, 1.4),
                       ylim = c(0, -1.4), zlim = c(0, 1.4),
                       type = 'l', box = FALSE, angle = 30,
                       axis = FALSE, grid = FALSE, asp = 1)
  s3d$points3d(c(0, 0), c(0, 0), c(0, 1.4), type = 'l',
               col = "black", lwd = 2)
  s3d$points3d(c(0, 0), c(0, -1.4), c(0, 0), type = 'l',
               col = "black", lwd = 2)
  s3d$points3d(c(0, 1.4), c(0, 0), c(0, 0), type = 'l',
               col = "black", lwd = 2)
  pp <- s3d$xyz.convert(0.62,0,0)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp1",
       adj = c(0, 1.5), cex = 1)
  pp <- s3d$xyz.convert(0, -1.4, 0)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp2",
       adj = c(0, 1.5), cex = 1, srt = 31)
  pp <- s3d$xyz.convert(0, 0, 0.75)
  text(x = pp$x,y = pp$y, labels = "Intrinsic growth rate sp3",
       adj = c(0, 1.5), cex = 1, srt = 90)
  lambda = c(0, 1.2)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  s3d$points3d(X, -Y, Z, type = 'l', col = "mediumseagreen",
               lwd = 4)
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  s3d$points3d(X2, -Y2, Z2, type = 'l', col = "mediumseagreen",
               lwd = 4)
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  s3d$points3d(X3, -Y3, Z3, type = 'l', col = "mediumseagreen",
               lwd = 4)
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  s3d$points3d(X4, -Y4, Z4, type = 'l', col = "blue4", lwd = 4)
  lambda <- c(1.2, 10)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
  s3d$points3d(X, -Y, Z, type = 'l', col = "mediumseagreen",
               lwd = 4, lty = 2)
  X2 <- v2[1] * lambda
  Y2 <- v2[2] * lambda
  Z2 <- v2[3] * lambda
  s3d$points3d(X2, -Y2, Z2, type = 'l', col = "mediumseagreen",
               lwd = 4, lty = 2)
  X3 <- v3[1] * lambda
  Y3 <- v3[2] * lambda
  Z3 <- v3[3] * lambda
  s3d$points3d(X3, -Y3, Z3, type = 'l', col = "mediumseagreen",
               lwd = 4, lty = 2)
  X4 <- vc[1] * lambda
  Y4 <- vc[2] * lambda
  Z4 <- vc[3] * lambda
  s3d$points3d(X4, -Y4, Z4, type = 'l', col = "blue4", lwd = 4,
               lty = 2)
  a <- seq(0, 1, by = 0.01)
  b <- sqrt(1-a^2)
  c <- rep(0, length(a))
  s3d$points3d(a*1.2, -b*1.2, c*1.2, type = 'l', col = "grey50",
               lwd = 2)
  s3d$points3d(c*1.2, -a*1.2, b*1.2, type = 'l', col = "grey50",
               lwd = 2)
  s3d$points3d(b*1.2, -c*1.2, a*1.2, type = 'l', col = "grey50",
               lwd = 2)
  mu <- seq(0, 1, by = 0.01)
  w1 <- t(t(v1)) %*% t(mu) + t(t(v2)) %*% t(1-mu)
  w1 <- w1 %*% diag(1/sqrt(colSums(w1^2)))
  w2 <- t(t(v2)) %*% t(mu) + t(t(v3)) %*% t(1-mu)
  w2 <- w2 %*% diag(1/sqrt(colSums(w2^2)))
  w3 <- t(t(v3)) %*% t(mu) + t(t(v1)) %*% t(1-mu)
  w3 <- w3 %*% diag(1/sqrt(colSums(w3^2)))
  wp1 <- s3d$xyz.convert(w1[1, ]*1.2, -w1[2, ]*1.2, w1[3, ]*1.2)
  wp2 <- s3d$xyz.convert(w2[1, ]*1.2, -w2[2, ]*1.2, w2[3, ]*1.2)
  wp3 <- s3d$xyz.convert(w3[1, ]*1.2, -w3[2, ]*1.2, w3[3, ]*1.2)
  XXX <- c(wp1$x, wp2$x, wp3$x)
  YYY <- c(wp1$y, wp2$y, wp3$y)
  color <- col2rgb("mediumseagreen")
  polygon(XXX, YYY,
          col = rgb(color[1, 1], color[2, 1], color[3, 1], 90,
                    maxColorValue = 255), border = FALSE)
  s3d$points3d(w1[1, ]*1.2, -w1[2, ]*1.2, w1[3, ]*1.2,
               type = 'l', col = "mediumseagreen", lwd = 4)
  s3d$points3d(w2[1, ]*1.2, -w2[2, ]*1.2, w2[3, ]*1.2,
               type = 'l', col = "mediumseagreen", lwd = 4)
  s3d$points3d(w3[1, ]*1.2, -w3[2, ]*1.2, w3[3, ]*1.2,
               type = 'l', col = "mediumseagreen", lwd = 4)
}

# projection_3sp_with_pairwise ====
# plot static projection for 3 specie community
projection_3sp_with_pairwise <- function(alpha, r){
  par(mar = c(0, 0, 0, 0))
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[, 1]
  v2 <- alpha_n[, 2]
  v3 <- alpha_n[, 3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  v3 <- v3 / sum(v3)
  vc <- vc / sum(vc)
  Xf <- sqrt(2)
  Yf <- sqrt(6) / 2
  XX <- c(-Xf / 2, Xf / 2, 0, -Xf / 2)
  YY <- c(0, 0, Yf, 0)
  plot(-XX, YY, axes = FALSE, xlab = '', ylab = '',
       xlim = c(-Xf/2-0.05, Xf/2+0.05),
       ylim = c(0-0.05, Yf+0.05),
       col = "grey50", type = "l", lwd = 2, asp = 1)
  v1P <- v1
  v1P[3] <- 0
  v1P <- v1P / sum(v1P)
  v2P <- v2
  v2P[3] <- 0
  v2P <- v2P / sum(v2P)
  vcP <- v1P/sqrt(sum(v1P^2)) + v2P/sqrt(sum(v2P^2))
  vcP[3] <- 0
  vcP <- vcP / sum(vcP)
  v1C <- c((0.5-0.5*v1P[3]-v1P[1])*Xf, v1P[3]*Yf)
  v2C <- c((0.5-0.5*v2P[3]-v2P[1])*Xf, v2P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf, vcP[3]*Yf)
  lines(-c(v1C[1], v2C[1]), c(v1C[2], v2C[2]),
        col = "mediumseagreen", lwd = 1)
  lines(-c(v1C[1], XX[3]), c(v1C[2], YY[3]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  lines(-c(v2C[1], XX[3]), c(v2C[2], YY[3]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  color <- col2rgb("mediumseagreen")
  polygon(-c(v1C[1], v2C[1], XX[3], v1C[1]),
          c(v1C[2], v2C[2], YY[3], v1C[2]),
          col = rgb(color[1, 1], color[2, 1], color[3, 1],
                    30, maxColorValue = 255), border = FALSE)
  points(-c(v1C[1], v2C[1]), c(v1C[2], v2C[2]),
         col = "dodgerblue", pch = 16, cex = 1.5)
  v1P <- v1
  v1P[2] <- 0
  v1P <- v1P / sum(v1P)
  v3P <- v3
  v3P[2] <- 0
  v3P <- v3P / sum(v3P)
  vcP <- v1P/sqrt(sum(v1P^2)) + v3P/sqrt(sum(v3P^2))
  vcP[2] <- 0
  vcP <- vcP / sum(vcP)
  v1C <- c((0.5-0.5*v1P[3]-v1P[1])*Xf, v1P[3]*Yf)
  v3C <- c((0.5-0.5*v3P[3]-v3P[1])*Xf, v3P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf, vcP[3]*Yf)
  lines(-c(v1C[1], v3C[1]), c(v1C[2], v3C[2]),
        col = "mediumseagreen", lwd = 1)
  lines(-c(v1C[1], XX[2]), c(v1C[2], YY[2]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  lines(-c(v3C[1], XX[2]), c(v3C[2], YY[2]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  polygon(-c(v1C[1], v3C[1], XX[2], v1C[1]),
          c(v1C[2], v3C[2], YY[2], v1C[2]),
          col = rgb(color[1, 1], color[2, 1], color[3, 1],
                    30, maxColorValue = 255), border = FALSE)
  points(-c(v1C[1], v3C[1]), c(v1C[2], v3C[2]),
         col = "dodgerblue", pch = 16, cex = 1.5)
  # points(-vcC[1], vcC[2], col = "blue4", pch = 16, cex = 1.5)
  v2P <- v2
  v2P[1] <- 0
  v2P <- v2P / sum(v2P)
  v3P <- v3
  v3P[1] <- 0
  v3P <- v3P / sum(v3P)
  vcP <- v2P/sqrt(sum(v2P^2)) + v3P/sqrt(sum(v3P^2))
  vcP[1] <- 0
  vcP <- vcP / sum(vcP)
  v2C <- c((0.5-0.5*v2P[3]-v2P[1])*Xf, v2P[3]*Yf)
  v3C <- c((0.5-0.5*v3P[3]-v3P[1])*Xf, v3P[3]*Yf)
  vcC <- c((0.5-0.5*vcP[3]-vcP[1])*Xf, vcP[3]*Yf)
  lines(-c(v2C[1], v3C[1]), c(v2C[2], v3C[2]),
        col = "mediumseagreen", lwd = 1)
  lines(-c(v2C[1], XX[1]), c(v2C[2], YY[1]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  lines(-c(v3C[1], XX[1]), c(v3C[2], YY[1]),
        col = "mediumseagreen", lty = 2, lwd = 1)
  polygon(-c(v2C[1], v3C[1], XX[1], v2C[1]),
          c(v2C[2], v3C[2], YY[1], v2C[2]),
          col = rgb(color[1, 1], color[2, 1], color[3, 1],
                    30, maxColorValue = 255), border = FALSE)
  points(-c(v2C[1], v3C[1]), c(v2C[2], v3C[2]),
         col = "dodgerblue", pch = 16, cex = 1.5)
  v1C <- c((0.5-0.5*v1[3]-v1[1])*Xf, v1[3]*Yf)
  v2C <- c((0.5-0.5*v2[3]-v2[1])*Xf, v2[3]*Yf)
  v3C <- c((0.5-0.5*v3[3]-v3[1])*Xf, v3[3]*Yf)
  vcC <- c((0.5-0.5*vc[3]-vc[1])*Xf, vc[3]*Yf)
  color <- col2rgb("green4")
  polygon(-c(v1C[1], v2C[1], v3C[1]), c(v1C[2], v2C[2], v3C[2]),
          col = rgb(color[1, 1], color[2, 1], color[3, 1],
                    90, maxColorValue = 255), border = FALSE)
  points(-c(v1C[1], v2C[1], v3C[1], v1C[1]),
         c(v1C[2], v2C[2], v3C[2], v1C[2]), col = "green4",
         type = 'l', cex = 1.5, lwd = 2)
  points(-c(v1C[1], v2C[1], v3C[1]), c(v1C[2], v2C[2], v3C[2]),
         col = "blue4", pch = 16, cex = 1.5)
  points(-vcC[1], vcC[2], col = "blue4", pch = 4, cex = 1.5) # plot centroid of D_F
  rX <- Xf*(0.5-0.5*((2*r[2]+r[3]) / (r[1]+r[2]+r[3])))
  rY <- Yf*(r[3] / (r[1]+r[2]+r[3]))
  points(rX, rY, col = "black", pch = 21, lwd = 4, cex = 3)
  text(-XX, YY, labels = c("sp1", "sp2", "sp3"), cex = 1.7,
       pos = c(1, 1, 3))
}

# plot_cone_3D ====
# plot interactive (rgl) cone for 3 species community
plot_cone_3D <- function(alpha, r = c(0, 0, 0),
                         sp_name = paste0("sp", 1:3)){
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[, 1]
  v2 <- alpha_n[, 2]
  v3 <- alpha_n[, 3]
  vc <- (v1 + v2 + v3)
  vc <- vc / sqrt(sum(vc^2))
  lambda <- c(0, 1.2)
  X <- v1[1] * lambda 
  Y <- v1[2] * lambda 
  Z <- v1[3] * lambda 
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
  rgl::plot3d(Sp1, -Sp2, Sp3, col = "mediumseagreen",
         xlab = "", ylab = "", zlab = "",
         type = 'l', lwd = 2.5, box = FALSE, axes = FALSE)
  rgl::lines3d(X2, -Y2, Z2, col = "mediumseagreen", lwd = 2.5)
  rgl::lines3d(X3, -Y3, Z3, col = "mediumseagreen", lwd = 2.5)
  rgl::lines3d(X4, -Y4, Z4, col = "orange", lwd = 2.5)
  lambda <- c(1.2, 2)
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
  rgl::lines3d(X,-Y,Z, col = "mediumseagreen", lwd = 1)
  rgl::lines3d(X2,-Y2,Z2, col = "mediumseagreen", lwd = 1)
  rgl::lines3d(X3,-Y3,Z3, col = "mediumseagreen", lwd = 1)
  rgl::lines3d(X4,-Y4,Z4, col = "orange", lwd = 1)
  #axes
  rgl::lines3d(c(0, 1.4), c(0, 0), c(0, 0), col = "black", lwd = 1)
  rgl::lines3d(c(0, 0), c(0, -1.4), c(0, 0), col = "black", lwd = 1)
  rgl::lines3d(c(0, 0), c(0, 0), c(0, 1.4), col = "black", lwd = 1)
  #arcs
  a <- seq(0, 1, by = 0.01)
  b <- sqrt(1-a^2)
  c <- rep(0, length(a))
  rgl::lines3d(a*1.2, -b*1.2, c*1.2, col = "grey", lwd = 2.5)
  rgl::lines3d(c*1.2, -a*1.2, b*1.2, col = "grey", lwd = 2.5)
  rgl::lines3d(b*1.2, -c*1.2, a*1.2, col = "grey", lwd = 2.5)
  mu <- seq(0, 1, by = 0.01)
  w1 <- t(t(v1)) %*% t(mu) + t(t(v2)) %*% t(1-mu)
  w1 <- w1 %*% diag(1 / sqrt(colSums(w1^2)))
  w2 <- t(t(v2)) %*% t(mu) + t(t(v3)) %*% t(1-mu)
  w2 <- w2 %*% diag(1 / sqrt(colSums(w2^2)))
  w3 <- t(t(v3)) %*% t(mu) + t(t(v1)) %*% t(1-mu)
  w3 <- w3 %*% diag(1 / sqrt(colSums(w3^2)))
  rgl::lines3d(w1[1, ]*1.2, -w1[2, ]*1.2, w1[3, ]*1.2,
               col = "mediumseagreen", lwd = 2) 
  rgl::lines3d(w2[1, ]*1.2, -w2[2, ]*1.2, w2[3, ]*1.2,
               col = "mediumseagreen", lwd = 2) 
  rgl::lines3d(w3[1, ]*1.2, -w3[2, ]*1.2, w3[3, ]*1.2,
               col = "mediumseagreen", lwd = 2) 
  # surface of the conical hull (yet to implement in rgl)
  # wp1 <- s3d$xyz.convert(w1[1, ]*1.2, -w1[2, ]*1.2, w1[3, ]*1.2)
  # wp2 <- s3d$xyz.convert(w2[1, ]*1.2, -w2[2, ]*1.2, w2[3, ]*1.2)
  # wp3 <- s3d$xyz.convert(w3[1, ]*1.2, -w3[2, ]*1.2, w3[3, ]*1.2)
  # XXX <- c(wp1$x, wp2$x, wp3$x)
  # YYY <- c(wp1$y, wp2$y, wp3$y)
  # color <- col2rgb("mediumseagreen")
  # polygon(XXX, YYY, col = rgb(color[1, 1], color[2, 1], color[3, 1],
  # 90, maxColorValue = 255), border = FALSE)
  # vector of growth rates
  rs <- 1.2*r / sqrt(r[1]^2+r[2]^2+r[3]^2)
  rgl::lines3d(c(0, rs[1]), c(0, -rs[2]), c(0, rs[3]), col = "red",
               lwd = 2)
  rs2 <- 1.8*r / sqrt(r[1]^2+r[2]^2+r[3]^2)
  rgl::lines3d(c(0, rs2[1]), c(0, -rs2[2]), c(0, rs2[3]), col = "red",
               lwd = 0.8)
  rgl::text3d(1.5, 0, 0, text = sp_name[1], col = "black", lwd = 0.8)
  rgl::text3d(0, -1.5, 0, text = sp_name[2], col = "black", lwd = 0.8)
  rgl::text3d(0, 0, 1.5, text = sp_name[3], col = "black", lwd = 0.8)
  rgl::aspect3d("iso")
}

# projection_4sp_3D ====
# 4 species representation projected on the 3-simplex
projection_4sp_3D <- function(alpha, r, sp_name = paste0("sp", 1:4)){
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  v1 <- alpha_n[,1]
  v2 <- alpha_n[,2]
  v3 <- alpha_n[,3]
  v4 <- alpha_n[,4]
  vc <- (v1 + v2 + v3 + v4)
  vc <- vc / sqrt(sum(vc^2))
  vr <- r / sqrt(sum(r^2))
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  v3 <- v3 / sum(v3)
  v4 <- v4 / sum(v4)
  vc <- vc / sum(vc)
  vr <- vr / sum(vr)
  t1 <- c(1, 0, 0, 0)
  t2 <- c(0, 1, 0, 0)
  t3 <- c(0, 0, 1, 0)
  t4 <- c(0, 0, 0, 1)
  e1 <- c(1, 0, 0, 0)
  e2 <- c(0, 1, 0, 0)
  e3 <- c(0, 0, 1, 0)
  e4 <- c(0, 0, 0, 1)
  w2 <- t(t(e2-e1))
  w3 <- t(t(e3-e1))
  w4 <- t(t(e4-e1))
  A <- pracma::gramSchmidt(cbind(w2, w3, w4))
  TT <- A$R
  v1 <- v1[-1]
  v2 <- v2[-1]
  v3 <- v3[-1]
  v4 <- v4[-1]
  vc <- vc[-1]
  vr <- vr[-1]
  t1 <- t1[-1]
  t2 <- t2[-1]
  t3 <- t3[-1]
  t4 <- t4[-1]
  v1 <- as.vector(TT %*% v1)
  v2 <- as.vector(TT %*% v2)
  v3 <- as.vector(TT %*% v3)
  v4 <- as.vector(TT %*% v4)
  vc <- as.vector(TT %*% vc)
  vr <- as.vector(TT %*% vr)
  t1 <- as.vector(TT %*% t1)
  t2 <- as.vector(TT %*% t2)
  t3 <- as.vector(TT %*% t3)
  t4 <- as.vector(TT %*% t4)
  e1 <- c(0, 0, 0)
  e2 <- c(1, 0, 0)
  e3 <- c(0, 1, 0)
  e4 <- c(0, 0, 1)
  e1 <- as.vector(TT %*% e1)
  e2 <- as.vector(TT %*% e2)
  e3 <- as.vector(TT %*% e3)
  e4 <- as.vector(TT %*% e4)
  X <- c(e1[1], e2[1])
  Y <- c(e1[2], e2[2])
  Z <- c(e1[3], e2[3])
  rgl::plot3d(X, Y, Z, col = "grey", 
              xlab = "", ylab = "", zlab = "", type = 'l', lwd = 2,
              box = FALSE, axes = FALSE)
  X <- c(e1[1], e3[1])
  Y <- c(e1[2], e3[2])
  Z <- c(e1[3], e3[3])
  rgl::lines3d(X, Y, Z, col = "grey", lwd = 2)
  X <- c(e1[1], e4[1])
  Y <- c(e1[2], e4[2])
  Z <- c(e1[3], e4[3])
  rgl::lines3d(X, Y, Z, col = "grey", lwd = 2)
  X <- c(e2[1], e3[1])
  Y <- c(e2[2], e3[2])
  Z <- c(e2[3], e3[3])
  rgl::lines3d(X, Y, Z, col = "grey", lwd = 2)
  X <- c(e2[1], e4[1])
  Y <- c(e2[2], e4[2])
  Z <- c(e2[3], e4[3])
  rgl::lines3d(X, Y, Z, col = "grey", lwd = 2)
  X <- c(e3[1], e4[1])
  Y <- c(e3[2], e4[2])
  Z <- c(e3[3], e4[3])
  rgl::lines3d(X, Y, Z, col = "grey", lwd = 2)
  X <- c(v1[1], v2[1])
  Y <- c(v1[2], v2[2])
  Z <- c(v1[3], v2[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- c(v1[1], v3[1])
  Y <- c(v1[2], v3[2])
  Z <- c(v1[3], v3[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- c(v1[1], v4[1])
  Y <- c(v1[2], v4[2])
  Z <- c(v1[3], v4[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- c(v2[1], v3[1])
  Y <- c(v2[2], v3[2])
  Z <- c(v2[3], v3[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- c(v2[1], v4[1])
  Y <- c(v2[2], v4[2])
  Z <- c(v2[3], v4[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- c(v3[1], v4[1])
  Y <- c(v3[2], v4[2])
  Z <- c(v3[3], v4[3])
  rgl::lines3d(X, Y, Z, col = "darkgreen", lwd = 2)
  X <- v1[1]
  Y <- v1[2]
  Z <- v1[3]
  rgl::points3d(X, Y, Z, col = "darkgreen", size = 8)
  X <- v2[1]
  Y <- v2[2]
  Z <- v2[3]
  rgl::points3d(X, Y, Z, col = "darkgreen", size = 8)
  X <- v3[1]
  Y <- v3[2]
  Z <- v3[3]
  rgl::points3d(X, Y, Z, col = "darkgreen", size = 8)
  X <- v4[1]
  Y <- v4[2]
  Z <- v4[3]
  rgl::points3d(X, Y, Z, col = "darkgreen", size = 8)
  X <- vc[1]
  Y <- vc[2]
  Z <- vc[3]
  rgl::points3d(X, Y, Z, col = "darkorange", size = 10)
  X <- vr[1]
  Y <- vr[2]
  Z <- vr[3]
  rgl::points3d(X, Y, Z, col = "darkred", size = 10)
  rgl::triangles3d(c(v1[1], v2[1], v3[1]), 
                   c(v1[2], v2[2], v3[2]),
                   c(v1[3], v2[3], v3[3]),
                   alpha = 0.2, col = "darkgreen")
  rgl::triangles3d(c(v1[1], v2[1], v4[1]), 
                   c(v1[2], v2[2], v4[2]),
                   c(v1[3], v2[3], v4[3]),
                   alpha = 0.2, col = "darkgreen")
  rgl::triangles3d(c(v1[1], v3[1], v4[1]), 
                   c(v1[2], v3[2], v4[2]),
                   c(v1[3], v3[3], v4[3]),
                   alpha = 0.2, col = "darkgreen")
  rgl::triangles3d(c(v3[1], v2[1], v4[1]), 
                   c(v3[2], v2[2], v4[2]),
                   c(v3[3], v2[3], v4[3]),
                   alpha = 0.2, col = "darkgreen")
  X <- t1[1]-.07  # offset species labels position
  Y <- t1[2]-.07
  Z <- t1[3]
  rgl::text3d(X, Y, Z, text = sp_name[1])
  X <- t2[1]+.07
  Y <- t2[2]-.07
  Z <- t2[3]
  rgl::text3d(X, Y, Z, text = sp_name[2])
  X <- t3[1]
  Y <- t3[2]+.1
  Z <- t3[3]
  rgl::text3d(X, Y, Z, text = sp_name[3])
  X <- t4[1]
  Y <- t4[2]
  Z <- t4[3]+.1
  rgl::text3d(X, Y, Z, text = sp_name[4])
  rgl::aspect3d("iso")
}

# test_feasibility_trips ====
# test feasibility of species triplets
test_feasibility_trips <- function(alpha, r) {
  out <- c()
  v <- 1:4
  for (i in v) {
    if (test_feasibility(alpha[-i, -i], r[-i])) {
      out <- paste(out, paste(v[!v == i], collapse = " + "), sep = ", ")
    }
  }
  return(substr(out, 3, nchar(out)))
}

#################################################-
# Shiny server logic ----
#################################################-
# This section takes inputs from the ui and displays the desired
# information using the code from Saavedra et al. above
shinyServer(function(input, output, session) {
  # toggle to near-neutral scenario ####
  observeEvent(input$neutral, {
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(0.999, 0.998, 0.998,
                                       0.998, 0.999, 0.998,
                                       0.998, 0.998, 0.999),
                                     nrow = 3,
                                     dimnames = list(c("α(1,_)", "α(2,_)", "α(3,_)"), c("α(_,1)", "α(_,2)", "α(_,3)")),
                                     byrow = TRUE))
    updateNumericInput(session, "r1", value = 1)
    updateNumericInput(session, "r2", value = 1)
    updateNumericInput(session, "r3", value = 1)
  })
  # toggle to intransient (rock-paper-scissors) scenario ####
  observeEvent(input$intransient, {
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(0.5, 0.05, 1.5,
                                       1.5, 0.5, 0.05,
                                       0.05, 1.5, 0.5),
                                     nrow = 3,
                                     dimnames = list(c("α(1,_)", "α(2,_)", "α(3,_)"), c("α(_,1)", "α(_,2)", "α(_,3)")),
                                     byrow = TRUE))
    updateNumericInput(session, "r1", value = 1)
    updateNumericInput(session, "r2", value = 0.6)
    updateNumericInput(session, "r3", value = 0.9)
  })
  # toggle to weak intraspecific interactions scenario ####
  observeEvent(input$weakintra, {
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(1, 0.05, 0.05,
                                       0.05, 1, 0.05,
                                       0.05, 0.05, 1),
                                     nrow = 3,
                                     dimnames = list(c("α(1,_)", "α(2,_)", "α(3,_)"), c("α(_,1)", "α(_,2)", "α(_,3)")),
                                     byrow = TRUE))
    updateNumericInput(session, "r1", value = 1)
    updateNumericInput(session, "r2", value = 0.15)
    updateNumericInput(session, "r3", value = 0.15)
  })
  # render table with coexistence metrics -- 3 species ####
  output$stats <- renderTable({
    r <- c(input$r1, input$r2, input$r3)
    feas_pairs <- test_feasibility_pairs(input$alphamat, r)
    feas_txt <- ifelse(sum(feas_pairs$feasibility) == 0L,
                       "none",
                       ifelse(sum(feas_pairs$feasibility) == 1L,
                              paste(feas_pairs$pairs[, which(feas_pairs$feasibility == 1L)],
                                    collapse = " + "),
                              paste(apply(feas_pairs$pairs[, which(feas_pairs$feasibility == 1L)],
                                          MARGIN = 2, FUN = paste,
                                          collapse = " + "),
                                    collapse = ", ")))
    data.frame(Metric = c("Structural niche difference (Ω)",
                          "Structural fitness difference (θ)",
                          HTML("Centroid of feasibility domain (r<sub>c</sub>)"),
                          "Feasible triplet?",
                          paste0("Feasible pairs (",
                                 ifelse(feas_txt == "none",
                                        0,
                                        (nchar(feas_txt)+2)/7),
                                 "/3)")),
               Value = c(paste0(round(10^(Omega(alpha=input$alphamat)),
                                      3), "㏛"),
                         paste0(round(theta(input$alphamat, r), 3), "°"),
                         paste(round(r_centroid(input$alphamat), 3),
                               collapse = ", "),
                         ifelse(test_feasibility(input$alphamat, r) == 1L,
                                "yes", "no"),
                         feas_txt))},
    sanitize.text.function = function(x) x)
  # render table with coexistence metrics -- 4 species ####
  output$stats4sp <- renderTable({
    rr <- c(input$rr1, input$rr2, input$rr3, input$rr4)
    feas_pairs <- test_feasibility_pairs(input$alphamat4, rr)
    feas_trips <- test_feasibility_trips(input$alphamat4, rr)
    feas_txt <- ifelse(sum(feas_pairs$feasibility) == 0L,
                       "none",
                       ifelse(sum(feas_pairs$feasibility) == 1L,
                              paste(feas_pairs$pairs[, which(feas_pairs$feasibility == 1L)],
                                    collapse = " + "),
                              paste(apply(feas_pairs$pairs[, which(feas_pairs$feasibility == 1L)],
                                          MARGIN = 2, FUN = paste,
                                          collapse = " + "),
                                    collapse = ", ")))
    data.frame(Metric = c("Structural niche difference (Ω)",
                          "Structural fitness difference (θ)",
                          HTML("Centroid of feasibility domain (r<sub>c</sub>)"),
                          "Feasible quadruplet?",
                          paste0("Feasible triplets (", 
                                 (nchar(feas_trips)+2)/11,
                                 "/4)"),
                          paste0("Feasible pairs (",
                                 (nchar(feas_txt)+2)/7,
                                 "/6)")),
               Value = c(paste0(round(10^(Omega(alpha = input$alphamat4)),
                                      3), "㏛"),
                         paste0(round(theta(input$alphamat4, rr), 3),
                                "°"),
                         paste(round(r_centroid(input$alphamat4), 3),
                               collapse = ", "),
                         ifelse(test_feasibility(input$alphamat4, rr) == 1L,
                                "yes", "no"),
                         feas_trips,
                         feas_txt))},
    sanitize.text.function = function(x) x)
  # render network plot -- 3 species ####
  output$network_3sp <- renderPlot({
    mk_graph_3sp(input$alphamat,
                 rs = c(input$r1, input$r2, input$r3))
  })
  # render network plot -- 4 species ####
  output$network_4sp <- renderPlot({
    mk_graph_4sp(input$alphamat4,
                 rs = c(input$rr1, input$rr2, input$rr3, input$rr4))
  })
  # render static 3D cone ####
  output$cone <- renderPlot({
    cone_3sp(input$alphamat)
  })
  # render 3 species plot ####
  output$proj <- renderPlot({
    projection_3sp_with_pairwise(input$alphamat,
                                 r = c(input$r1, input$r2, input$r3))
  })
  # render interactive cone -- 3 species ####
  output$cone3d <- rgl::renderRglwidget({
    rgl::open3d(useNULL = TRUE)
    plot_cone_3D(input$alphamat, r = c(input$r1, input$r2, input$r3))
    rgl::rglwidget()
  })
  # render interactive cone -- 4 species ####
  output$proj4sp <- rgl::renderRglwidget({
    rgl::open3d(useNULL = TRUE)
    projection_4sp_3D(input$alphamat4, r = c(input$rr1, input$rr2,
                                             input$rr3, input$rr4))
    rgl::rglwidget()
  })
})
