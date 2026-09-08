library(Rcpp); library(RcppArmadillo); library(coda)
library(lattice); library(sna); library(stats)
library(mvtnorm); library(MASS)

wkdir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),   # RStudio: folder of this script
  error = function(e) getwd()                            # fallback: current working dir
)
setwd(wkdir)

simulation.1 <- function(n = 20) {
  beta  <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
  p     <- length(beta)
  sigma <- 3
  Cov   <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
  X     <- mvrnorm(n, rep(0, p), Cov)
  X     <- scale(X, TRUE, TRUE)
  y     <- X %*% beta + sigma * rnorm(n)
  y     <- scale(y, TRUE, FALSE)
  list(y = y, X = X, V = Cov, beta = beta)
}

mydata <- simulation.1(20)
Y <- as.vector(mydata$y); Y <- Y - mean(Y)
X <- scale(mydata$X)              # same as the kronecker thing, cleaner
n.data <- length(Y)
p.data <- ncol(X)

Beta_init    <- rnorm(p.data, 0, 0)
sigma2_init  <- var(Y)
Lambda_init  <- diag(1, p.data)
a_init <- 0.1; b_init <- 0.1; r_init <- 0.1
ahyper <- 0.01; bhyper <- 0.01; rhyper <- 0.01
Nburn  <- 10000; Niter <- 20000; Nchain <- 1; Nthin <- 1

source("BVSG.r")                  # this in turn sourceCpp's BVSG.cpp

MCMCResult <- BVSGR(X, Y, a_init, b_init, r_init,
                    Beta_init, Lambda_init, sigma2_init,
                    ahyper, bhyper, rhyper,
                    Nburn, Niter, Nchain, Nthin)   # <-- Nthin added

betahat <- MCMCResult[[1]]$Betahat
print(betahat)