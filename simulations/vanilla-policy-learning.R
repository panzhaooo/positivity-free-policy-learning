library(grf)
library(rgenoud)
library(MASS)


simData <- function(n = 1000) {
  X <- matrix(data = c(runif(n * 2),
                       mvrnorm(n, mu = rep(0, 2), Sigma = matrix(c(1, 0.3, 0.3, 1), 2, 2))),
              ncol = 4)
  
  propensity_score <- plogis(0.35 - 0.3 * (X[, 1] + X[, 2] + X[, 3] + X[, 4]))
  A <- as.numeric(runif(n) < propensity_score)
  
  err <- rnorm(n)
  Y <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + 
    A * 25 * (3 - 5 * X[, 1] + 2 * X[, 2] - 3 * X[, 3] + X[, 4]) + 
    err * 20
  
  out <- list(X = X, A = A, Y = Y)
  return(out)
}

compute_ps <- function(X, A) {
  A <- as.factor(A)
  grfA <- probability_forest(X, A, 
                             honesty = FALSE, num.trees = 2000,
                             sample.fraction = 0.8, ci.group.size = 1)
  out <- list()
  out$ps <- c(grfA$predictions[, 2])
  out$ps.eval <- c(predict(grfA, Xeval)$predictions[, 2])
  return(out)
}

compute_or <- function(Y, X, A) {
  grfY <- regression_forest(X = cbind(X, A), Y = Y)
  mu0 <- predict(grfY, cbind(X, 0))$predictions
  mu1 <- predict(grfY, cbind(X, 1))$predictions
  out <- list(mu0 = mu0,
              mu1 = mu1)
  return(out)
}

ITR_estimators <- function(data) {
  x <- data$X
  a <- data$A
  y <- data$Y
  
  value <- numeric(6)
  
  # standard methods
  prob <- compute_ps(x, a)$ps
  mu <- compute_or(y, x, a)
  
  # IPW
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, x) %*% eta) > 0)
    mean(y * as.numeric(a == d.x) / (a * prob + (1 - a) * (1 - prob)))
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  d.eval <- as.numeric(c(cbind(1, Xeval) %*% eta) > 0)
  value[1] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # OR
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, x) %*% eta) > 0)
    mean(d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  d.eval <- as.numeric(c(cbind(1, Xeval) %*% eta) > 0)
  value[2] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # AIPW
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, x) %*% eta) > 0)
    mean(as.numeric(a == d.x) * (y - a * mu$mu1 - (1 - a) * mu$mu0) / (a * prob + (1 - a) * (1 - prob)) + d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  d.eval <- as.numeric(c(cbind(1, Xeval) %*% eta) > 0)
  value[3] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # proposed methods
  prob <- compute_ps(x, a)
  ps.eval <- prob$ps.eval
  prob <- prob$ps
  mu <- compute_or(y, x, a)
  
  # IPW-IPS
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, x) %*% eta))
    mean(y * (delta.x * a + 1 - a) / (delta.x * prob + 1 - prob))
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  delta.x.eval <- c(exp(cbind(1, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  d.eval[is.na(d.eval)] <- 1
  value[4] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # OR-IPS
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, x) %*% eta))
    mean((delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob))
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  delta.x.eval <- c(exp(cbind(1, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  d.eval[is.na(d.eval)] <- 1
  value[5] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # One-step
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, x) %*% eta))
    mean((a * delta.x * (y - mu$mu1) + (1 - a) * (y - mu$mu0) + delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob) 
         + (delta.x * (mu$mu1 - mu$mu0) * (a - prob)) / (delta.x * prob + 1 - prob)^2)
  }
  
  Vopt <- genoud(mean_est, default.domains = 1000, nvars = ncol(x) + 1,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  eta <- Vopt$par
  delta.x.eval <- c(exp(cbind(1, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  d.eval[is.na(d.eval)] <- 1
  value[6] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  out <- list(value = value)
  return(out)
}

compute.true.value <- function(n = 1e6) {
  X <- matrix(data = c(runif(n * 2),
                       mvrnorm(n, mu = rep(0, 2), Sigma = matrix(c(1, 0.3, 0.3, 1), 2, 2))),
              ncol = 4)
  
  err <- rnorm(n)
  Y0 <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + err * 20
  Y1 <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + 
    25 * (3 - 5 * X[, 1] + 2 * X[, 2] - 3 * X[, 3] + X[, 4]) + err * 20
  
  out <- mean(pmax(Y0, Y1))
  return(out)
}

##########################
## vanilla policy learning

set.seed(1752)

N <- 1e5
Xeval <- matrix(data = c(runif(N * 2),
                         mvrnorm(N, mu = rep(0, 2), Sigma = matrix(c(1, 0.3, 0.3, 1), 2, 2))),
                ncol = 4)
err <- rnorm(N)
Y0 <- 20 * (1 + Xeval[, 1] - Xeval[, 2] + Xeval[, 3]^2 + exp(Xeval[, 2])) + err * 20
Y1 <- 20 * (1 + Xeval[, 1] - Xeval[, 2] + Xeval[, 3]^2 + exp(Xeval[, 2])) + 
  25 * (3 - 5 * Xeval[, 1] + 2 * Xeval[, 2] - 3 * Xeval[, 3] + Xeval[, 4]) + err * 20

mean(pmax(Y0, Y1))
mean(Y1 > Y0)

opt.value <- compute.true.value()


num_mc <- 100
value <- matrix(0, nrow = num_mc, ncol = 6)

for (i in 1:num_mc) {
  cat(i, as.character(Sys.time()), "\n")
  
  data <- simData(n = 2000)
  res <- ITR_estimators(data)
  value[i, ] <- res$value
}



pdf("vanilla-policy-learning.pdf",
    height = 6, width = 8)
boxplot(value[, 1], value[, 2], value[, 3], value[, 4], value[, 5], value[, 6],
        ylab="True value", ylim = c(125, 130),
        names = c("IPW", "OR", "AIPW", "IPW-IPS", "OR-IPS", "One-step"))
abline(h = opt.value, col="blue")
dev.off()


