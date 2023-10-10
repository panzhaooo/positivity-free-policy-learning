library(nloptr)


simData <- function(n = 1000) {
  S <- rbinom(n, 1, 0.5)
  X <- matrix(data = runif(n * 3), ncol = 3)
  
  propensity_score <- plogis(-1 - X[, 1] + 1.5 * X[, 2] - 0.25 * X[, 3] - 3.1 * S)
  A <- as.numeric(runif(n) < propensity_score)
  
  err <- rnorm(n)
  Y <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + 
    A * 25 * (3 - 5 * X[, 1] + 2 * X[, 2] - 3 * X[, 3] + S) + 
    err * 20
  
  out <- list(S = S, X = X, A = A, Y = Y)
  return(out)
}

compute_ps <- function(S, X, A) {
  psData <- data.frame(A = A, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3], S = S)
  glmA <- glm(A ~ X1 + X2 + X3 + S, 
              family = binomial(link = "logit"),
              data = psData)
  
  evalData <- data.frame(X1 = Xeval[, 1], X2 = Xeval[, 2], X3 = Xeval[, 3], S = Seval)
  out <- list()
  out$ps <- glmA$fitted.values
  out$ps.eval <- predict(glmA, evalData, type = "response")
  return(out)
}

compute_or <- function(Y, S, X, A) {
  orData <- data.frame(Y = Y, A = A, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3], S = S)
  glmY <- glm(Y ~ X1 + X2 + X3^2 + exp(X2) +
                A + A:X1 + A:X2 + A:X3 + A:S, 
              family = gaussian(),
              data = orData)
  
  orData0 <- data.frame(A = 0, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3], S = S)
  orData1 <- data.frame(A = 1, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3], S = S)
  mu0 <- predict(glmY, orData0)
  mu1 <- predict(glmY, orData1)
  out <- list(mu0 = mu0,
              mu1 = mu1)
  return(out)
}

ITR_estimators <- function(data, tau = 0.01) {
  s <- data$S
  x <- data$X
  a <- data$A
  y <- data$Y
  
  value <- numeric(6)
  
  # standard methods
  prob <- compute_ps(s, x, a)$ps
  mu <- compute_or(y, s, x, a)
  
  dp_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    tau - sqrt((mean(d.x[s == 1]) - mean(d.x))^2 + (mean(d.x[s == 0]) - mean(d.x))^2)
  }
  
  # IPW
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    mean(-y * as.numeric(a == d.x) / (a * prob + (1 - a) * (1 - prob)))
  }
  
  test.cobyla <- try(do.call(cobyla, list(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                                          nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))))
  if(inherits(test.cobyla, "try-error")) {
    value[1] <- NA
  } else {
    eta <- test.cobyla$par
    d.eval <- as.numeric(c(cbind(1, Seval, Xeval) %*% eta) > 0)
    value[1] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  }
  
  # OR
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    -mean(d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  res <- cobyla(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  d.eval <- as.numeric(c(cbind(1, Seval, Xeval) %*% eta) > 0)
  value[2] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # AIPW
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    -mean(as.numeric(a == d.x) * (y - a * mu$mu1 - (1 - a) * mu$mu0) / (a * prob + (1 - a) * (1 - prob)) + d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  test.cobyla <- try(do.call(cobyla, list(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                                          nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))))
  if(inherits(test.cobyla, "try-error")) {
    value[3] <- NA
  } else {
    eta <- test.cobyla$par
    d.eval <- as.numeric(c(cbind(1, Seval, Xeval) %*% eta) > 0)
    value[3] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  }
  
  
  # proposed methods
  prob <- compute_ps(s, x, a)
  ps.eval <- prob$ps.eval
  prob <- prob$ps
  mu <- compute_or(y, s, x, a)
  
  dp_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    d.x <- delta.x * prob / (delta.x * prob + 1 - prob)
    tau - sqrt((mean(d.x[s == 1]) - mean(d.x))^2 + (mean(d.x[s == 0]) - mean(d.x))^2)
  }
  
  # IPW-IPS
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    -mean(y * (delta.x * a + 1 - a) / (delta.x * prob + 1 - prob))
  }
  
  res <- cobyla(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, Seval, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[4] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # OR-IPS
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    -mean((delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob))
  }
  
  res <- cobyla(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, Seval, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[5] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  # One-step
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    -mean((a * delta.x * (y - mu$mu1) + (1 - a) * (y - mu$mu0) + delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob) 
          + (delta.x * (mu$mu1 - mu$mu0) * (a - prob)) / (delta.x * prob + 1 - prob)^2)
  }
  
  res <- cobyla(x0 = rnorm(ncol(x) + 2), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, Seval, Xeval) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[6] <- mean(d.eval * Y1 + (1 - d.eval) * Y0)
  
  out <- list(value = value)
  return(out)
}

compute.true.value <- function(n = 1e6,
                               tau = 0.01) {
  S <- rbinom(n, 1, 0.5)
  X <- matrix(data = runif(n * 3), ncol = 3)
  
  propensity_score <- plogis(-1 - X[, 1] + 1.5 * X[, 2] - 0.25 * X[, 3] - 3.1 * S)

  err <- rnorm(n)
  Y0 <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + err * 20
  Y1 <- 20 * (1 + X[, 1] - X[, 2] + X[, 3]^2 + exp(X[, 2])) + 
    25 * (3 - 5 * X[, 1] + 2 * X[, 2] - 3 * X[, 3] + S) + 
    err * 20
  
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, S, X) %*% eta))
    d.x <- delta.x * propensity_score / (delta.x * propensity_score + 1 - propensity_score)
    -mean(d.x * Y1 + (1 - d.x) * Y0)
  }
  
  dp_est <- function(eta) {
    delta.x <- c(exp(cbind(1, S, X) %*% eta))
    d.x <- delta.x * propensity_score / (delta.x * propensity_score + 1 - propensity_score)
    tau - sqrt((mean(d.x[S == 1]) - mean(d.x))^2 + (mean(d.x[S == 0]) - mean(d.x))^2)
  }
  
  res <- cobyla(x0 = rnorm(ncol(X) + 2), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  
  out <- -res$value
  return(out)
}


#######################
## fair policy learning
## demographic parity
## parametric models

set.seed(9920)

N <- 1e5
Seval <- rbinom(N, 1, 0.5)
Xeval <- matrix(data = runif(N * 3), ncol = 3)

err <- rnorm(N)
Y0 <- 20 * (1 + Xeval[, 1] - Xeval[, 2] + Xeval[, 3]^2 + exp(Xeval[, 2])) + err * 20
Y1 <- 20 * (1 + Xeval[, 1] - Xeval[, 2] + Xeval[, 3]^2 + exp(Xeval[, 2])) + 
  25 * (3 - 5 * Xeval[, 1] + 2 * Xeval[, 2] - 3 * Xeval[, 3] + Seval) + 
  err * 20 * (1 + 1 * (1 + Xeval[, 1] + Xeval[, 2] + Xeval[, 3] + Seval))

mean(pmax(Y0, Y1))
mean(Y1 > Y0)

opt.value <- compute.true.value()


num_mc <- 100
value <- matrix(0, nrow = num_mc, ncol = 6)

for (i in 1:num_mc) {
  cat(i, as.character(Sys.time()), "\n")
  
  data <- simData(n = 500)
  res <- ITR_estimators(data, tau = 0.01)
  value[i, ] <- res$value
}



pdf("fairness-demographic-parity-para.pdf",
    height = 6, width = 8)
boxplot(value[, 1], value[, 2], value[, 3], value[, 4], value[, 5], value[, 6],
        ylab="True value", ylim = c(55, 90),
        names = c("IPW", "OR", "AIPW", "IPW-IPS", "OR-IPS", "One-step"))
abline(h = opt.value, col="blue")
dev.off()


