library(tidyverse)
library(grf)
library(nloptr)
library(ggplot2)


load("Diabetes.RData")


compute_ps <- function(S, X, A,
                       Xtest, Stest) {
  A <- as.factor(A)
  X <- cbind(X, S)
  Xeval <- cbind(Xtest, Stest)
  
  grfA <- probability_forest(X, A, 
                             honesty = FALSE, num.trees = 2000,
                             sample.fraction = 0.8, ci.group.size = 1)
  out <- list()
  out$ps <- c(grfA$predictions[, 2])
  out$ps.eval <- c(predict(grfA, Xeval)$predictions[, 2])
  return(out)
}

compute_or <- function(Y, S, X, A) {
  X <- cbind(X, S)
  grfY <- regression_forest(X = cbind(X, A), Y = Y)
  mu0 <- predict(grfY, cbind(X, 0))$predictions
  mu1 <- predict(grfY, cbind(X, 1))$predictions
  out <- list(mu0 = mu0,
              mu1 = mu1)
  return(out)
}


#########################
## Diabetes data analysis
#########################

X <- DATA %>% select(race, gender, age,
                     time_in_hospital, num_lab_procedures,
                     num_medications, number_diagnoses)

X$time_in_hospital <- X$time_in_hospital / max(X$time_in_hospital)
X$num_lab_procedures <- X$num_lab_procedures / max(X$num_lab_procedures)
X$num_medications <- X$num_medications / max(X$num_medications)
X$number_diagnoses <- X$number_diagnoses / max(X$number_diagnoses)

A <- DATA$diabetesMed
Y <- DATA$readmit

## check propensity score
grf.A <- probability_forest(X = X, Y = as.factor(A))
plot(ecdf(grf.A$predictions[, 2]))


err <- rnorm(nrow(X))
Y0 <- 20 * (1 + X$gender - X$age + X$time_in_hospital + X$num_lab_procedures + X$num_medications + X$num_medications^2 + exp(X$number_diagnoses)) + err * 20
Y1 <- 20 * (1 + X$gender - X$age + X$time_in_hospital + X$num_lab_procedures + X$num_medications + X$num_medications^2 + exp(X$number_diagnoses)) + 
  25 * (3 - 5 * X$age + 2 * X$time_in_hospital - 3 * X$num_medications + X$race) + 
  err * 20

Y <- A * Y1 + (1 - A) * Y0

data <- cbind(X, A, Y, Y0, Y1)


num.MonteCarlo <- 50
value <- matrix(nrow = num.MonteCarlo, ncol = 6)

n.total <- 2500
n.train <- 500
n.test <- n.total - n.train
tau <- 0.01

set.seed(156)
for (i in 1:num.MonteCarlo) {
  cat(i, as.character(Sys.time()), "\n")
  
  data.train.test <- data[sample(nrow(X), n.total), ]
  data.train <- data.train.test[1:n.train, ]
  data.test <- data.train.test[(n.train + 1):n.total, ]
  
  s <- data.train$race
  x <- as.matrix(data.train[, c("gender", "age", "time_in_hospital", "num_lab_procedures", "num_medications", "number_diagnoses")])
  a <- data.train$A
  y <- data.train$Y
  
  s.test <- data.test$race
  x.test <- as.matrix(data.test[, c("gender", "age", "time_in_hospital", "num_lab_procedures", "num_medications", "number_diagnoses")])
  y0 <- data.test$Y0
  y1 <- data.test$Y1
  
  # standard methods
  prob <- compute_ps(s, x, a, x.test, s.test)$ps
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
  
  test.cobyla <- try(do.call(cobyla, list(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_est, hin = dp_est,
                                          nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))))
  if(inherits(test.cobyla, "try-error")) {
    value[i, 1] <- NA
  } else {
    eta <- test.cobyla$par
    d.eval <- as.numeric(c(cbind(1, s.test, x.test) %*% eta) > 0)
    value[i, 1] <- mean(d.eval * y1 + (1 - d.eval) * y0)
  }
  
  # OR
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    -mean(d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  test.cobyla <- try(do.call(cobyla, list(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_est, hin = dp_est,
                                          nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))))
  if(inherits(test.cobyla, "try-error")) {
    value[i, 2] <- NA
  } else {
    eta <- test.cobyla$par
    d.eval <- as.numeric(c(cbind(1, s.test, x.test) %*% eta) > 0)
    value[i, 2] <- mean(d.eval * y1 + (1 - d.eval) * y0)
  }
  
  # AIPW
  mean_est <- function(eta) {
    d.x <- as.numeric(c(cbind(1, s, x) %*% eta) > 0)
    -mean(as.numeric(a == d.x) * (y - a * mu$mu1 - (1 - a) * mu$mu0) / (a * prob + (1 - a) * (1 - prob)) + d.x * mu$mu1 + (1 - d.x) * mu$mu0)
  }
  
  test.cobyla <- try(do.call(cobyla, list(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_est, hin = dp_est,
                                          nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))))
  if(inherits(test.cobyla, "try-error")) {
    value[i, 3] <- NA
  } else {
    eta <- test.cobyla$par
    d.eval <- as.numeric(c(cbind(1, s.test, x.test) %*% eta) > 0)
    value[i, 3] <- mean(d.eval * y1 + (1 - d.eval) * y0)
  }
  
  
  # proposed methods
  prob <- compute_ps(s, x, a, x.test, s.test)
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
  
  res <- cobyla(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, s.test, x.test) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[i, 4] <- mean(d.eval * y1 + (1 - d.eval) * y0)
  
  # OR-IPS
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    -mean((delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob))
  }
  
  res <- cobyla(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_est, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, s.test, x.test) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[i, 5] <- mean(d.eval * y1 + (1 - d.eval) * y0)
  
  # One-step
  mean_est <- function(eta) {
    delta.x <- c(exp(cbind(1, s, x) %*% eta))
    -mean((a * delta.x * (y - mu$mu1) + (1 - a) * (y - mu$mu0) + delta.x * prob * mu$mu1 + (1 - prob) * mu$mu0) / (delta.x * prob + 1 - prob) 
          + (delta.x * (mu$mu1 - mu$mu0) * (a - prob)) / (delta.x * prob + 1 - prob)^2)
  }
  
  res <- cobyla(x0 = runif(ncol(x) + 2, -1, 1), fn = mean_esto, hin = dp_est,
                nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))
  eta <- res$par
  delta.x.eval <- c(exp(cbind(1, s.test, x.test) %*% eta))
  d.eval <- delta.x.eval * ps.eval / (delta.x.eval * ps.eval + 1 - ps.eval)
  value[i, 6] <- mean(d.eval * y1 + (1 - d.eval) * y0)
}



pdf("Diabetes_data.pdf",
    height = 6, width = 8)
boxplot(value[, 1], value[, 2], value[, 3], value[, 4], value[, 5], value[, 6],
        ylab="True value",
        names = c("IPW", "OR", "AIPW", "IPW-IPS", "OR-IPS", "One-step"))
dev.off()


