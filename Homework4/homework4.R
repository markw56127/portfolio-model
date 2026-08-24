# Homework 4, Exercise 1
# PRESS decomposition, Klemkosky and Martin (1975), equation (2).
# Period 1 betas are used as predictions of the period 2 betas
# for the 30 stocks in the project data.
#
# Run with the Homework4 folder as the working directory.
# Uses base R only.

# Estimate each stock's beta by regressing its return on the S&P 500.
betas <- function(file) {
  a <- read.csv(file, header = TRUE, check.names = FALSE)
  p <- as.matrix(a[, -c(1, 2)])        # drop the row index and Date
  n <- nrow(p)
  r <- p[-1, ] / p[-n, ] - 1           # simple monthly returns
  Rm <- r[, "^GSPC"]
  stocks <- setdiff(colnames(r), "^GSPC")
  sapply(stocks, function(s) unname(coef(lm(r[, s] ~ Rm))[2]))
}

P <- betas("../stockData_train.csv")   # period 1, the predictions
A <- betas("../stockData_test.csv")    # period 2, what happened
n <- length(P)

# Regress the actual betas on the predicted ones.
fit  <- lm(A ~ P)
b    <- unname(coef(fit)[2])
Ahat <- fitted(fit)

# The three components of equation (2).
bias   <- n * (mean(A) - mean(P))^2
ineff  <- (b - 1)^2 * sum((P - mean(P))^2)
random <- sum((A - Ahat)^2)

PRESS <- sum((A - P)^2)

# PRESS and its decomposition.
round(c(bias = bias, inefficiency = ineff, random = random,
        sum = bias + ineff + random, PRESS = PRESS), 6)
