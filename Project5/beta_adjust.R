# Project 5, parts a-c: optimisation under the single index model, and the
# Blume and Vasicek beta adjustments.

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "beta_adjust.R")
    setwd(dirname(normalizePath(f)))
})

DATA <- ".."

# ---- data -------------------------------------------------------------------
# Period 1 is the training file, the same 59 returns Projects 1, 2 and 4 use, so
# the betas here are Project 4's betas. Period 2 is the test file.

read_ret <- function(f) {
  a <- read.csv(file.path(DATA, f), sep = ",", header = TRUE, check.names = FALSE)
  stopifnot(!any(duplicated(names(a))))
  P <- as.matrix(a[, -1][, -1])
  stopifnot(ncol(P) == 31, !anyNA(P))
  n <- nrow(P)
  P[-1, ] / P[-n, ] - 1
}

ret1 <- read_ret("stockData_train.csv")     # 01-Jan-2017 to 01-Jan-2022
ret2 <- read_ret("stockData_test.csv")      # 01-Jan-2022 to 31-Jul-2026
stopifnot(nrow(ret1) == 59, nrow(ret2) == 54)

stk <- setdiff(colnames(ret1), "^GSPC")
stopifnot(length(stk) == 30, identical(stk, setdiff(colnames(ret2), "^GSPC")))

# ---- single index fits in each period ---------------------------------------

fit_period <- function(r) {
  Rm  <- as.numeric(r[, "^GSPC"])
  r30 <- r[, stk]
  fit <- lapply(stk, function(s) lm(r30[, s] ~ Rm))
  names(fit) <- stk
  list(
    Rm    = Rm,
    Rbar  = colMeans(r30),
    alpha = vapply(fit, function(f) unname(coef(f)[1]), numeric(1)),
    beta  = vapply(fit, function(f) unname(coef(f)[2]), numeric(1)),
    se    = vapply(fit, function(f) unname(summary(f)$coefficients[2, 2]), numeric(1)),
    sig2e = vapply(fit, function(f) summary(f)$sigma^2, numeric(1)),
    sd    = apply(r30, 2, sd),
    Tn    = nrow(r30)
  )
}

p1 <- fit_period(ret1)
p2 <- fit_period(ret2)

stopifnot(max(abs(p1$beta - apply(ret1[, stk], 2,
                  function(x) cov(x, p1$Rm) / var(p1$Rm)))) < 1e-12)
stopifnot(max(abs(p1$alpha - (p1$Rbar - p1$beta * mean(p1$Rm)))) < 1e-15)

# ---- (a) optimal portfolio, Z = Sigma^-1 R ----------------------------------

Rf  <- 0.001
pos <- p1$beta > 0

cat("(a) stocks with positive beta:", sum(pos), "of", length(pos), "\n")

sel   <- stk[pos]
beta  <- p1$beta[pos]
sig2e <- p1$sig2e[pos]
Rbar  <- p1$Rbar[pos]
sig2m <- var(p1$Rm)

Sigma <- outer(beta, beta) * sig2m
diag(Sigma) <- beta^2 * sig2m + sig2e
stopifnot(isSymmetric(Sigma), min(eigen(Sigma, only.values = TRUE)$values) > 0)

Rex <- Rbar - Rf
Z   <- solve(Sigma, Rex)
x   <- Z / sum(Z)

# Sherman-Morrison: Sigma is diagonal plus rank one, so its inverse is explicit
D_inv <- diag(1 / sig2e)
sm    <- D_inv - sig2m * (D_inv %*% outer(beta, beta) %*% D_inv) /
                 as.numeric(1 + sig2m * t(beta) %*% D_inv %*% beta)
stopifnot(max(abs(sm %*% Rex - Z)) < 1e-9)

# The cut-off form of the same solution, handout #13
C_star <- sig2m * sum(Rex * beta / sig2e) / (1 + sig2m * sum(beta^2 / sig2e))
z_egp  <- (beta / sig2e) * (Rex / beta - C_star)
stopifnot(max(abs(z_egp - Z)) < 1e-9)

E_G      <- sum(x * Rbar)
sd_G     <- sqrt(as.numeric(t(x) %*% Sigma %*% x))
sharpe_G <- (E_G - Rf) / sd_G
beta_G   <- sum(x * beta)

# No reweighting beats it
set.seed(183)
stopifnot(sharpe_G >= max(replicate(4000, {
  w <- x + rnorm(length(x), 0, 0.02); w <- w / sum(w)
  (sum(w * Rbar) - Rf) / sqrt(as.numeric(t(w) %*% Sigma %*% w))
})))

cat(sprintf("    Rf = %.4f, C* = %.6f\n", Rf, C_star))
cat(sprintf("    E %.4f%%, sd %.4f%%, Sharpe %.4f, portfolio beta %.4f\n",
            E_G * 100, sd_G * 100, sharpe_G, beta_G))
cat(sprintf("    %d short positions, sum|x| = %.2f\n\n", sum(x < 0), sum(abs(x))))

# ---- (b) Blume ---------------------------------------------------------------
#   beta2 = a + b beta1, then push the period 2 betas through the fitted line

blume     <- lm(p2$beta ~ p1$beta)
b_int     <- unname(coef(blume)[1])
b_slope   <- unname(coef(blume)[2])
b_r2      <- summary(blume)$r.squared
beta_next <- b_int + b_slope * p2$beta

# OLS passes through the means
stopifnot(abs((b_int + b_slope * mean(p1$beta)) - mean(p2$beta)) < 1e-12)

cat(sprintf("(b) Blume: beta2 = %.4f + %.4f beta1   (R^2 = %.3f)\n",
            b_int, b_slope, b_r2))
cat(sprintf("    beta mean/sd:  p1 %.4f/%.4f   p2 %.4f/%.4f   forecast %.4f/%.4f\n\n",
            mean(p1$beta), sd(p1$beta), mean(p2$beta), sd(p2$beta),
            mean(beta_next), sd(beta_next)))

# ---- (b) Vasicek -------------------------------------------------------------
#   beta_adj = (var_bi/(var_bi+var_b1)) beta_bar + (var_b1/(var_bi+var_b1)) beta_i

var_bi   <- p1$se^2                      # squared standard error of each beta
var_b1   <- var(p1$beta)                 # cross-sectional variance of the betas
beta_bar <- mean(p1$beta)
w_prior  <- var_bi / (var_bi + var_b1)
beta_vas <- w_prior * beta_bar + (1 - w_prior) * p1$beta

stopifnot(all(w_prior > 0 & w_prior < 1))
stopifnot(all(abs(beta_vas - beta_bar) <= abs(p1$beta - beta_bar) + 1e-12))

cat("(b) Vasicek\n")
cat(sprintf("    beta_bar %.4f, cross-sectional variance %.5f\n", beta_bar, var_b1))
cat(sprintf("    weight on the mean: %.3f to %.3f (median %.3f)\n",
            min(w_prior), max(w_prior), median(w_prior)))
cat(sprintf("    beta sd: raw %.4f -> adjusted %.4f\n\n", sd(p1$beta), sd(beta_vas)))

# ---- (c) PRESS ---------------------------------------------------------------

# Sum of squared forecast errors, the Klemkosky and Martin (1975) form used in
# Homework 4 exercise 1.
press_vas <- sum((beta_vas - p2$beta)^2)
press_raw <- sum((p1$beta  - p2$beta)^2)

# The unadjusted figure is the PRESS Homework 4 computes, so it must agree
stopifnot(abs(press_raw - 3.778242) < 1e-5)

cat("(c) PRESS against the realised period 2 betas\n")
cat(sprintf("    Vasicek       %.4f   (mean square %.6f)\n",
            press_vas, press_vas / length(stk)))
cat(sprintf("    unadjusted    %.4f   (mean square %.6f)\n\n",
            press_raw, press_raw / length(stk)))

betas <- data.frame(stock = stk, beta1 = p1$beta, beta2 = p2$beta,
                    vasicek = beta_vas, blume_next = beta_next, row.names = NULL)

# ---- plotting helpers -------------------------------------------------------

surface <- "#fcfcfb"; ink <- "#0b0b0b"; muted <- "#898781"; grid <- "#e1e0d9"
c_front <- "#2a78d6"; c_eq <- "#eb6834"; c_spx <- "#1baf7a"; c_arb <- "#6f6d67"

open_png <- function(f, h = 4.7) {
  png(f, width = 6.5, height = h, units = "in", res = 200)
  par(bg = surface, mar = c(4.2, 4.6, 3.4, 1.4), family = "sans")
}
ax <- function(side, at, lab) {
  axis(side, at = at, labels = lab, col = "#c3c2b7", col.axis = muted,
       cex.axis = 0.75, lwd = 1, las = if (side == 2) 1 else 0)
}
pt <- function(x, y, col, lab, pos = 4, cex = 0.72, off = 0.6) {
  points(x, y, pch = 21, bg = col, col = surface, lwd = 1.8, cex = 1.5)
  text(x, y, lab, pos = pos, offset = off, cex = cex, font = 2, col = col)
}

# ---- (a) the optimal portfolio and the capital allocation line --------------

open_png("fig_tangency.png", h = 5.0)
xr <- c(-0.012, max(p1$sd) * 1.04)
yr <- c(-0.004, max(c(p1$Rbar, E_G)) * 1.14)
xtick <- pretty(c(0, max(p1$sd)))
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = xtick, col = grid, lwd = 1)

cal_x <- sd_G * 1.15
segments(0, Rf, cal_x, Rf + sharpe_G * cal_x, col = c_eq, lwd = 2)
text(cal_x, Rf + sharpe_G * cal_x, "CAL", pos = 4, offset = 0.3, cex = 0.75,
     font = 2, col = c_eq)

points(p1$sd, p1$Rbar, pch = 21, bg = "#d9d8d2", col = surface, lwd = 1, cex = 0.85)
pt(sd(p1$Rm), mean(p1$Rm), c_spx, "S&P 500", 1)
pt(0, Rf, c_eq, expression(R[f]), 2, off = 0.45)
pt(sd_G, E_G, c_front, "optimal portfolio", 2)

ax(1, xtick, sprintf("%.0f%%", xtick * 100))
ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(a) Optimal portfolio under the single index model", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.05, font = 2)
mtext("Grey points are the 30 stocks over period 1, short sales allowed.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

# ---- (b) Blume ---------------------------------------------------------------

open_png("fig_blume.png")
rr <- range(c(p1$beta, p2$beta)) + c(-0.12, 0.12)
plot(NA, xlim = rr, ylim = rr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(rr), v = pretty(rr), col = grid, lwd = 1)
abline(0, 1, col = "#c9c8c0", lwd = 1, lty = 2)
text(rr[2], rr[2], "no change", pos = 2, offset = 0.3, cex = 0.7, col = "#b4b3ab")
abline(b_int, b_slope, col = c_eq, lwd = 2.4)
points(p1$beta, p2$beta, pch = 21, bg = c_arb, col = surface, lwd = 1.2, cex = 1.05)

ax(1, pretty(rr), sprintf("%.1f", pretty(rr)))
ax(2, pretty(rr), sprintf("%.1f", pretty(rr)))
mtext("Beta, period 1", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Beta, period 2", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(b) Blume regression", 3, line = 1.8, adj = 0, col = ink, cex = 1.05, font = 2)
mtext(sprintf("Fitted line beta2 = %.4f + %.4f beta1, against the 45 degree line.",
              b_int, b_slope), 3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

# ---- (c) Vasicek forecasts against the realised betas ----------------------

open_png("fig_vasicek.png")
rr <- range(c(p1$beta, p2$beta, beta_vas)) + c(-0.12, 0.12)
plot(NA, xlim = rr, ylim = rr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(rr), v = pretty(rr), col = grid, lwd = 1)
abline(0, 1, col = "#c9c8c0", lwd = 1.4, lty = 2)
text(rr[2], rr[2], "perfect forecast", pos = 2, offset = 0.3, cex = 0.7, col = "#b4b3ab")
segments(p1$beta, p2$beta, beta_vas, p2$beta, col = "#dcdbd4", lwd = 1)
points(p1$beta, p2$beta, pch = 21, bg = c_arb, col = surface, lwd = 1.1, cex = 0.95)
points(beta_vas, p2$beta, pch = 21, bg = c_front, col = surface, lwd = 1.1, cex = 0.95)

legend("topleft", bty = "n", cex = 0.72, text.col = muted,
       pch = 21, pt.bg = c(c_arb, c_front), col = surface, pt.lwd = 1.1,
       legend = c(sprintf("unadjusted   (PRESS %.3f)", press_raw),
                  sprintf("Vasicek      (PRESS %.3f)", press_vas)))

ax(1, pretty(rr), sprintf("%.1f", pretty(rr)))
ax(2, pretty(rr), sprintf("%.1f", pretty(rr)))
mtext("Forecast beta", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Beta, period 2", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(c) Vasicek forecasts against the realised betas", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.05, font = 2)
mtext("Each grey segment joins a stock's raw period 1 beta to its adjusted value.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

cat("wrote fig_tangency.png, fig_blume.png, fig_vasicek.png\n")
