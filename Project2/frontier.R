local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "frontier.R")
    setwd(dirname(normalizePath(f)))
})

DATA <- ".."

# ---- data -------------------------------------------------------------------

a <- read.csv(file.path(DATA, "stockData_train.csv"), sep = ",", header = TRUE,
              check.names = FALSE)
stopifnot(!any(duplicated(names(a))))

prices <- a[, -1]
dates  <- as.Date(prices$Date)
P      <- as.matrix(prices[, -1])
stopifnot(nrow(P) == 60, ncol(P) == 31, !anyNA(P))

n   <- nrow(P)
ret <- P[-1, ] / P[-n, ] - 1              # 59 x 31 simple monthly returns

stk    <- setdiff(colnames(ret), "^GSPC")
ret30  <- ret[, stk]
Rbar   <- colMeans(ret30)                 # 30 x 1 mean vector
S      <- cov(ret30)                      # 30 x 30 variance-covariance matrix
stdev  <- sqrt(diag(cov(ret)))            # all 31, for the plot
means  <- colMeans(ret)

stopifnot(ncol(S) == 30, qr(S)$rank == 30)

# ---- (a) A, B, C, D ---------------------------------------------------------
# A = 1' S^-1 Rbar,  B = Rbar' S^-1 Rbar,  C = 1' S^-1 1,  D = BC - A^2

Sinv <- solve(S)
one  <- rep(1, 30)

A <- as.numeric(t(one)  %*% Sinv %*% Rbar)
B <- as.numeric(t(Rbar) %*% Sinv %*% Rbar)
C <- as.numeric(t(one)  %*% Sinv %*% one)
D <- B * C - A^2

stopifnot(B > 0, C > 0, D > 0)            # Merton footnotes 4 and 5

cat("(a) A =", format(A, digits = 6), "\n")
cat("    B =", format(B, digits = 6), "\n")
cat("    C =", format(C, digits = 6), "\n")
cat("    D = BC - A^2 =", format(D, digits = 6), "\n")

# The global minimum variance portfolio sits at E = A/C with variance 1/C
E_gmv   <- A / C
sd_gmv  <- sqrt(1 / C)
cat(sprintf("    => global min variance portfolio: E = A/C = %.6f (%.2f%%), ",
            E_gmv, E_gmv * 100))
cat(sprintf("sd = sqrt(1/C) = %.6f (%.2f%%)\n\n", sd_gmv, sd_gmv * 100))

# ---- (b) Lagrange multipliers -----------------------------------------------

E0 <- 0.02        # chosen prescribed return: 2% per month, above A/C so efficient
stopifnot(E0 > E_gmv)

lam1 <- (C * E0 - A) / D
lam2 <- (B - A * E0) / D

cat(sprintf("(b) choosing E = %.4f (%.0f%% per month)\n", E0, E0 * 100))
cat("    lambda1 = (CE - A)/D =", format(lam1, digits = 6), "\n")
cat("    lambda2 = (B - AE)/D =", format(lam2, digits = 6), "\n\n")

# ---- (c) composition of the efficient portfolio with return E ---------------
# x = S^-1 [lambda1 Rbar + lambda2 1]

x_E <- Sinv %*% (lam1 * Rbar + lam2 * one)

# Cross-check against the g + hE form from the Results handout, part (a)
g <- as.numeric(Sinv %*% (B * one - A * Rbar) / D)
h <- as.numeric(Sinv %*% (C * Rbar - A * one) / D)
x_gh <- g + h * E0

stopifnot(all.equal(as.numeric(x_E), x_gh))     # the two forms agree
stopifnot(abs(sum(g) - 1) < 1e-9)               # sum(g) = (BC - A^2)/D = 1
stopifnot(abs(sum(h))     < 1e-9)               # sum(h) = (CA - AC)/D = 0
stopifnot(abs(sum(x_E) - 1) < 1e-9)             # fully invested
stopifnot(abs(as.numeric(t(x_E) %*% Rbar) - E0) < 1e-12)   # hits E

sd_E <- sqrt(as.numeric(t(x_E) %*% S %*% x_E))

cat("(c) efficient portfolio for E =", E0, "\n")
cat(sprintf("    mean = %.6f (%.2f%%)   sd = %.6f (%.2f%%)\n",
            E0, E0 * 100, sd_E, sd_E * 100))
cat(sprintf("    weights sum to %.6f, range %.3f to %.3f, %d short positions\n",
            sum(x_E), min(x_E), max(x_E), sum(x_E < 0)))
cat("    largest 5 weights:\n")
print(round(sort(setNames(as.numeric(x_E), stk), decreasing = TRUE)[1:5], 4))
cat("\n")

# Sanity: the closed form sigma^2 = (CE^2 - 2AE + B)/D must match x'Sx
stopifnot(abs(sd_E^2 - (C * E0^2 - 2 * A * E0 + B) / D) < 1e-14)

# ---- portfolios carried over from Project 1 ---------------------------------

w_eq  <- rep(1 / 30, 30)
mu_eq <- as.numeric(t(w_eq) %*% Rbar)
sd_eq <- sqrt(as.numeric(t(w_eq) %*% S %*% w_eq))

x_min  <- Sinv %*% one / C                      # = S^-1 1 / (1' S^-1 1)
mu_min <- as.numeric(t(x_min) %*% Rbar)
sd_min <- sqrt(as.numeric(t(x_min) %*% S %*% x_min))
stopifnot(abs(mu_min - E_gmv) < 1e-12,
          abs(sd_min - sd_gmv) < 1e-12)         # agrees with A/C and sqrt(1/C)

# ---- (g) three arbitrary portfolios -----------------------------------------

set.seed(183)
arb <- list(
  "G1" = local({ w <- numeric(30); w[1:6] <- 1 / 6; w }),        # technology only
  "G2" = local({ w <- runif(30); w / sum(w) }),                  # random, all long
  "G3" = local({ w <- 1 / 30 + rnorm(30, 0, 0.05); w / sum(w) }) # random long/short
)
arb_mu <- vapply(arb, function(w) as.numeric(t(w) %*% Rbar), numeric(1))
arb_sd <- vapply(arb, function(w) sqrt(as.numeric(t(w) %*% S %*% w)), numeric(1))
stopifnot(all(abs(vapply(arb, sum, numeric(1)) - 1) < 1e-12))

cat("(g) three arbitrary portfolios (weights sum to 1):\n")
print(round(data.frame(mean_pct = arb_mu * 100, sd_pct = arb_sd * 100), 3))
cat("\n")

# ---- frontier curves --------------------------------------------------------

Egrid   <- seq(min(means) - 0.012, max(means) + 0.012, length.out = 1200)
var_par <- (C * Egrid^2 - 2 * A * Egrid + B) / D
sd_hyp  <- sqrt(1 / C + (C / D) * (Egrid - E_gmv)^2)
stopifnot(max(abs(sd_hyp - sqrt(var_par))) < 1e-12)   # same curve, two algebras

# ---- plotting helpers -------------------------------------------------------

surface <- "#fcfcfb"; ink <- "#0b0b0b"; muted <- "#898781"; grid <- "#e1e0d9"
c_front <- "#2a78d6"   # frontier and the portfolios that sit on it
c_eq    <- "#eb6834"   # equal allocation
c_spx   <- "#1baf7a"   # S&P 500
c_arb   <- "#6f6d67"   # arbitrary portfolios

open_png <- function(f, h = 4.7) {
  png(f, width = 6.5, height = h, units = "in", res = 200)
  par(bg = surface, mar = c(4.2, 4.6, 3.4, 1.4), family = "sans")
}
ax <- function(side, at, lab) {
  axis(side, at = at, labels = lab, col = "#c3c2b7", col.axis = muted,
       cex.axis = 0.75, lwd = 1, las = if (side == 2) 1 else 0)
}

# ---- (d) parabola in mean-variance space ------------------------------------

open_png("fig_parabola.png")
xr <- range(Egrid); yr <- c(0, max(var_par) * 1.04)
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)
lines(Egrid, var_par, col = c_front, lwd = 2.2)
points(E_gmv, 1 / C, pch = 21, bg = c_front, col = surface, lwd = 2, cex = 1.5)
text(E_gmv, 1 / C, "minimum variance, (A/C, 1/C)", pos = 1, offset = 0.7,
     cex = 0.72, font = 2, col = c_front)
ax(1, pretty(xr), sprintf("%.0f%%", pretty(xr) * 100))
ax(2, pretty(yr), format(pretty(yr), digits = 2))
mtext("Expected return E (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext(expression(paste("Variance  ", sigma^2)), 2, line = 3.2, col = muted, cex = 0.85)
mtext("(d) Frontier in mean-variance space", 3, line = 1.8, adj = 0,
      col = ink, cex = 1.05, font = 2)
mtext(expression(paste(sigma^2, " = (C", E^2, " - 2AE + B)/D  -  a parabola in E")),
      3, line = 0.5, adj = 0, col = muted, cex = 0.8)
dev.off()

# ---- (e), (f), (g) hyperbola in mean-sd space -------------------------------

open_png("fig_frontier.png", h = 5.2)
is_g <- colnames(ret) == "^GSPC"
xr <- c(0, max(stdev) * 1.06)
yr <- c(min(means) - 0.009, max(means) + 0.007)

plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)

# asymptotes E = A/C +- sqrt(D/C) sigma   (Merton eq. 16)
abline(a = E_gmv, b =  sqrt(D / C), col = "#c9c8c0", lwd = 1, lty = 2)
abline(a = E_gmv, b = -sqrt(D / C), col = "#c9c8c0", lwd = 1, lty = 2)

eff <- Egrid >= E_gmv
lines(sd_hyp[!eff], Egrid[!eff], col = c_front, lwd = 1.1, lty = 3)   # inefficient
lines(sd_hyp[ eff], Egrid[ eff], col = c_front, lwd = 2.4)            # efficient

points(stdev[!is_g], means[!is_g], pch = 21, bg = "#d9d8d2",
       col = surface, lwd = 1, cex = 0.8)

pt <- function(x, y, col, lab, pos = 4, cex = 0.72, off = 0.6) {
  points(x, y, pch = 21, bg = col, col = surface, lwd = 1.8, cex = 1.5)
  text(x, y, lab, pos = pos, offset = off, cex = cex, font = 2, col = col)
}
points(arb_sd, arb_mu, pch = 21, bg = c_arb, col = surface, lwd = 1.6, cex = 1.3)
text(arb_sd, arb_mu, names(arb), pos = c(4, 1, 3), offset = 0.5, cex = 0.7,
     font = 2, col = c_arb)
pt(stdev[is_g], means[is_g], c_spx,   "S&P 500", 1)
pt(sd_eq,  mu_eq,  c_eq,    "equal allocation", 4)
pt(sd_min, mu_min, c_front, "minimum risk", 1)
pt(sd_E,   E0,     c_front, "E = 2%", 2)

ax(1, pretty(xr), sprintf("%.0f%%", pretty(xr) * 100))
ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(e)-(g) Frontier in mean-standard deviation space", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.05, font = 2)
mtext("Hyperbola, with the 30 stocks, the index, and four portfolios. Dashed lines are the asymptotes.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

# ---- (h) mutual fund theorem ------------------------------------------------
# Two efficient portfolios span the whole frontier. Take the one from (c) and a
# second efficient portfolio, then combine them as in handout #4 page 5:
#   E_p    = a E_1 + (1-a) E_2
#   sd_p^2 = a^2 s1^2 + (1-a)^2 s2^2 + 2 a (1-a) s12

E1 <- E0;   x1 <- as.numeric(x_E)
E2 <- 0.035; lam1b <- (C * E2 - A) / D; lam2b <- (B - A * E2) / D
x2 <- as.numeric(Sinv %*% (lam1b * Rbar + lam2b * one))
stopifnot(abs(sum(x2) - 1) < 1e-9, abs(as.numeric(t(x2) %*% Rbar) - E2) < 1e-12)

s1  <- sqrt(as.numeric(t(x1) %*% S %*% x1))
s2  <- sqrt(as.numeric(t(x2) %*% S %*% x2))
s12 <- as.numeric(t(x1) %*% S %*% x2)

# The Results handout, part (c), gives the covariance of two frontier
# portfolios in closed form. It must agree with x1' S x2.
s12_closed <- (C / D) * (E1 - A / C) * (E2 - A / C) + 1 / C
stopifnot(abs(s12 - s12_closed) < 1e-14)

agrid <- seq(-1.5, 2.5, length.out = 800)
mf_E  <- agrid * E1 + (1 - agrid) * E2
mf_sd <- sqrt(agrid^2 * s1^2 + (1 - agrid)^2 * s2^2 + 2 * agrid * (1 - agrid) * s12)

# Every combination must land exactly on the frontier: check sigma against the
# closed form for its own expected return.
stopifnot(max(abs(mf_sd^2 - (C * mf_E^2 - 2 * A * mf_E + B) / D)) < 1e-15)

cat("(h) mutual fund theorem\n")
cat(sprintf("    fund 1: E = %.4f, sd = %.4f   fund 2: E = %.4f, sd = %.4f\n",
            E1, s1, E2, s2))
cat(sprintf("    cov(1,2) = %.8f from x1'Sx2, %.8f from (C/D)(E1-A/C)(E2-A/C)+1/C\n",
            s12, s12_closed))
cat(sprintf("    max |combination - frontier| over %d mixes: %.2e\n",
            length(agrid), max(abs(mf_sd - sqrt((C * mf_E^2 - 2 * A * mf_E + B) / D)))))

open_png("fig_mutualfund.png")
sel <- mf_E >= yr[1] & mf_E <= yr[2]
xr2 <- c(0, max(mf_sd[sel]) * 1.08)
plot(NA, xlim = xr2, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = pretty(xr2), col = grid, lwd = 1)
lines(sd_hyp, Egrid, col = c_front, lwd = 4)
points(mf_sd[sel][seq(1, sum(sel), by = 14)], mf_E[sel][seq(1, sum(sel), by = 14)],
       pch = 21, bg = surface, col = c_eq, lwd = 1.6, cex = 0.85)
pt(s1, E1, c_front, "fund 1  (E = 2%)", 2)
pt(s2, E2, c_front, "fund 2  (E = 3.5%)", 2)
ax(1, pretty(xr2), sprintf("%.0f%%", pretty(xr2) * 100))
ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(h) Two funds trace the whole frontier", 3, line = 1.8, adj = 0,
      col = ink, cex = 1.05, font = 2)
mtext("Thick line: frontier from (e). Circles: combinations of funds 1 and 2 only.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.75)
dev.off()

cat("\nwrote fig_parabola.png, fig_frontier.png, fig_mutualfund.png\n")
