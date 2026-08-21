local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "tangency.R")
    setwd(dirname(normalizePath(f)))
})

# ---- (a) prices to returns --------------------------------------------------

a <- read.table("statc183c283_5stocks.txt", header = TRUE)
stopifnot(nrow(a) == 216, names(a)[1] == "date")

a <- a[order(a$date), ]                   # oldest -> newest
dates <- as.Date(as.character(a$date), format = "%Y%m%d")
P <- as.matrix(a[, -1])
colnames(P) <- c("XOM", "GM", "HPQ", "MCD", "BA")
stopifnot(!anyNA(P), !is.unsorted(dates))

n   <- nrow(P)
ret <- P[-1, ] / P[-n, ] - 1              # 215 x 5 simple monthly returns
stopifnot(nrow(ret) == 215, ncol(ret) == 5)

# ---- (b) mean returns and variance-covariance matrix ------------------------

Rbar <- colMeans(ret)
S    <- cov(ret)
stopifnot(qr(S)$rank == 5)


cat("(b) mean monthly return, in percent:\n")
print(round(Rbar * 100, 4))
cat("\n    variance-covariance matrix:\n")
print(round(S, 6))
cat("\n")

# ---- (c) minimum risk portfolio, Exxon-Mobil and Boeing ---------------------
# Two assets, so the weight follows from setting d(sigma^2)/dx = 0:
#   x_XOM = (sd_BA^2 - cov) / (sd_XOM^2 + sd_BA^2 - 2 cov)

stk2  <- c("XOM", "BA")
S2    <- S[stk2, stk2]
Rbar2 <- Rbar[stk2]

vx <- S2["XOM", "XOM"]; vb <- S2["BA", "BA"]; cxb <- S2["XOM", "BA"]
x2min <- c(XOM = (vb - cxb) / (vx + vb - 2 * cxb))
x2min["BA"] <- 1 - x2min["XOM"]

# The general form Sigma^-1 1 / (1' Sigma^-1 1) must give the same answer
stopifnot(all.equal(as.numeric(x2min),
                    as.numeric(solve(S2) %*% rep(1, 2) / sum(solve(S2)))))

E2min  <- as.numeric(t(x2min) %*% Rbar2)
sd2min <- sqrt(as.numeric(t(x2min) %*% S2 %*% x2min))

# It has to be at least as good as holding either stock alone
stopifnot(sd2min < sqrt(vx), sd2min < sqrt(vb))

cat("(c) minimum risk portfolio, Exxon-Mobil and Boeing\n")
cat(sprintf("    x_XOM = %.6f   x_BA = %.6f\n", x2min["XOM"], x2min["BA"]))
cat(sprintf("    E = %.6f (%.4f%%)   sd = %.6f (%.4f%%)\n\n",
            E2min, E2min * 100, sd2min, sd2min * 100))

# ---- (d) portfolio possibilities curve for the two stocks -------------------

x2grid <- seq(-0.6, 1.6, length.out = 900)
W2     <- cbind(x2grid, 1 - x2grid)
E2c    <- as.numeric(W2 %*% Rbar2)
sd2c   <- sqrt(rowSums((W2 %*% S2) * W2))
stopifnot(min(sd2c) >= sd2min - 1e-12)     # nothing beats the minimum

# ---- (e) cloud of points, Exxon-Mobil / McDonalds / Boeing ------------------
# 2499 combinations of (xa, xb, xc) summing to 1, short sales allowed

abc <- read.table("statc183c283_abc.txt", header = TRUE)
# weights are stored to 6 decimals, so they sum to 1 only to ~1e-9
stopifnot(ncol(abc) == 3, max(abs(rowSums(abc) - 1)) < 1e-8)

stk3  <- c("XOM", "MCD", "BA")
Rbar3 <- Rbar[stk3]
S3    <- S[stk3, stk3]
S3inv <- solve(S3)
one3  <- rep(1, 3)

X3       <- as.matrix(abc)
E_cloud  <- as.numeric(X3 %*% Rbar3)
sd_cloud <- sqrt(rowSums((X3 %*% S3) * X3))

cat("(e) cloud of", nrow(X3), "portfolios\n")
cat(sprintf("    E  ranges %.5f to %.5f\n", min(E_cloud), max(E_cloud)))
cat(sprintf("    sd ranges %.5f to %.5f\n\n", min(sd_cloud), max(sd_cloud)))

# ---- (f) point of tangency G, Rf = 0.001 ------------------------------------
# Maximising the slope (E_p - Rf)/sigma_p gives z = Sigma^-1 (Rbar - Rf 1),
# then rescale to weights: x = z / sum(z)

Rf1 <- 0.001

tangency <- function(Rf) {
  z <- S3inv %*% (Rbar3 - Rf * one3)
  x <- as.numeric(z / sum(z))
  setNames(x, stk3)
}

xG      <- tangency(Rf1)
E_G     <- as.numeric(t(xG) %*% Rbar3)
sd_G    <- sqrt(as.numeric(t(xG) %*% S3 %*% xG))
slope_G <- (E_G - Rf1) / sd_G

stopifnot(abs(sum(xG) - 1) < 1e-12)
# G must have the highest Sharpe ratio in the whole cloud
stopifnot(slope_G >= max((E_cloud - Rf1) / sd_cloud))

cat("(f) point of tangency G, Rf =", Rf1, "\n")
print(round(xG, 6))
cat(sprintf("    E = %.6f (%.4f%%)   sd = %.6f (%.4f%%)   slope = %.6f\n\n",
            E_G, E_G * 100, sd_G, sd_G * 100, slope_G))

# ---- (g) 60% G, 40% risk free asset -----------------------------------------

wG   <- 0.6
E_g  <- wG * E_G + (1 - wG) * Rf1
sd_g <- wG * sd_G                          # the risk free asset adds no variance

stopifnot(abs(E_g - (Rf1 + slope_G * sd_g)) < 1e-15)   # sits on the CAL

cat("(g) 60% G + 40% risk free\n")
cat(sprintf("    E = %.6f (%.4f%%)   sd = %.6f (%.4f%%)\n\n",
            E_g, E_g * 100, sd_g, sd_g * 100))

# ---- (h) x = (E - Rf) Sigma^-1 (Rbar - Rf 1) / (Rbar - Rf 1)' Sigma^-1 (...) -

excess <- Rbar3 - Rf1 * one3
Hq     <- as.numeric(t(excess) %*% S3inv %*% excess)
x_h    <- as.numeric((E_g - Rf1) * S3inv %*% excess / Hq)
names(x_h) <- stk3

# x is the holding of each *stock* in the combined portfolio of (g): it is
# exactly 0.6 (the weight on G) times G's own composition, and what is left
# over, 1 - sum(x) = 0.4, is the risk free asset.
stopifnot(all.equal(x_h, wG * xG))
stopifnot(abs(sum(x_h) - wG) < 1e-12)
stopifnot(abs(as.numeric(t(x_h) %*% Rbar3) + (1 - wG) * Rf1 - E_g) < 1e-15)
stopifnot(abs(sqrt(as.numeric(t(x_h) %*% S3 %*% x_h)) - sd_g) < 1e-15)

cat("(h) x for E =", format(E_g, digits = 6), "\n")
print(round(x_h, 6))
cat(sprintf("    sum(x) = %.6f, so %.0f%% in the three stocks and %.0f%% risk free\n",
            sum(x_h), sum(x_h) * 100, (1 - sum(x_h)) * 100))
cat(sprintf("    max |x - 0.6 * xG| = %.2e\n\n", max(abs(x_h - wG * xG))))

# ---- (i) short sales allowed, no risk free asset ----------------------------

Rf2 <- 0.002

# (i)1 two portfolios tangent to the frontier
xA <- tangency(Rf1)                        # same portfolio as G in (f)
xB <- tangency(Rf2)
stopifnot(identical(xA, xG))

E_A  <- as.numeric(t(xA) %*% Rbar3); sd_A <- sqrt(as.numeric(t(xA) %*% S3 %*% xA))
E_B  <- as.numeric(t(xB) %*% Rbar3); sd_B <- sqrt(as.numeric(t(xB) %*% S3 %*% xB))
stopifnot(abs(sum(xB) - 1) < 1e-12)

cat("(i)1 two tangency portfolios\n")
print(round(rbind(A = xA, B = xB), 6))
cat(sprintf("     A: E = %.6f, sd = %.6f   (Rf = %.3f)\n", E_A, sd_A, Rf1))
cat(sprintf("     B: E = %.6f, sd = %.6f   (Rf = %.3f)\n\n", E_B, sd_B, Rf2))

# (i)2 covariance between A and B
sig_AB <- as.numeric(t(xA) %*% S3 %*% xB)

# Merton's constants for the three stocks, used both to check sig_AB against
# the closed form from handout #11 and to check the traced frontier below.
A3 <- as.numeric(t(one3)  %*% S3inv %*% Rbar3)
B3 <- as.numeric(t(Rbar3) %*% S3inv %*% Rbar3)
C3 <- as.numeric(t(one3)  %*% S3inv %*% one3)
D3 <- B3 * C3 - A3^2
stopifnot(B3 > 0, C3 > 0, D3 > 0)

sig_AB_closed <- (C3 / D3) * (E_A - A3 / C3) * (E_B - A3 / C3) + 1 / C3
stopifnot(abs(sig_AB - sig_AB_closed) < 1e-14)

cat("(i)2 covariance between A and B\n")
cat(sprintf("     xA' Sigma xB                           = %.10f\n", sig_AB))
cat(sprintf("     (C/D)(E_A - A/C)(E_B - A/C) + 1/C      = %.10f\n", sig_AB_closed))
cat(sprintf("     correlation = %.6f\n\n", sig_AB / (sd_A * sd_B)))

# (i)3 trace the frontier out of A and B
#   E_p    = a E_A + (1-a) E_B
#   sd_p^2 = a^2 sd_A^2 + (1-a)^2 sd_B^2 + 2 a (1-a) sig_AB

Elo <- min(E_cloud) - 0.001; Ehi <- max(E_cloud) + 0.001
agrid <- sort(c((Elo - E_B) / (E_A - E_B), (Ehi - E_B) / (E_A - E_B)))
agrid <- seq(agrid[1], agrid[2], length.out = 2000)

fr_E  <- agrid * E_A + (1 - agrid) * E_B
fr_sd <- sqrt(agrid^2 * sd_A^2 + (1 - agrid)^2 * sd_B^2 +
              2 * agrid * (1 - agrid) * sig_AB)

# Every mix of A and B must land on the frontier of the three stocks
fr_closed <- sqrt((C3 * fr_E^2 - 2 * A3 * fr_E + B3) / D3)
stopifnot(max(abs(fr_sd - fr_closed)) < 1e-14)

# and the frontier must sit to the left of every point in the cloud
cloud_min_sd <- sqrt((C3 * E_cloud^2 - 2 * A3 * E_cloud + B3) / D3)
stopifnot(all(sd_cloud >= cloud_min_sd - 1e-12))

cat("(i)3 frontier traced from A and B\n")
cat(sprintf("     max |two-fund mix - closed form| over %d mixes: %.2e\n",
            length(agrid), max(abs(fr_sd - fr_closed))))
cat(sprintf("     min slack of the cloud over the frontier: %.2e\n\n",
            min(sd_cloud - cloud_min_sd)))

# (i)4 minimum risk portfolio of the three stocks
x3min  <- setNames(as.numeric(S3inv %*% one3 / C3), stk3)
E3min  <- A3 / C3
sd3min <- sqrt(1 / C3)

stopifnot(abs(sum(x3min) - 1) < 1e-12)
stopifnot(abs(as.numeric(t(x3min) %*% Rbar3) - E3min) < 1e-15)
stopifnot(abs(sqrt(as.numeric(t(x3min) %*% S3 %*% x3min)) - sd3min) < 1e-15)
stopifnot(sd3min <= min(sd_cloud))         # the leftmost attainable point

cat("(i)4 minimum risk portfolio, three stocks\n")
print(round(x3min, 6))
cat(sprintf("     E = A/C = %.6f (%.4f%%)   sd = sqrt(1/C) = %.6f (%.4f%%)\n\n",
            E3min, E3min * 100, sd3min, sd3min * 100))
# ---- plotting helpers -------------------------------------------------------
# Same palette and helpers as Project 2, so the figures match across projects.

surface <- "#fcfcfb"; ink <- "#0b0b0b"; muted <- "#898781"; grid <- "#e1e0d9"
c_front <- "#2a78d6"   # frontier and the portfolios that sit on it
c_eq    <- "#eb6834"   # tangency portfolio, CAL
c_spx   <- "#1baf7a"   # the traced frontier in (i)3
c_arb   <- "#6f6d67"   # individual stocks

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

# ---- (c) and (d): the two stock possibilities curve -------------------------

open_png("fig_two_stock.png")
xr <- c(0.05, max(sd2c) * 1.02); yr <- range(E2c) + c(-0.0006, 0.0006)
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)

eff2 <- E2c >= E2min
lines(sd2c[!eff2], E2c[!eff2], col = c_front, lwd = 1.2, lty = 3)   # inefficient
lines(sd2c[ eff2], E2c[ eff2], col = c_front, lwd = 2.6)            # efficient

pt(sqrt(vx), Rbar2["XOM"], c_arb, "Exxon-Mobil", 2)
pt(sqrt(vb), Rbar2["BA"],  c_arb, "Boeing", 2)
pt(sd2min, E2min, c_front, "minimum risk", 4)

ax(1, pretty(xr), sprintf("%.0f%%", pretty(xr) * 100))
ax(2, pretty(yr), sprintf("%.1f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(c)-(d) Portfolio possibilities curve, Exxon-Mobil and Boeing", 3,
      line = 1.8, adj = 0, col = ink, cex = 1.02, font = 2)
mtext("Solid arm is the efficient frontier; dotted arm is dominated by it.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.75)
dev.off()

# ---- (e), (f), (g): cloud, tangency portfolio, capital allocation line ------

open_png("fig_cal.png", h = 5.2)
# left padding so the labels on the low-risk portfolios are not clipped
xr <- c(-0.04, max(sd_cloud) * 1.03)
xtick <- pretty(c(0, max(sd_cloud)))
yr <- c(min(E_cloud) - 0.0012, max(E_cloud) + 0.0018)
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = xtick, col = grid, lwd = 1)

points(sd_cloud, E_cloud, pch = 16, col = "#c9c8c0", cex = 0.42)

# the capital allocation line: E = Rf + slope * sigma
abline(a = Rf1, b = slope_G, col = c_eq, lwd = 2)
text(xr[2] * 0.62, Rf1 + slope_G * xr[2] * 0.62, "CAL", pos = 3, offset = 0.4,
     cex = 0.75, font = 2, col = c_eq)

pt(0, Rf1, c_eq, expression(R[f] == 0.001), 4, off = 0.5)
pt(sd_G, E_G, c_eq, "G (tangency)", 4)
pt(sd_g, E_g, c_front, "60% G + 40% Rf", 2)

ax(1, xtick, sprintf("%.0f%%", xtick * 100))
ax(2, pretty(yr), sprintf("%.1f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(e)-(g) Cloud, point of tangency, and the capital allocation line", 3,
      line = 1.8, adj = 0, col = ink, cex = 1.02, font = 2)
mtext(sprintf("%d combinations of Exxon-Mobil, McDonalds and Boeing with weights summing to 1.",
              nrow(X3)), 3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

# ---- (i)3: the frontier traced from A and B, drawn over the cloud -----------

open_png("fig_frontier3.png", h = 5.2)
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = xtick, col = grid, lwd = 1)

points(sd_cloud, E_cloud, pch = 16, col = "#333331", cex = 0.42)

eff3 <- fr_E >= E3min
lines(fr_sd[!eff3], fr_E[!eff3], col = c_spx, lwd = 2.4)
lines(fr_sd[ eff3], fr_E[ eff3], col = c_spx, lwd = 3.2)

pt(sd_A, E_A, c_eq,    "A", 2, off = 0.5)
pt(sd_B, E_B, c_eq,    "B", 2, off = 0.5)
pt(sd3min, E3min, c_front, "minimum risk", 2)

ax(1, xtick, sprintf("%.0f%%", xtick * 100))
ax(2, pretty(yr), sprintf("%.1f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(i) Efficient frontier traced from portfolios A and B", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.02, font = 2)
mtext("Green curve is combinations of A and B only; it envelopes the cloud from the left.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.72)
dev.off()

cat("wrote fig_two_stock.png, fig_cal.png, fig_frontier3.png\n")
