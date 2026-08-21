# Project 4, parts 1-3: the single index model
#
#   R_i = alpha_i + beta_i R_m + epsilon_i

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "single_index.R")
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
ret <- P[-1, ] / P[-n, ] - 1              # 59 x 31

stk    <- setdiff(colnames(ret), "^GSPC")
ret30  <- ret[, stk]
Rm     <- as.numeric(ret[, "^GSPC"])      # the index, our single factor
Rbar   <- colMeans(ret30)
S_hist <- cov(ret30)                      # Project 2's historical matrix

Tn <- nrow(ret30); N <- ncol(ret30)
stopifnot(Tn == 59, N == 30)

# ---- 1. regress each stock on the S&P 500 -----------------------------------

fit   <- lapply(stk, function(s) lm(ret30[, s] ~ Rm))
names(fit) <- stk
alpha <- vapply(fit, function(f) unname(coef(f)[1]), numeric(1))
beta  <- vapply(fit, function(f) unname(coef(f)[2]), numeric(1))
sig2e <- vapply(fit, function(f) summary(f)$sigma^2, numeric(1))   # RSS/(T-2)
r2    <- vapply(fit, function(f) summary(f)$r.squared, numeric(1))

stopifnot(max(abs(beta - apply(ret30, 2, function(x) cov(x, Rm) / var(Rm)))) < 1e-12)
stopifnot(max(abs(alpha - (Rbar - beta * mean(Rm)))) < 1e-15)
stopifnot(max(abs(sig2e / vapply(fit, function(f) sum(residuals(f)^2) / (Tn - 2),
                                 numeric(1)) - 1)) < 1e-14)

stopifnot(max(abs(alpha + beta * mean(Rm) - Rbar)) < 1e-16)

cat("1. single index model estimates (monthly)\n")
cat(sprintf("   market mean %.4f%%, market sd %.4f%%\n",
            mean(Rm) * 100, sd(Rm) * 100))
cat(sprintf("   beta  range %.3f to %.3f (mean %.3f)\n",
            min(beta), max(beta), mean(beta)))
cat(sprintf("   alpha range %.4f%% to %.4f%% per month\n",
            min(alpha) * 100, max(alpha) * 100))
cat(sprintf("   R^2   range %.3f to %.3f (mean %.3f)\n\n",
            min(r2), max(r2), mean(r2)))

sim_tbl <- data.frame(alpha = alpha, beta = beta, sigma2_e = sig2e, R2 = r2)

# ---- 2. the single index variance-covariance matrix -------------------------

sig2m <- var(Rm)
S_sim <- outer(beta, beta) * sig2m
diag(S_sim) <- beta^2 * sig2m + sig2e
dimnames(S_sim) <- list(stk, stk)

stopifnot(isSymmetric(S_sim))
stopifnot(all.equal(diag(S_sim), beta^2 * sig2m + sig2e))

ev_sim  <- eigen(S_sim,  only.values = TRUE)$values
ev_hist <- eigen(S_hist, only.values = TRUE)$values
stopifnot(min(ev_sim) > 0)

cat("2. variance-covariance matrix from the single index model\n")
cat(sprintf("   dim %d x %d, symmetric, smallest eigenvalue %.3e\n",
            nrow(S_sim), ncol(S_sim), min(ev_sim)))
cat(sprintf("   condition number: single index %.1f   historical %.1f\n",
            max(ev_sim) / min(ev_sim), max(ev_hist) / min(ev_hist)))
cat(sprintf("   free parameters: single index %d, historical %d, from %d months\n",
            2 * N + 1, N * (N + 1) / 2, Tn))
cat(sprintf("   mean |off-diagonal|: single index %.6f, historical %.6f\n\n",
            mean(abs(S_sim[upper.tri(S_sim)])),
            mean(abs(S_hist[upper.tri(S_hist)]))))

# ---- 3. frontier from each matrix -------------------------------------------

frontier <- function(Sig, mu) {
  Si  <- solve(Sig); one <- rep(1, length(mu))
  A <- as.numeric(t(one) %*% Si %*% mu)
  B <- as.numeric(t(mu)  %*% Si %*% mu)
  C <- as.numeric(t(one) %*% Si %*% one)
  D <- B * C - A^2
  stopifnot(B > 0, C > 0, D > 0)
  list(A = A, B = B, C = C, D = D, Emin = A / C, sdmin = sqrt(1 / C),
       xmin = as.numeric(Si %*% one / C))
}
f_hist <- frontier(S_hist, Rbar)
f_sim  <- frontier(S_sim,  Rbar)

# The historical numbers must reproduce Project 2 part (a) exactly
stopifnot(abs(f_hist$A - 36.9864) < 1e-3, abs(f_hist$B - 1.25629) < 1e-4,
          abs(f_hist$C - 2330.03) < 1e-1, abs(f_hist$D - 1559.19) < 1e-1)

sd_at <- function(f, E) sqrt(1 / f$C + (f$C / f$D) * (E - f$A / f$C)^2)

# The hyperbola and parabola forms must agree, for both matrices
for (fr in list(f_hist, f_sim)) {
  Eg <- seq(-0.02, 0.09, length.out = 500)
  stopifnot(max(abs(sd_at(fr, Eg)^2 -
                    (fr$C * Eg^2 - 2 * fr$A * Eg + fr$B) / fr$D)) < 1e-15)
}
stopifnot(abs(sum(f_hist$xmin) - 1) < 1e-9, abs(sum(f_sim$xmin) - 1) < 1e-9)

cat("3. frontier inputs\n")
cat(sprintf("   %-13s %9s %9s %9s %9s %9s %9s\n",
            "matrix", "A", "B", "C", "D", "E min", "sd min"))
cat(sprintf("   %-13s %9.4f %9.5f %9.2f %9.2f %8.3f%% %8.3f%%\n",
            "historical", f_hist$A, f_hist$B, f_hist$C, f_hist$D,
            f_hist$Emin * 100, f_hist$sdmin * 100))
cat(sprintf("   %-13s %9.4f %9.5f %9.2f %9.2f %8.3f%% %8.3f%%\n\n",
            "single index", f_sim$A, f_sim$B, f_sim$C, f_sim$D,
            f_sim$Emin * 100, f_sim$sdmin * 100))

E_cmp <- 0.02
gap   <- sd_at(f_sim, E_cmp) - sd_at(f_hist, E_cmp)
cat(sprintf("   at E = %.0f%%: historical %.4f%%, single index %.4f%% (+%.4f pp)\n",
            E_cmp * 100, sd_at(f_hist, E_cmp) * 100,
            sd_at(f_sim, E_cmp) * 100, gap * 100))

dgap  <- function(E) sd_at(f_sim, E) - sd_at(f_hist, E)
xhi   <- uniroot(dgap, c(0.02, 0.30), tol = 1e-12)$root
xlo   <- uniroot(dgap, c(-0.30, 0.01), tol = 1e-12)$root
stopifnot(dgap(mean(c(xlo, xhi))) > 0, dgap(xhi + 0.01) < 0)

cat(sprintf("   asymptote slope sqrt(D/C): historical %.4f, single index %.4f\n",
            sqrt(f_hist$D / f_hist$C), sqrt(f_sim$D / f_sim$C)))
cat(sprintf("   the curves cross at E = %.2f%% and E = %.2f%%; between them the\n",
            xlo * 100, xhi * 100))
cat(sprintf("   single index frontier is to the right, outside them to the left\n\n"))

real <- function(x) sqrt(as.numeric(t(x) %*% S_hist %*% x))
stopifnot(real(f_sim$xmin) >= real(f_hist$xmin))   # historical is the true floor

cat("   minimum risk portfolio, scored under the historical matrix:\n")
cat(sprintf("     historical weights   %.4f%%   (the in-sample floor)\n",
            real(f_hist$xmin) * 100))
cat(sprintf("     single index weights %.4f%%   (single index claims %.4f%%,\n",
            real(f_sim$xmin) * 100, f_sim$sdmin * 100))
cat(sprintf("       so it understates its own risk by %.4f pp)\n",
            (real(f_sim$xmin) - f_sim$sdmin) * 100))
cat(sprintf("   but the weights are far less extreme: sum|w| %.2f vs %.2f,\n",
            sum(abs(f_sim$xmin)), sum(abs(f_hist$xmin))))
cat(sprintf("     largest %.3f vs %.3f, %d shorts vs %d\n\n",
            max(f_sim$xmin), max(f_hist$xmin),
            sum(f_sim$xmin < 0), sum(f_hist$xmin < 0)))

# ---- plots ------------------------------------------------------------------

surface <- "#fcfcfb"; ink <- "#0b0b0b"; muted <- "#898781"; grid <- "#e1e0d9"
c_hist <- "#2a78d6"    # historical covariance matrix (Project 2)
c_sim  <- "#eb6834"    # single index covariance matrix
c_spx  <- "#1baf7a"

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

means_all <- colMeans(ret); sd_all <- sqrt(diag(cov(ret)))
is_g <- colnames(ret) == "^GSPC"

draw_both <- function(xr, yr, title, sub) {
  Eg <- seq(yr[1], yr[2], length.out = 1500)
  plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)
  for (fr in list(list(f_hist, c_hist), list(f_sim, c_sim))) {
    s <- sd_at(fr[[1]], Eg); eff <- Eg >= fr[[1]]$Emin
    lines(s[!eff], Eg[!eff], col = fr[[2]], lwd = 1.2, lty = 3)
    lines(s[ eff], Eg[ eff], col = fr[[2]], lwd = 2.6)
  }
  ax(1, pretty(xr), sprintf("%.0f%%", pretty(xr) * 100))
  ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
  mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
  mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
  mtext(title, 3, line = 1.8, adj = 0, col = ink, cex = 1.02, font = 2)
  mtext(sub, 3, line = 0.5, adj = 0, col = muted, cex = 0.72)
}

open_png("fig_frontiers.png", h = 5.2)
xr <- c(0, max(sd_all) * 1.06)
yr <- c(min(means_all) - 0.009, max(means_all) + 0.007)
draw_both(xr, yr, "Two frontiers: same means, different covariance matrices",
          "Solid: efficient half. Dotted: inefficient half. Grey: the 30 stocks.")
points(sd_all[!is_g], means_all[!is_g], pch = 21, bg = "#d9d8d2",
       col = surface, lwd = 1, cex = 0.8)
points(sd_at(f_hist, xhi), xhi, pch = 4, col = ink, lwd = 2, cex = 1.1)
text(sd_at(f_hist, xhi), xhi, sprintf("curves cross, E = %.2f%%", xhi * 100),
     pos = 4, offset = 0.5, cex = 0.66, font = 2, col = ink)
pt(sd_all[is_g], means_all[is_g], c_spx, "S&P 500", 1)
pt(f_hist$sdmin, f_hist$Emin, c_hist, "historical", 2)
pt(f_sim$sdmin,  f_sim$Emin,  c_sim,  "single index", 4)
dev.off()

open_png("fig_frontiers_zoom.png")
xr2 <- c(0, max(f_sim$sdmin, f_hist$sdmin) * 2.3)
yr2 <- c(min(f_hist$Emin, f_sim$Emin) - 0.014, 0.05)
draw_both(xr2, yr2, "The same two frontiers, near the vertex",
          "Over the range the data occupies, the single index frontier sits to the right.")
segments(sd_at(f_hist, E_cmp), E_cmp, sd_at(f_sim, E_cmp), E_cmp,
         col = ink, lwd = 1.4)
text(sd_at(f_sim, E_cmp), E_cmp, sprintf("+%.2f pp at E = 2%%", gap * 100),
     pos = 4, offset = 0.4, cex = 0.68, font = 2, col = ink)
pt(f_hist$sdmin, f_hist$Emin, c_hist, "historical", 2)
pt(f_sim$sdmin,  f_sim$Emin,  c_sim,  "single index", 4)
dev.off()

cat("wrote fig_frontiers.png, fig_frontiers_zoom.png\n")
