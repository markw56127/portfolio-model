# Project 6, parts a-e: the Elton-Gruber-Padberg ranking procedure under the
# single index model and under the constant correlation model, with and without
# short sales.

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "optimization.R")
    setwd(dirname(normalizePath(f)))
})

DATA <- ".."

# ---- data -------------------------------------------------------------------
# The training file, the same 59 returns Projects 1, 2, 4 and 5 use.

a <- read.csv(file.path(DATA, "stockData_train.csv"), sep = ",", header = TRUE,
              check.names = FALSE)
stopifnot(!any(duplicated(names(a))))

P <- as.matrix(a[, -1][, -1])
stopifnot(nrow(P) == 60, ncol(P) == 31, !anyNA(P))

n     <- nrow(P)
ret   <- P[-1, ] / P[-n, ] - 1
stk   <- setdiff(colnames(ret), "^GSPC")
ret30 <- ret[, stk]
Rm    <- as.numeric(ret[, "^GSPC"])
stopifnot(nrow(ret) == 59, length(stk) == 30)

Rbar   <- colMeans(ret30)
sdv    <- apply(ret30, 2, sd)
S_hist <- cov(ret30)

fit   <- lapply(stk, function(s) lm(ret30[, s] ~ Rm))
names(fit) <- stk
beta  <- vapply(fit, function(f) unname(coef(f)[2]), numeric(1))
sig2e <- vapply(fit, function(f) summary(f)$sigma^2, numeric(1))
sig2m <- var(Rm)

S_sim <- outer(beta, beta) * sig2m
diag(S_sim) <- beta^2 * sig2m + sig2e
dimnames(S_sim) <- list(stk, stk)

Rf <- 0.001                          # the same Rf as Project 5

pos <- beta > 0                      # "use only the stocks with positive betas"
stopifnot(all(pos))                  # all 30 qualify, so nothing is dropped

# Project 2's frontier, for the plots
frontier <- function(Sig, mu) {
  Si <- solve(Sig); one <- rep(1, length(mu))
  A <- as.numeric(t(one) %*% Si %*% mu); B <- as.numeric(t(mu) %*% Si %*% mu)
  C <- as.numeric(t(one) %*% Si %*% one); D <- B * C - A^2
  list(A = A, B = B, C = C, D = D, Emin = A / C, sdmin = sqrt(1 / C))
}
f_hist <- frontier(S_hist, Rbar)
stopifnot(abs(f_hist$A - 36.9864) < 1e-3, abs(f_hist$C - 2330.03) < 1e-1)
sd_at <- function(f, E) sqrt(1 / f$C + (f$C / f$D) * (E - f$A / f$C)^2)

# round the numeric columns of a table for printing
rnd <- function(d, k = 6) {
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], round, k)
  d
}

# portfolio statistics under either covariance matrix
stats <- function(x, Sig) {
  w <- rep(0, 30); names(w) <- stk; w[names(x)] <- x
  c(E = sum(w * Rbar), sd = sqrt(as.numeric(t(w) %*% Sig %*% w)))
}

# ---- (a) ranking by excess return to beta, handout #28 ----------------------

egp_table <- function(Rf) {
  ratio <- (Rbar - Rf) / beta
  ord   <- order(ratio, decreasing = TRUE)
  num   <- (Rbar - Rf) * beta / sig2e
  den   <- beta^2 / sig2e
  cn    <- cumsum(num[ord])
  cd    <- cumsum(den[ord])
  data.frame(
    stock     = stk[ord],
    Rbar_Rf   = (Rbar - Rf)[ord],
    beta      = beta[ord],
    ratio     = ratio[ord],
    sig2e     = sig2e[ord],
    num       = num[ord],
    cum_num   = cn,
    den       = den[ord],
    cum_den   = cd,
    Ci        = sig2m * cn / (1 + sig2m * cd),
    row.names = NULL
  )
}

tblA <- egp_table(Rf)
stopifnot(!is.unsorted(rev(tblA$ratio)))          # ranked high to low

cat("(a) excess return to beta ranking, Rf =", Rf, "\n")
print(rnd(head(tblA[, c("stock", "ratio", "beta", "sig2e", "Ci")], 5)))
cat("    ...\n")
cat(sprintf("    C_n (all 30 stocks) = %.6f\n\n", tblA$Ci[30]))

# ---- (b) tangency portfolio, single index, with and without short sales -----

egp_sim <- function(Rf, short) {
  t0    <- egp_table(Rf)
  if (short) {
    k <- nrow(t0)
  } else {
    ok <- which(t0$ratio > t0$Ci)
    if (!length(ok)) return(NULL)
    k <- max(ok)
  }
  Cstar <- t0$Ci[k]
  inc   <- t0$stock[seq_len(k)]
  z     <- (beta[inc] / sig2e[inc]) * ((Rbar[inc] - Rf) / beta[inc] - Cstar)
  list(x = z / sum(z), Cstar = Cstar, cutoff = k, table = t0)
}

ss   <- egp_sim(Rf, short = TRUE)
noss <- egp_sim(Rf, short = FALSE)

stopifnot(abs(sum(ss$x) - 1) < 1e-12, abs(sum(noss$x) - 1) < 1e-12)
stopifnot(all(noss$x > 0))                        # no short positions

# the short sales answer must reproduce Project 5 part (a), which solved
# Z = Sigma^-1 (Rbar - Rf 1) directly instead of ranking
Z5 <- solve(S_sim, Rbar - Rf)
x5 <- Z5 / sum(Z5)
stopifnot(max(abs(ss$x[stk] - x5[stk])) < 1e-9)

st_ss   <- stats(ss$x,   S_hist)
st_noss <- stats(noss$x, S_hist)
st_ss_sim   <- stats(ss$x,   S_sim)
st_noss_sim <- stats(noss$x, S_sim)

cat("(b) point of tangency under the single index model\n")
cat(sprintf("    short sales allowed : C* = %.6f, all %d stocks, %d short\n",
            ss$Cstar, length(ss$x), sum(ss$x < 0)))
cat(sprintf("                          E %.4f%%  sd %.4f%% (historical)  sd %.4f%% (single index)\n",
            st_ss["E"] * 100, st_ss["sd"] * 100, st_ss_sim["sd"] * 100))
cat(sprintf("    short sales not     : C* = %.6f, %d stocks kept of 30\n",
            noss$Cstar, noss$cutoff))
cat(sprintf("                          E %.4f%%  sd %.4f%% (historical)  sd %.4f%% (single index)\n",
            st_noss["E"] * 100, st_noss["sd"] * 100, st_noss_sim["sd"] * 100))
cat("    weights without short sales:\n")
print(round(noss$x, 4))
cat("\n")

# ---- (c) efficient frontier with no short sales, by varying Rf --------------

Rf_grid <- seq(-0.03, max(Rbar) - 1e-4, length.out = 400)
noss_fr <- lapply(Rf_grid, function(rf) {
  p <- egp_sim(rf, short = FALSE)
  if (is.null(p)) return(NULL)
  s <- stats(p$x, S_hist)
  c(Rf = rf, E = s["E"], sd = s["sd"], k = p$cutoff)
})
noss_fr <- as.data.frame(do.call(rbind, noss_fr[!vapply(noss_fr, is.null, logical(1))]))
names(noss_fr) <- c("Rf", "E", "sd", "k")

# every one of these portfolios is long only and fully invested, and none can
# sit inside the unconstrained frontier
stopifnot(all(noss_fr$sd >= sd_at(f_hist, noss_fr$E) - 1e-9))

cat("(c) no short sales frontier traced over", nrow(noss_fr), "values of Rf\n")
cat(sprintf("    Rf from %.4f to %.4f\n", min(noss_fr$Rf), max(noss_fr$Rf)))
cat(sprintf("    E  from %.4f%% to %.4f%%\n", min(noss_fr$E) * 100, max(noss_fr$E) * 100))
cat(sprintf("    stocks held: %d down to %d\n\n", max(noss_fr$k), min(noss_fr$k)))

# ---- (d) constant correlation model, handout #33 ---------------------------

CC   <- cor(ret30)
rho  <- mean(CC[upper.tri(CC)])

cc_table <- function(Rf) {
  ratio <- (Rbar - Rf) / sdv
  ord   <- order(ratio, decreasing = TRUE)
  i     <- seq_along(ord)
  cum   <- cumsum(ratio[ord])
  data.frame(
    stock     = stk[ord],
    Rbar_Rf   = (Rbar - Rf)[ord],
    sd        = sdv[ord],
    ratio     = ratio[ord],
    cum_ratio = cum,
    mult      = rho / (1 - rho + i * rho),
    Ci        = (rho / (1 - rho + i * rho)) * cum,
    row.names = NULL
  )
}

tblD <- cc_table(Rf)
stopifnot(!is.unsorted(rev(tblD$ratio)))

cat("(d) constant correlation model, Rf =", Rf, "\n")
cat(sprintf("    average pairwise correlation rho = %.6f\n", rho))
print(rnd(head(tblD[, c("stock", "ratio", "cum_ratio", "mult", "Ci")], 5)))
cat("    ...\n")
cat(sprintf("    C_n (all 30 stocks) = %.6f\n\n", tblD$Ci[30]))

# ---- (e) tangency portfolio, constant correlation --------------------------

egp_cc <- function(Rf, short) {
  t0 <- cc_table(Rf)
  if (short) {
    k <- nrow(t0)
  } else {
    ok <- which(t0$ratio > t0$Ci)
    if (!length(ok)) return(NULL)
    k <- max(ok)
  }
  Cstar <- t0$Ci[k]
  inc   <- t0$stock[seq_len(k)]
  z     <- (1 / ((1 - rho) * sdv[inc])) * ((Rbar[inc] - Rf) / sdv[inc] - Cstar)
  list(x = z / sum(z), Cstar = Cstar, cutoff = k, table = t0)
}

cc_ss   <- egp_cc(Rf, short = TRUE)
cc_noss <- egp_cc(Rf, short = FALSE)

stopifnot(abs(sum(cc_ss$x) - 1) < 1e-12, abs(sum(cc_noss$x) - 1) < 1e-12)
stopifnot(all(cc_noss$x > 0))

# the constant correlation matrix implied by rho, as a cross-check on the
# short sales solution: it must equal Sigma_cc^-1 (Rbar - Rf 1) rescaled
S_cc <- rho * outer(sdv, sdv)
diag(S_cc) <- sdv^2
Zcc <- solve(S_cc, Rbar - Rf)
stopifnot(max(abs(cc_ss$x[stk] - (Zcc / sum(Zcc))[stk])) < 1e-9)

st_cc_ss   <- stats(cc_ss$x,   S_hist)
st_cc_noss <- stats(cc_noss$x, S_hist)

cat("(e) point of tangency under the constant correlation model\n")
cat(sprintf("    short sales allowed : C* = %.6f, all %d stocks, %d short\n",
            cc_ss$Cstar, length(cc_ss$x), sum(cc_ss$x < 0)))
cat(sprintf("                          E %.4f%%  sd %.4f%%\n",
            st_cc_ss["E"] * 100, st_cc_ss["sd"] * 100))
cat(sprintf("    short sales not     : C* = %.6f, %d stocks kept of 30\n",
            cc_noss$Cstar, cc_noss$cutoff))
cat(sprintf("                          E %.4f%%  sd %.4f%%\n",
            st_cc_noss["E"] * 100, st_cc_noss["sd"] * 100))
cat("    weights without short sales:\n")
print(round(cc_noss$x, 4))
cat("\n")
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

# the common backdrop: 30 stocks, the index, and Project 2's frontier
backdrop <- function(xr, yr, xtick) {
  plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), v = xtick, col = grid, lwd = 1)
  Eg  <- seq(yr[1], yr[2], length.out = 800)
  sdg <- sd_at(f_hist, Eg)
  eff <- Eg >= f_hist$Emin
  lines(sdg[!eff], Eg[!eff], col = c_front, lwd = 1.1, lty = 3)
  lines(sdg[ eff], Eg[ eff], col = c_front, lwd = 2.2)
  points(sdv, Rbar, pch = 21, bg = "#d9d8d2", col = surface, lwd = 1, cex = 0.85)
  points(sd(Rm), mean(Rm), pch = 21, bg = c_spx, col = surface, lwd = 1.8, cex = 1.5)
  text(sd(Rm), mean(Rm), "S&P 500", pos = 1, offset = 0.6, cex = 0.72,
       font = 2, col = c_spx)
  ax(1, xtick, sprintf("%.0f%%", xtick * 100))
  ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
  mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
  mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
}

# ---- (b) single index tangency portfolios ----------------------------------

open_png("fig_sim_tangency.png", h = 5.0)
xr <- c(-0.014, max(c(sdv, st_ss["sd"])) * 1.05)
yr <- c(-0.004, max(c(Rbar, st_ss["E"])) * 1.12)
backdrop(xr, yr, pretty(c(0, max(sdv))))
pt(st_ss["sd"],   st_ss["E"],   c_eq,    "short sales allowed", 2)
pt(st_noss["sd"], st_noss["E"], c_front, "no short sales", 4)
mtext("(b) Point of tangency under the single index model", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.02, font = 2)
mtext("Blue curve is Project 2's frontier. Risk is measured with the historical covariance matrix.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

# ---- (c) the no short sales frontier ---------------------------------------

open_png("fig_noshort_frontier.png", h = 5.0)
o  <- order(noss_fr$E)
xr <- c(-0.010, max(sdv) * 1.05)
yr <- c(-0.004, max(c(Rbar, noss_fr$E)) * 1.10)
backdrop(xr, yr, pretty(c(0, max(sdv))))
lines(noss_fr$sd[o], noss_fr$E[o], col = c_eq, lwd = 2.8)
pt(st_noss["sd"], st_noss["E"], c_front, sprintf("Rf = %.3f", Rf), 4)
legend("bottomright", bty = "n", cex = 0.72, text.col = muted,
       lwd = c(2.8, 2.2), col = c(c_eq, c_front),
       legend = c("no short sales", "short sales allowed (Project 2)"))
mtext("(c) Efficient frontier with no short sales", 3, line = 1.8, adj = 0,
      col = ink, cex = 1.02, font = 2)
mtext(sprintf("Traced by varying Rf over %d values and re-optimising each time.",
              nrow(noss_fr)), 3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

# ---- (e) constant correlation tangency portfolios --------------------------

open_png("fig_cc_tangency.png", h = 5.0)
xr <- c(-0.017, max(c(sdv, st_cc_ss["sd"])) * 1.05)
yr <- c(-0.004, max(c(Rbar, st_cc_ss["E"])) * 1.10)
backdrop(xr, yr, pretty(c(0, max(c(sdv, st_cc_ss["sd"])))))
pt(st_cc_ss["sd"],   st_cc_ss["E"],   c_eq,    "short sales allowed", 2)
pt(st_cc_noss["sd"], st_cc_noss["E"], c_front, "no short sales", 4)
mtext("(e) Point of tangency under the constant correlation model", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.02, font = 2)
mtext(sprintf("Average pairwise correlation rho = %.4f.", rho),
      3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

cat("wrote fig_sim_tangency.png, fig_noshort_frontier.png, fig_cc_tangency.png\n")
