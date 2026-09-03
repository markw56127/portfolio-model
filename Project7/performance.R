# Project 7: the multigroup model, and out of sample evaluation of every
# portfolio built in Projects 1 to 7.

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "performance.R")
    setwd(dirname(normalizePath(f)))
})

DATA <- ".."

# ---- data -------------------------------------------------------------------
# Training period builds the portfolios, testing period evaluates them.

read_prices <- function(f) {
  a <- read.csv(file.path(DATA, f), sep = ",", header = TRUE, check.names = FALSE)
  stopifnot(!any(duplicated(names(a))))
  list(P = as.matrix(a[, -1][, -1]), d = as.Date(a$Date))
}

tr <- read_prices("stockData_train.csv")
te <- read_prices("stockData_test.csv")

rets <- function(P) { n <- nrow(P); P[-1, ] / P[-n, ] - 1 }
ret1 <- rets(tr$P)                      # 59 x 31
ret2 <- rets(te$P)                      # 54 x 31
dat2 <- te$d[-1]
stopifnot(nrow(ret1) == 59, nrow(ret2) == 54)

stk   <- setdiff(colnames(ret1), "^GSPC")
r1    <- ret1[, stk]; r2 <- ret2[, stk]
Rm1   <- as.numeric(ret1[, "^GSPC"]); Rm2 <- as.numeric(ret2[, "^GSPC"])
N     <- length(stk)
stopifnot(N == 30)

Rbar   <- colMeans(r1)
sdv    <- apply(r1, 2, sd)
S_hist <- cov(r1)
Rf     <- 0.001                         # same Rf as Projects 5 and 6

# single index inputs
fit   <- lapply(stk, function(s) lm(r1[, s] ~ Rm1))
names(fit) <- stk
beta  <- vapply(fit, function(f) unname(coef(f)[2]), numeric(1))
sig2e <- vapply(fit, function(f) summary(f)$sigma^2, numeric(1))
sig2m <- var(Rm1)
S_sim <- outer(beta, beta) * sig2m
diag(S_sim) <- beta^2 * sig2m + sig2e
dimnames(S_sim) <- list(stk, stk)

# Project 2's frontier, for the plot
one <- rep(1, N)
Si  <- solve(S_hist)
fA  <- as.numeric(t(one) %*% Si %*% Rbar); fB <- as.numeric(t(Rbar) %*% Si %*% Rbar)
fC  <- as.numeric(t(one) %*% Si %*% one);  fD <- fB * fC - fA^2
stopifnot(abs(fA - 36.9864) < 1e-3, abs(fC - 2330.03) < 1e-1)
sd_at <- function(E) sqrt(1 / fC + (fC / fD) * (E - fA / fC)^2)

# ---- (a) the multigroup model -----------------------------------------------
# Five sectors of six stocks each, in the order the data files are laid out.
# Correlation is assumed constant within a group and constant between any two
# groups, so the 30 x 30 correlation matrix collapses to a 5 x 5 one.

grp    <- rep(1:5, each = 6)
gnames <- c("Technology", "Communication", "Financial", "Consumer Cyclical",
            "Healthcare")
names(grp) <- stk
K <- 5; Nk <- as.numeric(table(grp))
stopifnot(all(Nk == 6))

CCm <- cor(r1)
rho <- matrix(0, K, K, dimnames = list(gnames, gnames))
for (k in 1:K) for (l in 1:K) {
  blk <- CCm[grp == k, grp == l]
  rho[k, l] <- if (k == l) mean(blk[upper.tri(blk)]) else mean(blk)
}
stopifnot(isSymmetric(rho), all(diag(rho) > 0))

# the multigroup covariance matrix
S_mg <- rho[grp, grp] * outer(sdv, sdv)
diag(S_mg) <- sdv^2
dimnames(S_mg) <- list(stk, stk)
stopifnot(isSymmetric(S_mg), min(eigen(S_mg, only.values = TRUE)$values) > 0)

# Solved directly, Z = Sigma^-1 (Rbar - Rf 1)
Z_mg <- solve(S_mg, Rbar - Rf)
x_mg <- Z_mg / sum(Z_mg)

# Handout #37 solves it through the group cut-offs instead. With
# A_k = sum of (Rbar_i - Rf)/sd_i over group k, the C_k satisfy
#   C_k + sum_l [rho_kl N_l/(1-rho_ll)] C_l = sum_l rho_kl A_l/(1-rho_ll)
# and then z_i = [(Rbar_i - Rf)/sd_i - C_k] / ((1-rho_kk) sd_i).
ratio <- (Rbar - Rf) / sdv
Ak    <- vapply(1:K, function(k) sum(ratio[grp == k]), numeric(1))
M     <- outer(rep(1, K), Nk / (1 - diag(rho))) * rho
bvec  <- as.numeric(rho %*% (Ak / (1 - diag(rho))))
Ck    <- solve(diag(K) + M, bvec)
names(Ck) <- gnames

z_h37 <- (ratio - Ck[grp]) / ((1 - diag(rho)[grp]) * sdv)
stopifnot(max(abs(z_h37 - Z_mg)) < 1e-9)          # the two routes agree
stopifnot(abs(sum(x_mg) - 1) < 1e-12)

cat("(a) multigroup model\n")
cat("    correlation within and between groups:\n")
print(round(rho, 4))
cat("\n    group cut-offs C_k:\n")
print(round(Ck, 6))
cat(sprintf("\n    max |handout 37 route - direct solve| = %.2e\n",
            max(abs(z_h37 - Z_mg))))
cat(sprintf("    %d short positions, sum|x| = %.2f\n\n", sum(x_mg < 0), sum(abs(x_mg))))

# ---- the portfolios from every project --------------------------------------

# Rebuilt exactly as in Project 6: rank, find the cut-off, weight.

egp_sim <- function(short) {
  ratio <- (Rbar - Rf) / beta
  ord   <- order(ratio, decreasing = TRUE)
  num   <- (Rbar - Rf) * beta / sig2e
  den   <- beta^2 / sig2e
  Ci    <- sig2m * cumsum(num[ord]) / (1 + sig2m * cumsum(den[ord]))
  k     <- if (short) N else max(which(ratio[ord] > Ci))
  inc   <- stk[ord][seq_len(k)]
  z     <- (beta[inc] / sig2e[inc]) * (ratio[inc] - Ci[k])
  z / sum(z)
}

rho_c <- mean(CCm[upper.tri(CCm)])

egp_cc <- function(short) {
  ratio <- (Rbar - Rf) / sdv
  ord   <- order(ratio, decreasing = TRUE)
  Ci    <- (rho_c / (1 - rho_c + seq_len(N) * rho_c)) * cumsum(ratio[ord])
  k     <- if (short) N else max(which(ratio[ord] > Ci))
  inc   <- stk[ord][seq_len(k)]
  z     <- (1 / ((1 - rho_c) * sdv[inc])) * (ratio[inc] - Ci[k])
  z / sum(z)
}

x_sim_ss   <- egp_sim(TRUE);  x_sim_noss <- egp_sim(FALSE)
x_cc_ss    <- egp_cc(TRUE);   x_cc_noss  <- egp_cc(FALSE)
stopifnot(all(x_sim_noss > 0), all(x_cc_noss > 0))

pad <- function(x) { w <- setNames(rep(0, N), stk); w[names(x)] <- x; w }

port <- list(
  "Equal allocation"      = setNames(rep(1 / N, N), stk),
  "Minimum risk"          = as.numeric(Si %*% one / fC),
  "Min risk (SIM)"        = as.numeric(solve(S_sim, one) / sum(solve(S_sim, one))),
  "Single index, SS"      = pad(x_sim_ss),
  "Single index, no SS"   = pad(x_sim_noss),
  "Const. corr., SS"      = pad(x_cc_ss),
  "Const. corr., no SS"   = pad(x_cc_noss),
  "Multigroup, SS"        = pad(x_mg)
)
port <- lapply(port, function(w) setNames(as.numeric(w), stk))
stopifnot(all(abs(vapply(port, sum, numeric(1)) - 1) < 1e-9))

# the short sales single index portfolio must still be Project 5 part (a)
Z5 <- solve(S_sim, Rbar - Rf)
stopifnot(max(abs(port[["Single index, SS"]] - Z5 / sum(Z5))) < 1e-9)

# training period statistics, for the plot
tr_E  <- vapply(port, function(w) sum(w * Rbar), numeric(1))
tr_sd <- vapply(port, function(w) sqrt(as.numeric(t(w) %*% S_hist %*% w)), numeric(1))

cat("(a) the multigroup portfolio on the training period\n")
cat(sprintf("    E = %.4f%%   sd = %.4f%%  (historical covariance matrix)\n\n",
            tr_E["Multigroup, SS"] * 100, tr_sd["Multigroup, SS"] * 100))

# ---- (b) out of sample evaluation, 01-Jan-2022 to 31-Jul-2026 ---------------

pr <- vapply(port, function(w) as.numeric(r2 %*% w), numeric(nrow(r2)))
pr <- cbind(pr, "S&P 500" = Rm2)
labels <- colnames(pr)

# b1: growth of one dollar
growth <- rbind(1, apply(1 + pr, 2, cumprod))

# b2: geometric mean
geo  <- apply(1 + pr, 2, function(z) prod(z)^(1 / length(z)) - 1)

# b3: Sharpe, differential excess return, Treynor, Jensen
mu   <- colMeans(pr)
sdp  <- apply(pr, 2, sd)
bet  <- apply(pr, 2, function(z) cov(z, Rm2) / var(Rm2))
mu_m <- mean(Rm2); sd_m <- sd(Rm2)

sharpe  <- (mu - Rf) / sdp
treynor <- (mu - Rf) / bet
jensen  <- mu - (Rf + bet * (mu_m - Rf))
diffex  <- mu - (Rf + (sdp / sd_m) * (mu_m - Rf))

stopifnot(abs(sharpe["S&P 500"] - (mu_m - Rf) / sd_m) < 1e-12)
stopifnot(abs(jensen["S&P 500"]) < 1e-12, abs(diffex["S&P 500"]) < 1e-12)

perf <- data.frame(portfolio = labels, geometric = geo, mean = mu, sd = sdp,
                   beta = bet, sharpe = sharpe, diff_excess = diffex,
                   treynor = treynor, jensen = jensen, row.names = NULL)

cat("(b) testing period", format(dat2[1]), "to", format(dat2[length(dat2)]),
    "-", nrow(pr), "months\n\n")
print(within(perf, {
  geometric <- round(geometric * 100, 4); mean <- round(mean * 100, 4)
  sd <- round(sd * 100, 4); beta <- round(beta, 4); sharpe <- round(sharpe, 4)
  diff_excess <- round(diff_excess * 100, 4); treynor <- round(treynor, 4)
  jensen <- round(jensen * 100, 4)
}))
cat("\n    (geometric, mean, sd, diff_excess and jensen in percent per month)\n\n")

# b4: Fama decomposition for the single index portfolio with no short sales
fp   <- "Single index, no SS"
f_mu <- mu[fp]; f_sd <- sdp[fp]; f_beta <- bet[fp]

fama <- c(
  overall      = unname(f_mu - Rf),
  risk         = unname(f_beta * (mu_m - Rf)),
  selectivity  = unname(f_mu - Rf - f_beta * (mu_m - Rf)),
  diversif     = unname((f_sd / sd_m - f_beta) * (mu_m - Rf)),
  net_select   = unname(f_mu - Rf - (f_sd / sd_m) * (mu_m - Rf))
)
stopifnot(abs(fama["overall"] - (fama["risk"] + fama["selectivity"])) < 1e-15)
stopifnot(abs(fama["selectivity"] - (fama["diversif"] + fama["net_select"])) < 1e-15)
stopifnot(abs(fama["selectivity"] - jensen[fp]) < 1e-15)

beta_eq <- unname(f_sd / sd_m)          # beta a fully diversified portfolio
                                        # with this total risk would carry

cat("(b4) Fama decomposition,", fp, "\n")
cat(sprintf("     beta %.4f, sd %.4f%%, beta equivalent of total risk %.4f\n",
            f_beta, f_sd * 100, beta_eq))
cat(sprintf("     overall performance   %8.4f%%\n", fama["overall"] * 100))
cat(sprintf("       from systematic risk%8.4f%%\n", fama["risk"] * 100))
cat(sprintf("       from selectivity    %8.4f%%\n", fama["selectivity"] * 100))
cat(sprintf("         diversification   %8.4f%%\n", fama["diversif"] * 100))
cat(sprintf("         net selectivity   %8.4f%%\n\n", fama["net_select"] * 100))
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

pal <- c("Equal allocation"    = "#6f6d67", "Minimum risk"        = "#2a78d6",
         "Min risk (SIM)"      = "#7fb2ea", "Single index, SS"    = "#eb6834",
         "Single index, no SS" = "#b8431a", "Const. corr., SS"    = "#8a5cd6",
         "Const. corr., no SS" = "#b79ae8", "Multigroup, SS"      = "#c9a227",
         "S&P 500"             = "#1baf7a")

# ---- (a) every portfolio on the training period axes ------------------------

open_png("fig_multigroup.png", h = 5.4)
xr <- c(-0.02, max(c(sdv, tr_sd)) * 1.05)
yr <- c(-0.004, max(c(Rbar, tr_E)) * 1.10)
xtick <- pretty(c(0, max(c(sdv, tr_sd))))
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = xtick, col = grid, lwd = 1)

Eg <- seq(yr[1], yr[2], length.out = 800); sdg <- sd_at(Eg)
eff <- Eg >= fA / fC
lines(sdg[!eff], Eg[!eff], col = c_front, lwd = 1.1, lty = 3)
lines(sdg[ eff], Eg[ eff], col = c_front, lwd = 2.2)

points(sdv, Rbar, pch = 21, bg = "#d9d8d2", col = surface, lwd = 1, cex = 0.8)
points(sd(Rm1), mean(Rm1), pch = 21, bg = c_spx, col = surface, lwd = 1.6, cex = 1.4)
text(sd(Rm1), mean(Rm1), "S&P 500", pos = 1, offset = 0.55, cex = 0.7,
     font = 2, col = c_spx)

for (nm in names(tr_E))
  points(tr_sd[nm], tr_E[nm], pch = 21, bg = pal[nm], col = surface,
         lwd = 1.6, cex = 1.5)
pt(tr_sd["Multigroup, SS"], tr_E["Multigroup, SS"], pal["Multigroup, SS"],
   "multigroup", 3, off = 0.7)

legend("bottomright", bty = "n", cex = 0.66, text.col = muted, ncol = 1,
       pch = 21, pt.bg = pal[names(tr_E)], col = surface, pt.lwd = 1.4,
       legend = names(tr_E))

ax(1, xtick, sprintf("%.0f%%", xtick * 100))
ax(2, pretty(yr), sprintf("%.0f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(a) The multigroup portfolio alongside every earlier portfolio", 3,
      line = 1.8, adj = 0, col = ink, cex = 1.02, font = 2)
mtext("Training period, historical covariance matrix. Blue curve is Project 2's frontier.",
      3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

# ---- (b1) growth of one dollar over the testing period ---------------------

open_png("fig_growth.png", h = 5.2)
gd <- c(dat2[1] - 31, dat2)
yr <- range(growth) * c(0.96, 1.06)
plot(NA, xlim = range(gd), ylim = yr, axes = FALSE, xlab = "", ylab = "")
yt <- pretty(yr)
abline(h = yt, col = grid, lwd = 1)
abline(h = 1, col = "#c3c2b7", lwd = 1.2, lty = 2)

for (nm in labels)
  lines(gd, growth[, nm], col = pal[nm], lwd = if (nm == "S&P 500") 2.8 else 1.7)

ax(2, yt, sprintf("%.1f", yt))
xt <- seq(as.Date("2022-01-01"), as.Date("2026-07-01"), by = "year")
axis(1, at = xt, labels = format(xt, "%Y"), col = "#c3c2b7", col.axis = muted,
     cex.axis = 0.75, lwd = 1)
legend("topleft", bty = "n", cex = 0.66, text.col = muted,
       lwd = 2, col = pal[labels], legend = labels)
mtext("Value of $1 invested", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(b1) Performance over the testing period", 3, line = 1.8, adj = 0,
      col = ink, cex = 1.02, font = 2)
mtext(sprintf("%s to %s, weights fixed at their training period values.",
              format(dat2[1], "%b %Y"), format(dat2[length(dat2)], "%b %Y")),
      3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

# ---- (b4) Fama decomposition on the return against beta plot ---------------

open_png("fig_fama.png", h = 5.2)
sml <- function(b) Rf + b * (mu_m - Rf)
xr  <- c(0, max(beta_eq, 1) * 1.22); yr <- c(0, max(f_mu, sml(xr[2])) * 1.12)
plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)

abline(a = Rf, b = mu_m - Rf, col = c_front, lwd = 2.2)
text(xr[2], sml(xr[2]), "SML", pos = 3, offset = 0.3, cex = 0.75, font = 2,
     col = c_front)

segments(f_beta, sml(f_beta), f_beta, f_mu, col = c_eq, lwd = 6)
segments(beta_eq, sml(f_beta), beta_eq, sml(beta_eq), col = "#c9a227", lwd = 6)
segments(beta_eq, sml(beta_eq), beta_eq, f_mu, col = "#8a5cd6", lwd = 6)
segments(f_beta, f_mu, beta_eq, f_mu, col = "#c9c8c0", lwd = 1, lty = 2)
segments(f_beta, sml(f_beta), beta_eq, sml(f_beta), col = "#c9c8c0", lwd = 1, lty = 2)

pt(f_beta, f_mu, "#b8431a", "portfolio", 2)
points(f_beta,  sml(f_beta),  pch = 21, bg = surface, col = c_front, lwd = 2, cex = 1.2)
points(beta_eq, sml(beta_eq), pch = 21, bg = surface, col = c_front, lwd = 2, cex = 1.2)
pt(1, mu_m, c_spx, "market", 4, off = 0.5)

text(f_beta,  (sml(f_beta) + f_mu) / 2, "selectivity", pos = 4, offset = 0.5,
     cex = 0.7, font = 2, col = c_eq)
text(beta_eq, (sml(f_beta) + sml(beta_eq)) / 2, "diversification", pos = 4,
     offset = 0.5, cex = 0.7, font = 2, col = "#c9a227")
text(beta_eq, (sml(beta_eq) + f_mu) / 2, "net selectivity", pos = 4, offset = 0.5,
     cex = 0.7, font = 2, col = "#8a5cd6")

ax(1, pretty(xr), sprintf("%.1f", pretty(xr)))
ax(2, pretty(yr), sprintf("%.1f%%", pretty(yr) * 100))
mtext(expression(paste("Beta   (", beta[p], " and ", sigma[p]/sigma[m], ")")),
      1, line = 2.6, col = muted, cex = 0.85)
mtext("Expected return (monthly)", 2, line = 3.2, col = muted, cex = 0.85)
mtext("(b4) Fama decomposition, single index with no short sales", 3, line = 1.8,
      adj = 0, col = ink, cex = 1.02, font = 2)
mtext(sprintf("Left tick is beta = %.3f; right tick is sd_p/sd_m = %.3f.",
              f_beta, beta_eq), 3, line = 0.5, adj = 0, col = muted, cex = 0.7)
dev.off()

cat("wrote fig_multigroup.png, fig_growth.png, fig_fama.png\n")
