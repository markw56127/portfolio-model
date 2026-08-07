# Portfolio project, parts (b)-(f)
# Monthly adjusted close prices, first 5 years only (2017-01 .. 2021-12).

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L && basename(f) == "portfolio.R")
    setwd(dirname(normalizePath(f)))
})


DATA <- ".."

# ---- (b) import prices and convert to returns -------------------------------

# check.names = FALSE keeps "^GSPC" intact instead of mangling it to "X.GSPC".
a <- read.csv(file.path(DATA, "stockData_train.csv"), sep = ",", header = TRUE,
              check.names = FALSE)

stopifnot(!any(duplicated(names(a))))

prices <- a[, -1]                     # drop the leading row-index column

dates <- as.Date(prices$Date)
P <- as.matrix(prices[, -1])          # 60 x 31, assets in columns

stopifnot(nrow(P) == 60)              # 5 years of monthly observations
stopifnot(ncol(P) == 31)              # ^GSPC + 30 stocks
stopifnot(!anyNA(P))

# Simple (arithmetic) returns: R_t = P_t / P_{t-1} - 1.
n <- nrow(P)
ret <- P[-1, ] / P[-n, ] - 1          # 59 x 31

cat("(b) prices ", nrow(P), "x", ncol(P),
    " -> returns", nrow(ret), "x", ncol(ret), "\n")
cat("    return period:", format(dates[2]), "to", format(dates[n]), "\n\n")

# ---- (c) means, standard deviations, variance-covariance matrix -------------

means  <- colMeans(ret)               # all 31 assets
covmat <- cov(ret)
stdev  <- sqrt(diag(covmat))

cat("(c) mean monthly return and sd, all 31 assets:\n")
print(round(data.frame(mean_pct = means * 100, sd_pct = stdev * 100), 3))
cat("\n    covariance matrix:", nrow(covmat), "x", ncol(covmat), "\n\n")

# The 30 stocks alone (drop ^GSPC) — this is what parts (e) and (f) use.
stk    <- setdiff(colnames(ret), "^GSPC")
ret30  <- ret[, stk]
mu30   <- colMeans(ret30)
S30    <- cov(ret30)                  # 30 x 30
stopifnot(ncol(S30) == 30)
stopifnot(qr(S30)$rank == 30)         # invertible, so part (f) is well posed

# ---- (e) equal allocation portfolio (30 stocks) -----------------------------

w_eq  <- rep(1 / 30, 30)
mu_eq <- as.numeric(t(w_eq) %*% mu30)
sd_eq <- sqrt(as.numeric(t(w_eq) %*% S30 %*% w_eq))

cat("(e) equal allocation portfolio (w = 1/30 each):\n")
cat(sprintf("    mean = %.4f (%.2f%% / month)   sd = %.4f (%.2f%%)\n\n",
            mu_eq, mu_eq * 100, sd_eq, sd_eq * 100))

# ---- (f) minimum risk portfolio ---------------------------------------------
# x = (S^-1 1) / (1' S^-1 1)

one    <- rep(1, 30)
Sinv   <- solve(S30)
x_min  <- Sinv %*% one / as.numeric(t(one) %*% Sinv %*% one)
mu_min <- as.numeric(t(x_min) %*% mu30)
sd_min <- sqrt(as.numeric(t(x_min) %*% S30 %*% x_min))

stopifnot(abs(sum(x_min) - 1) < 1e-10)   # weights sum to 1 by construction

cat("(f) minimum risk portfolio:\n")
cat(sprintf("    mean = %.4f (%.2f%% / month)   sd = %.4f (%.2f%%)\n",
            mu_min, mu_min * 100, sd_min, sd_min * 100))
cat(sprintf("    weights sum to %.6f, range %.3f to %.3f, %d short positions\n",
            sum(x_min), min(x_min), max(x_min), sum(x_min < 0)))
cat("    largest 5 weights:\n")
print(round(sort(x_min[, 1], decreasing = TRUE)[1:5], 4))
cat("\n")

# ---- (d) plot: expected return vs standard deviation ------------------------

surface <- "#fcfcfb"; ink <- "#0b0b0b"; muted <- "#898781"; grid <- "#e1e0d9"
c_gspc <- "#2a78d6"; c_eq <- "#eb6834"; c_min <- "#1baf7a"

repel <- function(x, y, lab, cex, fixed, iter = 600) {
  w <- strwidth(lab,  cex = cex) * 1.08
  h <- strheight(lab, cex = cex) * 1.85
  lx <- x; ly <- y + h * 0.85
  ly[fixed] <- y[fixed]
  for (it in seq_len(iter)) {
    moved <- FALSE
    for (i in seq_along(lx)) for (j in seq_along(lx)) if (i < j) {
      if (fixed[i] && fixed[j]) next
      dx <- lx[j] - lx[i]; dy <- ly[j] - ly[i]
      ox <- (w[i] + w[j]) / 2 - abs(dx)
      oy <- (h[i] + h[j]) / 2 - abs(dy)
      if (ox > 0 && oy > 0) {
        moved <- TRUE
        # split the correction, or give it all to the free one of the pair
        fi <- if (fixed[i]) 0 else if (fixed[j]) 1 else 0.5
        if (oy / h[i] < ox / w[i]) {
          s <- if (dy >= 0) oy else -oy
          ly[i] <- ly[i] - s * fi; ly[j] <- ly[j] + s * (1 - fi)
        } else {
          s <- if (dx >= 0) ox else -ox
          lx[i] <- lx[i] - s * fi; lx[j] <- lx[j] + s * (1 - fi)
        }
      }
    }
    if (!moved) break
  }
  list(x = lx, y = ly)
}

png("risk_return.png", width = 2100, height = 1450, res = 190)
op <- par(bg = surface, mar = c(4.6, 4.8, 4.2, 1.8), family = "sans")

# Pad left and right so the highlight labels have room to sit beside their marks
xr <- c(min(c(stdev, sd_min)) - 0.016, max(stdev) + 0.012)
yr <- range(c(means, mu_eq, mu_min)) + c(-0.006, 0.006)

plot(NA, xlim = xr, ylim = yr, xlab = "", ylab = "", axes = FALSE)
abline(h = pretty(yr), v = pretty(xr), col = grid, lwd = 1)
axis(1, col = "#c3c2b7", col.axis = muted, cex.axis = 0.85, lwd = 1,
     at = pretty(xr), labels = sprintf("%.0f%%", pretty(xr) * 100))
axis(2, col = "#c3c2b7", col.axis = muted, cex.axis = 0.85, lwd = 1, las = 1,
     at = pretty(yr), labels = sprintf("%.0f%%", pretty(yr) * 100))
mtext("Standard deviation (monthly)", 1, line = 2.9, col = muted, cex = 0.92)
mtext("Expected return (monthly)",    2, line = 3.4, col = muted, cex = 0.92)
mtext("Risk and return, 30 stocks vs. the S&P 500", 3, line = 2.4,
      adj = 0, col = ink, cex = 1.22, font = 2)
mtext("Monthly returns, Jan 2017 - Dec 2021 (60 months of prices, 59 returns)",
      3, line = 1.1, adj = 0, col = muted, cex = 0.9)

# The 30 individual stocks recede - they are context, not the finding.
is_g <- colnames(ret) == "^GSPC"
px <- c(stdev[!is_g], stdev[is_g], sd_eq, sd_min)
py <- c(means[!is_g], means[is_g], mu_eq, mu_min)
pl <- c(colnames(ret)[!is_g], "S&P 500", "Equal allocation (1/30)",
        "Minimum risk")
pf <- c(rep(FALSE, 30), TRUE, TRUE, TRUE)
pc <- c(rep(0.52, 30), 0.86, 0.86, 0.86)

points(px[1:30], py[1:30], pch = 21, bg = "#d9d8d2",
       col = surface, lwd = 1.4, cex = 1.25)

# Pin the highlight labels beside their marks, then let tickers repel around them
lp <- list(x = px, y = py)
lp$x[31] <- px[31]; lp$y[31] <- py[31] - 0.0052            # S&P 500, below
lp$x[32] <- px[32] - 0.0088; lp$y[32] <- py[32]            # equal alloc, left
lp$x[33] <- px[33] + 0.0072; lp$y[33] <- py[33]            # min risk, right
lp <- repel(px, py, pl, pc, pf)
lp$x[31:33] <- c(px[31], px[32] - 0.0088, px[33] + 0.0072)
lp$y[31:33] <- c(py[31] - 0.0052, py[32], py[33])

# Leader lines where a ticker drifted away from its point
d <- sqrt(((lp$x[1:30] - px[1:30]) / diff(xr))^2 +
          ((lp$y[1:30] - py[1:30]) / diff(yr))^2)
seg <- which(d > 0.018)
segments(px[seg], py[seg], lp$x[seg], lp$y[seg], col = "#dcdbd4", lwd = 0.8)
text(lp$x[1:30], lp$y[1:30], pl[1:30], cex = 0.52, col = muted)

# Three highlighted entities, each with a surface ring and a direct label.
for (k in 31:33) {
  col <- c(c_gspc, c_eq, c_min)[k - 30]
  points(px[k], py[k], pch = 21, bg = col, col = surface, lwd = 2.4, cex = 2.0)
  text(lp$x[k], lp$y[k], pl[k], cex = 0.86, font = 2, col = col,
       adj = c(if (k == 32) 1 else if (k == 33) 0 else 0.5, 0.5))
}

par(op); dev.off()
cat("(d) wrote risk_return.png\n")

# ---- save tables -------------------------------------------------------------

write.csv(data.frame(Date = dates[-1], ret, check.names = FALSE),
          "returns_train.csv", row.names = FALSE)
write.csv(data.frame(asset = names(means), mean = means, sd = stdev),
          "asset_stats.csv", row.names = FALSE)
write.csv(covmat, "covmat.csv")
cat("    wrote returns_train.csv, asset_stats.csv, covmat.csv\n")
