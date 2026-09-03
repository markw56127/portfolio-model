# Homework 7
# Option payoff and profit diagrams. Base R only.
# Run with the Homework7 folder as the working directory.

# ---- Exercise 1(a): protective put = long call + bond ----------------------

ST <- seq(0, 100, by = 0.5)
E <- 50

stock  <- ST
put    <- pmax(E - ST, 0)
prot   <- stock + put              # protective put payoff
fiduc  <- pmax(ST - E, 0) + E      # long call plus a bond worth E at expiry

png("hw7_ex1a.png", width = 1300, height = 550, res = 150)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))

plot(ST, prot, type = "l", col = "black", lwd = 3, ylim = c(0, 100),
     xlab = "ST", ylab = "Payoff", main = "Protective put")
lines(ST, stock, col = "red", lty = 2)
lines(ST, put, col = "blue", lty = 2)
abline(h = 0, lty = 2, col = "gray")
legend("topleft", legend = c("long stock", "long put", "protective put"),
       col = c("red", "blue", "black"), lty = c(2, 2, 1), lwd = c(1, 1, 3),
       bty = "n", cex = 0.8)

plot(ST, fiduc, type = "l", col = "darkgreen", lwd = 4, ylim = c(0, 100),
     xlab = "ST", ylab = "Payoff", main = "Long call + bond")
lines(ST, prot, col = "black", lty = 3, lwd = 2)
lines(ST, pmax(ST - E, 0), col = "blue", lty = 2)
abline(h = E, col = "red", lty = 2)
legend("topleft", legend = c("long call", "bond paying E", "call + bond",
                             "protective put"),
       col = c("blue", "red", "darkgreen", "black"), lty = c(2, 2, 1, 3),
       lwd = c(1, 1, 4, 2), bty = "n", cex = 0.8)
dev.off()

stopifnot(max(abs(prot - fiduc)) < 1e-12)

# ---- Exercise 1(b): the target payoff from calls and the stock -------------

ST <- seq(0, 140, by = 0.5)

stock    <- ST
short50  <- -pmax(ST - 50, 0)
short60  <- -pmax(ST - 60, 0)
long110  <-  pmax(ST - 110, 0)
total    <- stock + short50 + short60 + long110

target <- ifelse(ST <= 50, ST, ifelse(ST <= 60, 50, pmax(110 - ST, 0)))
stopifnot(max(abs(total - target)) < 1e-12)

png("hw7_graph.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 2, 1))
plot(ST[ST <= 115], target[ST <= 115], type = "l", lwd = 2, ylim = c(0, 50),
     xlab = "ST", ylab = "Payoff", xaxt = "n")
axis(1, at = seq(0, 110, by = 10))
dev.off()

png("hw7_ex1b.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, total, type = "l", col = "black", lwd = 3, ylim = c(-90, 140),
     xlab = "ST", ylab = "Payoff",
     main = "Stock - call(50) - call(60) + call(110)")
lines(ST, stock,   col = "red",       lty = 2)
lines(ST, short50, col = "blue",      lty = 2)
lines(ST, short60, col = "darkgreen", lty = 2)
lines(ST, long110, col = "purple",    lty = 2)
abline(h = 0, lty = 2, col = "gray")
legend("bottomleft", legend = c("long stock", "short call E=50", "short call E=60",
                                "long call E=110", "total"),
       col = c("red", "blue", "darkgreen", "purple", "black"),
       lty = c(2, 2, 2, 2, 1), lwd = c(1, 1, 1, 1, 3), bty = "n", cex = 0.8)
dev.off()

# ---- Exercise 3: short call, E = 50, premium 4 -----------------------------

ST <- seq(30, 70, by = 0.5)
profit <- 4 - pmax(ST - 50, 0)

png("hw7_ex3.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, profit, type = "l", col = "blue", lwd = 2,
     xlab = "Stock Price at Expiration (ST)", ylab = "Profit ($)",
     main = "Short Call Profit (E = 50, C = 4)")
abline(h = 0, lty = 2, col = "gray")
abline(v = c(50, 54), lty = 3, col = "red")
points(54, 0, pch = 19, col = "red")
text(54, 0.8, "Breakeven = 54", col = "red", pos = 4)
dev.off()

# ---- Exercise 4: long put, E = 40, premium 3 -------------------------------

ST <- seq(20, 60, by = 0.5)
profit <- pmax(40 - ST, 0) - 3

png("hw7_ex4.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, profit, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Stock Price at Expiration (ST)", ylab = "Profit ($)",
     main = "Long Put Profit (E = 40, P = 3)")
abline(h = 0, lty = 2, col = "gray")
abline(v = c(40, 37), lty = 3, col = "red")
points(37, 0, pch = 19, col = "red")
text(37, 0.8, "Breakeven = 37", col = "red", pos = 2)
dev.off()

# ---- Exercise 5: two puts and one call, E = 50 -----------------------------

ST <- seq(20, 80, by = 0.5)
profit_puts  <- 2 * pmax(50 - ST, 0) - 12
profit_call  <- pmax(ST - 50, 0) - 5
profit_total <- profit_puts + profit_call

png("hw7_ex5.png", width = 1300, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))

plot(ST, profit_puts, type = "l", col = "red", lwd = 2,
     xlab = "ST", ylab = "Profit ($)", main = "(a) 2 Long Puts (E = 50)")
abline(h = 0, lty = 2, col = "gray"); abline(v = 44, lty = 3, col = "red")
text(44, min(profit_puts) + 5, "BE = 44", col = "red", pos = 4, cex = 0.8)

plot(ST, profit_call, type = "l", col = "blue", lwd = 2,
     xlab = "ST", ylab = "Profit ($)", main = "(b) 1 Long Call (E = 50)")
abline(h = 0, lty = 2, col = "gray"); abline(v = 55, lty = 3, col = "red")
text(55, max(profit_call) - 3, "BE = 55", col = "red", pos = 4, cex = 0.8)

plot(ST, profit_total, type = "l", col = "purple", lwd = 2, ylim = c(-20, 70),
     xlab = "ST", ylab = "Profit ($)", main = "(c) Strip: 2 Puts + 1 Call")
lines(ST, profit_puts, col = "red", lty = 2)
lines(ST, profit_call, col = "blue", lty = 2)
abline(h = 0, lty = 2, col = "gray")
abline(v = c(41.5, 67), lty = 3, col = "red")
legend("top", legend = c("Combined", "2 Puts", "1 Call"),
       col = c("purple", "red", "blue"), lty = c(1, 2, 2), lwd = 2,
       bty = "n", cex = 0.8)
dev.off()

# ---- Exercise 6: ratio call spread -----------------------------------------

ST <- seq(25, 65, by = 0.5)
profit <- pmax(ST - 40, 0) - 2 * pmax(ST - 45, 0) + 2

png("hw7_ex6.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, profit, type = "l", col = "blue", lwd = 2,
     xlab = "Stock Price at Expiration (ST)", ylab = "Profit ($)",
     main = "Buy 1 Call (E=40) + Sell 2 Calls (E=45)")
abline(h = 0, lty = 2, col = "gray")
abline(v = c(40, 45, 52), lty = 3, col = "red")
points(52, 0, pch = 19, col = "red")
text(52, 0.8, "Breakeven = 52", col = "red", pos = 4)
dev.off()

# ---- Exercise 7: put call parity, four rearrangements ----------------------
# Take r = 0 so that PV(E) = E and C = P = 5 satisfies P + S0 = C + E exactly.

ST <- seq(20, 80, by = 0.5)
S0 <- 50; E <- 50; C <- 5; P <- 5
stopifnot(abs((P + S0) - (C + E)) < 1e-12)

png("hw7_ex7.png", width = 1300, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))

# (a) long put + long stock = long call + bond
payoff_a <- ST + pmax(E - ST, 0)
profit_a <- payoff_a - (S0 + P)
plot(ST, payoff_a, type = "l", col = "blue", lwd = 3, ylim = c(-40, 100),
     main = "(a) Long Put + Long Stock", xlab = "ST", ylab = "Payoff / Profit")
lines(ST, pmax(ST - E, 0) + E, col = "darkgreen", lty = 3, lwd = 2)
lines(ST, profit_a, col = "red", lwd = 2)
abline(h = 0, lty = 2, col = "gray")
legend("topleft", legend = c("payoff", "= call + bond", "profit"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 3, 1), lwd = 2,
       bty = "n", cex = 0.75)

# (b) short put + short stock
payoff_b <- -payoff_a
profit_b <- payoff_b + (S0 + P)
plot(ST, payoff_b, type = "l", col = "blue", lwd = 3, ylim = c(-100, 40),
     main = "(b) Short Put + Short Stock", xlab = "ST", ylab = "Payoff / Profit")
lines(ST, -(pmax(ST - E, 0) + E), col = "darkgreen", lty = 3, lwd = 2)
lines(ST, profit_b, col = "red", lwd = 2)
abline(h = 0, lty = 2, col = "gray")
legend("bottomleft", legend = c("payoff", "= short call + borrow", "profit"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 3, 1), lwd = 2,
       bty = "n", cex = 0.75)

# (c) long call + short stock = long put + borrow PV(E)
payoff_c <- pmax(ST - E, 0) - ST
profit_c <- payoff_c + (S0 - C)
plot(ST, payoff_c, type = "l", col = "blue", lwd = 3, ylim = c(-80, 20),
     main = "(c) Long Call + Short Stock", xlab = "ST", ylab = "Payoff / Profit")
lines(ST, pmax(E - ST, 0) - E, col = "darkgreen", lty = 3, lwd = 2)
lines(ST, profit_c, col = "red", lwd = 2)
abline(h = 0, lty = 2, col = "gray")
legend("bottomleft", legend = c("payoff", "= put - bond", "profit"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 3, 1), lwd = 2,
       bty = "n", cex = 0.75)

# (d) short call + long stock = bond + short put (covered call)
payoff_d <- ST - pmax(ST - E, 0)
profit_d <- payoff_d - (S0 - C)
plot(ST, payoff_d, type = "l", col = "blue", lwd = 3, ylim = c(-20, 80),
     main = "(d) Short Call + Long Stock", xlab = "ST", ylab = "Payoff / Profit")
lines(ST, E - pmax(E - ST, 0), col = "darkgreen", lty = 3, lwd = 2)
lines(ST, profit_d, col = "red", lwd = 2)
abline(h = 0, lty = 2, col = "gray")
legend("topleft", legend = c("payoff", "= bond - put", "profit"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 3, 1), lwd = 2,
       bty = "n", cex = 0.75)
dev.off()

stopifnot(max(abs(payoff_a - (pmax(ST - E, 0) + E))) < 1e-12)
stopifnot(max(abs(payoff_c - (pmax(E - ST, 0) - E))) < 1e-12)
stopifnot(max(abs(payoff_d - (E - pmax(E - ST, 0)))) < 1e-12)

# ---- Exercise 8: box spread -------------------------------------------------

ST <- seq(30, 80, by = 0.5)
bull_call <- pmax(ST - 50, 0) - pmax(ST - 60, 0)
bear_put  <- pmax(60 - ST, 0) - pmax(50 - ST, 0)
total_box <- bull_call + bear_put
stopifnot(max(abs(total_box - 10)) < 1e-12)

png("hw7_ex8.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, total_box, type = "l", col = "purple", lwd = 3, ylim = c(-2, 14),
     xlab = "Stock Price at Expiration (ST)", ylab = "Payoff ($)",
     main = "Box Spread Payoff Diagram")
lines(ST, bull_call, col = "blue", lty = 2, lwd = 1.5)
lines(ST, bear_put, col = "red", lty = 2, lwd = 1.5)
abline(h = 0, lty = 2, col = "gray")
legend("topleft", legend = c("Box Spread (Total = 10)", "Bull Call Spread",
                             "Bear Put Spread"),
       col = c("purple", "blue", "red"), lty = c(1, 2, 2), lwd = 2, bty = "n")
dev.off()

# ---- Exercise 9: bear spread with puts -------------------------------------

ST <- seq(30, 70, by = 0.5)
E1 <- 45; E2 <- 55; P1 <- 2; P2 <- 6

payoff_buy   <- pmax(E2 - ST, 0)
payoff_sell  <- -pmax(E1 - ST, 0)
payoff_total <- payoff_buy + payoff_sell
profit_total <- payoff_total - (P2 - P1)

png("hw7_ex9.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, profit_total, type = "l", col = "darkgreen", lwd = 2.5, ylim = c(-12, 22),
     xlab = "Stock Price at Expiration (ST)", ylab = "Profit / Payoff ($)",
     main = "Bear Put Spread (E1 = 45, E2 = 55)")
lines(ST, payoff_total, col = "purple", lty = 2, lwd = 1.5)
lines(ST, payoff_buy - P2, col = "blue", lty = 3)
lines(ST, payoff_sell + P1, col = "red", lty = 3)
abline(h = 0, lty = 2, col = "gray")
abline(v = E2 - (P2 - P1), lty = 3, col = "gray40")
legend("topright", legend = c("Total Profit", "Total Payoff", "Long Put Profit",
                              "Short Put Profit"),
       col = c("darkgreen", "purple", "blue", "red"), lty = c(1, 2, 3, 3),
       lwd = 2, bty = "n", cex = 0.85)
dev.off()

# ---- Exercise 10: bear spread with calls -----------------------------------

ST <- seq(30, 70, by = 0.5)
E1 <- 45; E2 <- 55; C1 <- 7; C2 <- 2

payoff_sell  <- -pmax(ST - E1, 0)
payoff_buy   <- pmax(ST - E2, 0)
payoff_total <- payoff_sell + payoff_buy
profit_total <- payoff_total + (C1 - C2)

png("hw7_ex10.png", width = 900, height = 600, res = 150)
par(mar = c(4.2, 4.2, 3, 1))
plot(ST, profit_total, type = "l", col = "darkred", lwd = 2.5, ylim = c(-22, 12),
     xlab = "Stock Price at Expiration (ST)", ylab = "Profit / Payoff ($)",
     main = "Bear Call Spread (E1 = 45, E2 = 55)")
lines(ST, payoff_total, col = "purple", lty = 2, lwd = 1.5)
lines(ST, payoff_sell + C1, col = "red", lty = 3)
lines(ST, payoff_buy - C2, col = "blue", lty = 3)
abline(h = 0, lty = 2, col = "gray")
abline(v = E1 + (C1 - C2), lty = 3, col = "gray40")
legend("topright", legend = c("Total Profit", "Total Payoff", "Short Call Profit",
                              "Long Call Profit"),
       col = c("darkred", "purple", "red", "blue"), lty = c(1, 2, 3, 3),
       lwd = 2, bty = "n", cex = 0.85)
dev.off()

cat("all checks passed; wrote hw7_graph.png and hw7_ex*.png\n")
